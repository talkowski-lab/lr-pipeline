#!/usr/bin/env python3
"""Apply the depth-based breakpoint refinement and filter to per-sample gCNV calls and write the passing calls as a
VCF with one record per (call, carrier) at the refined breakpoints.

Decision per sample-call (rows of split_depth_vcf_calls.py --out-bed):
  SD-involved calls (in --sd-members, from refine_sd_cnv_breakpoints.py): kept when the call's SD event is SD_PASS;
      breakpoints = the event's L_bp / R_bp
  other calls (in --scored, from score_sharpness_gbm.py): kept when p_gbm >= --gbm-cutoff; breakpoints = the refined
      start / end of refine_cnv_breakpoints.py
  calls in neither table (no usable depth): removed
  a kept call whose refined interval overlaps --gaps is removed
Output VCF: the header of --vcf; for each kept sample-call, a copy of its source record with ID <ID>__<sample>,
POS / END at the refined breakpoints, the carrier's FORMAT values unchanged and every other sample set to GT 0/0 with
CN = ECN. Records are sorted and the VCF is bgzipped and indexed.
Output table: one row per sample-call with class, decision, original and refined coordinates, p_gbm and SD status.

Usage:
  filter_depth_cnv_calls.py --vcf merged.vcf.gz --bed calls.per_sample.bed.gz --out-prefix out \\
      [--scored calls.scored.tsv.gz --sd-members sd.members.tsv.gz --sd-events sd.events.tsv.gz] [--gaps gaps.bed]
      [--gbm-cutoff 0.5]
"""

import argparse
import bisect
import gzip
import sys
from collections import defaultdict

import pandas as pd
import pysam

COLUMNS = [
    "chrom",
    "start",
    "end",
    "SVID",
    "orig_id",
    "sample",
    "SVTYPE",
    "SVLEN",
    "cls",
    "kept",
    "new_start",
    "new_end",
    "p_gbm",
    "sd_status",
]


def load_intervals(path):
    iv = defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or line.startswith(("#", "track", "browser")) or not f[1].isdigit():
                continue
            iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in iv.items():
        v.sort()
        out = []
        for s, e in v:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[c] = out
    return merged


def overlaps(iv, c, s, e):
    v = iv.get(c, [])
    i = bisect.bisect_right(v, (e, e)) - 1
    while i >= 0 and v[i][1] > s:
        if v[i][0] < e:
            return True
        i -= 1
    return False


def decide(bed, scored, members, events, gaps, cutoff):
    rows = []
    for r in bed.itertuples(index=False):
        p_gbm, sd_status = float("nan"), "."
        if r.SVID in members.index:
            ev = events.loc[members.at[r.SVID, "event"]]
            cls = "SD:" + members.at[r.SVID, "group"]
            sd_status = ev["status"]
            keep = sd_status == "SD_PASS"
            ns = int(ev["L_bp"]) if pd.notna(ev["L_bp"]) else r.start
            ne = int(ev["R_bp"]) if pd.notna(ev["R_bp"]) else r.end
            if ns >= ne:
                ns, ne = r.start, r.end
        elif r.SVID in scored.index:
            x = scored.loc[r.SVID]
            cls = "SD nearby" if bool(x["sd_overlap"]) else "no SD"
            p_gbm = float(x["p_gbm"])
            keep = p_gbm >= cutoff
            ns, ne = int(x["start"]), int(x["end"])
        else:
            cls, keep, ns, ne = "not scored", False, r.start, r.end
        if keep and overlaps(gaps, r.chrom, ns, ne):
            cls, keep = cls + " (refined into gap)", False
        row = [r.chrom, r.start, r.end, r.SVID, r.orig_id, r.carriers, r.SVTYPE, r.end - r.start]
        rows.append(row + [cls, keep, ns, ne, p_gbm, sd_status])
    return pd.DataFrame(rows, columns=COLUMNS)


def write_vcf(vcf_path, dec, out_vcf):
    keep = dec[dec["kept"]]
    want = defaultdict(list)
    for r in keep.itertuples(index=False):
        want[r.orig_id].append(r)
    with pysam.VariantFile(vcf_path) as vin:
        recs = []
        for rec in vin:
            for r in want.get(rec.id, []):
                new = rec.copy()
                new.id = r.SVID
                new.start = int(r.new_start)
                new.stop = int(r.new_end)
                for s in new.samples:
                    if s == r.sample:
                        continue
                    smp = new.samples[s]
                    smp["GT"] = (0, 0)
                    if "CN" in smp and "ECN" in smp and smp["ECN"] is not None:
                        smp["CN"] = smp["ECN"]
                recs.append(new)
        recs.sort(key=lambda x: (x.start, x.stop, x.id))
        with pysam.VariantFile(out_vcf, "wz", header=vin.header) as vout:
            for rec in recs:
                vout.write(rec)
    pysam.tabix_index(out_vcf, preset="vcf", force=True)
    return len(recs)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="merged gCNV depth VCF the BED was made from")
    p.add_argument("--bed", required=True, help="split_depth_vcf_calls.py output")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--scored", help="score_sharpness_gbm.py output")
    p.add_argument("--sd-members", help="refine_sd_cnv_breakpoints.py .members.tsv.gz")
    p.add_argument("--sd-events", help="refine_sd_cnv_breakpoints.py .events.tsv.gz")
    p.add_argument("--gaps", help="BED of assembly gaps, centromeres and telomeres")
    p.add_argument("--gbm-cutoff", type=float, default=0.5, help="minimum p_gbm to keep a call without SD (default 0.5)")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    bed = pd.read_csv(a.bed, sep="\t", dtype={"carriers": str, "orig_id": str}, keep_default_na=False)
    bed = bed.rename(columns={bed.columns[0]: "chrom"})
    empty = pd.DataFrame(columns=["SVID"]).set_index("SVID")
    scored = pd.read_csv(a.scored, sep="\t").set_index("SVID") if a.scored else empty
    members = pd.read_csv(a.sd_members, sep="\t").set_index("SVID") if a.sd_members else empty
    events = pd.read_csv(a.sd_events, sep="\t").set_index("event") if a.sd_events else empty
    gaps = load_intervals(a.gaps) if a.gaps else {}
    dec = decide(bed, scored, members, events, gaps, a.gbm_cutoff)
    dec.to_csv(a.out_prefix + ".decisions.tsv.gz", sep="\t", index=False, compression="gzip")
    n = write_vcf(a.vcf, dec, a.out_prefix + ".vcf.gz")
    summary = dec.groupby("cls")["kept"].agg(["size", "sum"]).to_dict("index") if len(dec) else {}
    print(f"{len(dec)} sample-calls, {n} kept; by class (n, kept): {summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
