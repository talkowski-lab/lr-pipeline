#!/usr/bin/env python3
"""Breakpoint sharpness metrics for refined CNV calls (refine_cnv_breakpoints.py output), plus segmental-duplication
and short-read annotations, so sharpness can be compared between groups of calls.

For every call in --calls (a refine_cnv_breakpoints.py .calls.tsv.gz) whose SVID is in --bed (the BED the refinement
was run on, for the carrier list):
  profile  each carrier: CN = RD / median coverage x 2, rescaled so its median over its copy-neutral flanks
           ([start - 2 kb - L, start - 2 kb) and [end + 2 kb, end + 2 kb + L), L = call length) is 2; the call
           profile is the per-bin median across carriers
  at each breakpoint (refined coordinate, or original when not refined), within +/- --window:
    step    median CN 1-2 kb inside - median CN 1-2 kb outside
    J1, J2  largest jump towards the inside level over 1 / 2 bins, as a fraction of |step|
    W10_90  bp between the 10% and 90% points of the step (3-bin running median)
    edge_z  largest one-bin jump towards the inside level / 99th percentile of |bin-to-bin change| in the
            profile's own flanks
  call-level: min_edge_z, min_abs_step, min_J2 = the weaker of the two breakpoints
Annotations: sd_overlap = any --segdup interval overlapping [start - window, end + window] (original or refined,
whichever is wider); sr_site_match = any --sr-bed record of a matching type (DEL<->DEL/CNV, DUP<->DUP/CNV) with
reciprocal overlap >= --min-ro to the refined call; sr_match_ids.

Usage:
  compute_breakpoint_sharpness.py --calls chr6.gCNV_2k.per_sample_calls.bp_refined.calls.tsv.gz \\
      --bed chr6.gCNV_2k.merged_cnvs.per_sample_calls.bed.gz --rd hprc_hgsvc.chr6.RD.txt.gz \\
      --samples hprc_hgsvc.medianCov.transposed.bed --medians hprc_hgsvc.medianCov.transposed.bed \\
      --segdup hg38.SegDup.sorted.merged.bed.gz --sr-bed 1KGP_sr_hg38.chr6.DEL_DUP_CNV.bed.gz \\
      --out chr6.gCNV_2k.per_sample_calls.bp_sharpness.tsv.gz [--status PASS,PARTIAL,FAIL]
"""

import argparse
import bisect
import gzip
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

TYPE_MATCH = {"DEL": {"DEL", "CNV"}, "DUP": {"DUP", "CNV"}}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls", required=True, help="refine_cnv_breakpoints.py .calls.tsv.gz")
    p.add_argument("--bed", required=True, help="BED with SVID and carriers columns (input of the refinement)")
    p.add_argument("--rd", required=True)
    p.add_argument("--samples", help="sample IDs in RD column order (if --rd has no '#' header)")
    p.add_argument("--medians", required=True)
    p.add_argument("--segdup", required=True, help="segmental duplication BED (chrom, start, end)")
    p.add_argument("--sr-bed", help="short-read SV BED: #chrom start end ID SVTYPE FILTER (optional; sr_* columns stay empty)")
    p.add_argument("--min-ro", type=float, default=0.5)
    p.add_argument("--window", type=int, default=2000)
    p.add_argument("--status", default="PASS,PARTIAL,FAIL", help="statuses to include (default all)")
    p.add_argument("--out", required=True)
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    return p.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def read_ids(path):
    lines = [x.rstrip("\n") for x in open(path) if x.strip()]
    return [x for x in lines[0].split("\t") if x] if "\t" in lines[0] else [x.strip() for x in lines]


def load_intervals(path):
    iv = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or not f[1].isdigit():
                continue
            iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in iv.items():
        v.sort()
        out = [list(v[0])]
        for a, b in v[1:]:
            if a <= out[-1][1]:
                out[-1][1] = max(out[-1][1], b)
            else:
                out.append([a, b])
        merged[c] = [tuple(x) for x in out]
    return merged


def overlaps_any(iv, chrom, s, e):
    """True if [s, e) overlaps any interval (intervals merged first in load_intervals)."""
    v = iv.get(chrom, [])
    i = bisect.bisect_left(v, (e, -1)) - 1  # last interval starting before e
    return i >= 0 and v[i][1] > s


def sharpness(st, x, bp, side, bg, window):
    """Metrics at one breakpoint; side 'L' (outside = lower coordinates) or 'R'."""
    if side == "L":
        out = x[(st >= bp - window) & (st < bp - window // 2)]
        inn = x[(st >= bp + window // 2) & (st < bp + window)]
        m = (st >= bp - window) & (st < bp + window)
        seq, pos = x[m], st[m]
    else:
        out = x[(st >= bp + window // 2) & (st < bp + window)]
        inn = x[(st >= bp - window) & (st < bp - window // 2)]
        m = (st >= bp - window) & (st < bp + window)
        seq, pos = x[m][::-1], st[m][::-1]
    if len(out) < 3 or len(inn) < 3 or len(seq) < 5:
        return [np.nan] * 5
    o, i = np.median(out), np.median(inn)
    step = i - o
    if abs(step) < 1e-6:
        return [0.0, np.nan, 0.0, 0.0, 0.0]
    sgn = np.sign(step)
    sm = np.array([np.median(seq[max(0, k - 1) : k + 2]) for k in range(len(seq))])
    frac = (sm - o) / step
    i10 = int(np.argmax(frac >= 0.1))
    i90 = int(np.argmax(frac >= 0.9)) if (frac >= 0.9).any() else len(frac) - 1
    w = abs(int(pos[i90]) - int(pos[i10]))
    d1 = np.diff(seq) * sgn
    d2 = (seq[2:] - seq[:-2]) * sgn
    return [step, w, d1.max() / abs(step), d2.max() / abs(step), d1.max() / bg]


def main():
    args = parse_args()
    keep_status = set(args.status.split(","))
    header = subprocess.run(["tabix", "-H", args.rd], capture_output=True, text=True).stdout.splitlines()
    samples = header[-1].split("\t")[3:] if header and header[-1].startswith("#") else read_ids(args.samples)
    med_lines = [x.rstrip("\n").split("\t") for x in open(args.medians) if x.strip()]
    med = dict(zip([x for x in med_lines[0] if x], [float(x) for x in med_lines[1] if x.strip()]))
    scale = np.array([2.0 / med[s] for s in samples], dtype=np.float32)
    col = {s: i for i, s in enumerate(samples)}

    calls = pd.read_csv(args.calls, sep="\t")
    calls = calls[calls["status"].isin(keep_status)]
    bed = pd.read_csv(args.bed, sep="\t", keep_default_na=False, dtype={"carriers": str})
    carriers = dict(zip(bed["SVID"], bed["carriers"]))
    sd = load_intervals(args.segdup)
    sr = defaultdict(list)
    with open_text(args.sr_bed) if args.sr_bed else open(os.devnull) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, s, e, i, t = line.rstrip("\n").split("\t")[:5]
            sr[(c, t)].append((int(s), int(e), i))
    for k in sr:
        sr[k].sort()

    rows = []
    for chrom, sub in calls.groupby("chrom", sort=False):
        print(f"loading {chrom} RD ...", file=sys.stderr)
        proc = subprocess.Popen(["tabix", args.rd, chrom], stdout=subprocess.PIPE)
        rd = pd.read_csv(proc.stdout, sep="\t", header=None, dtype={0: str})
        proc.wait()
        starts = rd[1].to_numpy()
        cn_all = rd.iloc[:, 3:].to_numpy(dtype=np.float32) * scale
        del rd
        for r in sub.itertuples():
            s = int(r.new_start) if r.start_q > 0 else int(r.orig_start)
            e = int(r.new_end) if r.end_q > 0 else int(r.orig_end)
            L = e - s
            lo = np.searchsorted(starts, s - args.window - L)
            hi = np.searchsorted(starts, e + args.window + L)
            st = starts[lo:hi]
            fl = ((st >= s - args.window - L) & (st < s - args.window)) | ((st >= e + args.window) & (st < e + args.window + L))
            prof = []
            for c in str(carriers.get(r.SVID, "")).split(","):
                if c not in col:
                    continue
                x = cn_all[lo:hi, col[c]].astype(float)
                base = np.median(x[fl]) if fl.any() else np.nan
                if np.isfinite(base) and base >= 0.3:
                    prof.append(x / base * 2)
            if not prof or fl.sum() < 10:
                continue
            x = np.median(np.vstack(prof), axis=0)
            bg = max(np.percentile(np.abs(np.diff(x[fl])), 99), 0.02)
            left = sharpness(st, x, s, "L", bg, args.window)
            right = sharpness(st, x, e, "R", bg, args.window)
            ws, we = min(s, int(r.orig_start)) - args.window, max(e, int(r.orig_end)) + args.window
            sd_hit = overlaps_any(sd, chrom, ws, we)
            m = []
            for t in TYPE_MATCH.get(r.SVTYPE, ()):
                for ss, se, i in sr.get((chrom, t), []):
                    if ss >= e:
                        break
                    ov = min(e, se) - max(s, ss)
                    if ov > 0 and min(ov / L, ov / (se - ss)) >= args.min_ro:
                        m.append(i)
            rows.append(
                [chrom, s, e, r.SVID, r.SVTYPE, L, r.n_carriers_used, r.status, r.bp_score, r.n_sharp]
                + left
                + right
                + [
                    np.nanmin([left[4], right[4]]),
                    np.nanmin([abs(left[0]), abs(right[0])]),
                    np.nanmin([left[3], right[3]]),
                    sd_hit,
                    bool(m),
                    ",".join(m) or ".",
                ]
            )
        del cn_all
    cols = ["chrom", "start", "end", "SVID", "SVTYPE", "SVLEN", "n_carriers_used", "status", "bp_score", "n_sharp"]
    cols += [f"{s}_{k}" for s in ("left", "right") for k in ("step", "W10_90", "J1", "J2", "edge_z")]
    cols += ["min_edge_z", "min_abs_step", "min_J2", "sd_overlap", "sr_site_match", "sr_match_ids"]
    out = pd.DataFrame(rows, columns=cols)
    out.round(4).to_csv(args.out, sep="\t", index=False, compression="gzip")
    print(
        f"{len(out)} calls -> {args.out}; sd_overlap {int(out.sd_overlap.sum())}, sr_site_match {int(out.sr_site_match.sum())}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
