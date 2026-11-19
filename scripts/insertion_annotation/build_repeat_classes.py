#!/usr/bin/env python3
"""Combine the MEI stage, the STR/VNTR stage and RepeatMasker satellite calls into one repeat class per insertion.

Satellite type: RepeatMasker elements of class Satellite* (and HSATII-like (CATTC)n simple repeats) covering
>= --min-sat-frac of the insertion; the type is the satellite name with the most bp, grouped: ALR/Alpha, HSATI, HSATII,
CER, ACRO1, SAR (HSat1A), TAR1, TIFr, SATR, LSAU/BSAT (LSAU and COMP-subunit_LSAU-BSAT families), HSAT5, BSR/Beta,
GSATII, other satellite.
repeat_class (first match wins):
  MEI full-length <family> / MEI truncated <family>   VCF MEIs (family from the VCF) and new MEI insertions (family from
                                                      RepeatMasker), split by completeness
  VNTR in MEI                                          VNTR variant within SVA, TR variant within L1
  segmental duplication paralog                        interspersed STR / VNTR (classify_tr_insertions.py) whose insertion site
                                                      and best-alignment source both lie in --sd-blocks (e.g. "LCR22A inserted
                                                      into LCR22D", or "LCR22D copy within LCR22D")
  MEI-family sequence without insertion hallmarks (annotate_mei_structure.py), by its local source
  (locate_local_source.py, --te-local):
    tandem duplication (in <family>)                  tandem / inverted duplication of the adjacent reference sequence
    interspersed duplication (in <family>)            nearby duplication whose source interval is >= --min-src-te of a
                                                      reference <family> element and covers >= --min-src-cov of the insertion
                                                      (with --te-local-aln: the gapped alignment's aligned fraction is used,
                                                      and its strand and identity are reported)
    <family> insertion into <site family>             no real local source ("partly local"); the site is inside a reference
                                                      element (e.g. Alu into Alu); "<family> insertion, no local source"
                                                      when the site is not inside one
    non-repeat (next stage)                           no local source, or a nearby duplication that fails the rule above
  satellite: <type>
  STR                                                  STR expansion (incl. diverged units), interspersed STR, local STR copy
  VNTR                                                 VNTR expansion (incl. diverged units), interspersed VNTR, local VNTR copy
  STR_novel / VNTR_novel                               uncatalogued tandem-repeat sequence ("_novel" = not in the current
                                                       TR catalog), typed by its period; local-source evidence in the detail
  MEI nested Alu                                       Alu chains / nested Alus (detail = RepeatMasker chain with orientation)
  other TR                                             in TR partly matching, in TR different tandem repeat
  non-repeat (next stage)                              not TR, in TR non-repeat sequence
"""
import argparse
import collections
import csv
import gzip
import re


def sat_group(name):
    if name.startswith("ALR"):
        return "ALR/Alpha"
    if name.startswith("HSATII") or re.match(r"^\((CATTC|GAATG|ATTCC|GGAAT|AATGG|CCATT|TTCCA|TGGAA)\)n$", name):
        return "HSATII"
    if name.startswith("HSATI"):
        return "HSATI"
    if name.startswith("HSAT5"):
        return "HSAT5"
    if name.startswith(("LSAU", "COMP-subunit_LSAU")):
        return "LSAU/BSAT"
    if name.startswith("SATR"):
        return "SATR"
    for k in ("CER", "ACRO1", "SAR", "TAR1", "TIFr", "BSR/Beta", "GSATII"):
        if name.startswith(k):
            return k
    return "other satellite"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mei-stage", required=True, help="annotate_mei_structure.py output")
    ap.add_argument("--tr-stage", required=True, help="classify_tr_insertions.py output")
    ap.add_argument("--rm-out", required=True, help="RepeatMasker .out of the inserted sequences")
    ap.add_argument("--sd-blocks", default=None,
                    help="Optional BED of named segmental-duplication blocks (e.g. LCR22A-H) for paralog insertions")
    ap.add_argument("--te-local", default=None,
                    help="Optional locate_local_source.py output for MEI-family sequence without hallmarks")
    ap.add_argument("--repeats", default=None, help="RepeatMasker-style reference BED.gz (needed with --te-local)")
    ap.add_argument("--te-local-aln", default=None,
                    help="Optional align_to_local_origin.py output for the --te-local IDs (aligned fraction, identity)")
    ap.add_argument("--min-src-te", type=float, default=0.9, help="Source fraction inside a same-family reference element [0.9]")
    ap.add_argument("--min-src-cov", type=float, default=0.9, help="Insertion fraction covered by the nearby source [0.9]")
    ap.add_argument("--tr-novel-local", default=None,
                    help="Optional locate_local_source.py output for the uncatalogued tandem repeats (evidence in the detail)")
    ap.add_argument("--min-sat-frac", type=float, default=0.5, help="Satellite fraction of the insertion [0.5]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    sat = collections.defaultdict(collections.Counter)
    with open(a.rm_out) as f:
        for i, line in enumerate(f):
            c = line.split()
            if i < 3 or len(c) < 15:
                continue
            name, cls = c[9], c[10]
            is_sat = cls.startswith("Satellite") or sat_group(name) == "HSATII"
            if is_sat:
                sat[c[4]][sat_group(name)] += int(c[6]) - int(c[5]) + 1
    blocks = []
    if a.sd_blocks:
        with open(a.sd_blocks) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                blocks.append((c[0], int(c[1]), int(c[2]), c[3]))

    def block_of(chrom, pos):
        for ch, s0, e0, name in blocks:
            if ch == chrom and s0 < pos <= e0:
                return name
        return ""

    loc = {}
    ref_el = collections.defaultdict(list)
    aln = {r["ID"]: r for r in csv.DictReader(open(a.te_local_aln), delimiter="\t")} if a.te_local_aln else {}
    if a.te_local:
        loc = {r["ID"]: r for r in csv.DictReader(open(a.te_local), delimiter="\t")}
        with gzip.open(a.repeats, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if c[0] == "chr22" or c[0] in {r["local_source"].split(":")[0] for r in loc.values() if r["local_source"]}:
                    fam = {"ERVK": "HERVK"}.get(c[5], c[5])
                    ref_el[(c[0], fam)].append((int(c[1]), int(c[2])))

    def te_frac_of(interval, fam):
        ch, rng = interval.split(":")
        s0, e0 = (int(v) for v in rng.split("-"))
        cov = 0
        for rs, re_ in ref_el.get((ch, fam), []):
            cov += max(0, min(e0, re_) - max(s0 - 1, rs))
        return min(1.0, cov / max(1, e0 - s0 + 1))

    novel_loc = {r["ID"]: r for r in csv.DictReader(open(a.tr_novel_local), delimiter="\t")} if a.tr_novel_local else {}
    mei = {r["ID"]: r for r in csv.DictReader(open(a.mei_stage), delimiter="\t")}
    tr = {r["ID"]: r for r in csv.DictReader(open(a.tr_stage), delimiter="\t")}

    cols = ["ID", "len", "repeat_class", "repeat_detail", "sat_type", "sat_frac", "mei_stage", "tr_stage",
            "completeness",
            "te_subfamily", "kept_parts", "rm_chain"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for iid, m in mei.items():
            t = tr.get(iid, {})
            n = int(m["len"])
            st, ms = t.get("tr_stage", ""), m["mei_stage"]
            sc = sat.get(iid, collections.Counter())
            stype, sbp = (sc.most_common(1)[0] if sc else ("", 0))
            sfrac = sum(sc.values()) / n if n else 0.0
            detail = ""
            te_note = (f"{m['te_subfamily']} {m['completeness']}, strand {m['te_strand']}, "
                       f"identity {m['te_identity']} to consensus") if m["te_subfamily"] else ""
            if ms.startswith("VCF MEI") or ms.startswith("new "):
                fam = ms.split()[2] if ms.startswith("VCF MEI") else ms.split()[1]
                kind = "full-length" if m["completeness"] == "full-length" else "truncated"
                cls = f"MEI {kind} {fam}"
                detail = f"{m['te_subfamily']} {m['completeness']} [{m['kept_parts']}], strand {m['te_strand']}, "
                detail += f"identity {m['te_identity']} to consensus"
                if not ms.startswith("VCF"):
                    detail += " (not in VCF)"
            elif ms in ("VNTR variant within SVA", "TR variant within L1"):
                cls, detail = "VNTR in MEI", ms
            elif blocks and st.startswith("interspersed") and t.get("best_hit"):
                bc, br_ = t["best_hit"].split(":")
                src = block_of(bc, (int(br_.split("-")[0]) + int(br_.split("-")[1])) // 2)
                site = block_of(iid.split("-")[0], int(iid.split("-")[1]))
                if src and site:
                    cls = "segmental duplication paralog"
                    detail = f"{src} inserted into {site}" if src != site else f"{src} copy within {site}"
                    detail += f" ({st}, source {t['best_hit']})"
                elif sfrac >= a.min_sat_frac:
                    cls, detail = f"satellite: {stype}", f"{stype} {sfrac:.2f} of insertion"
                else:
                    cls, detail = ("STR" if "STR" in st else "VNTR"), st
            elif ms.endswith("without insertion hallmarks") and iid in loc:
                fam = ms.split()[0]
                lc = loc[iid]
                note = te_note
                al = aln.get(iid)
                cov = float(al["aligned_frac"]) if al and al["aligned_frac"] else float(lc["local_cov"])
                ev = (f", identity {al['identity']} over {100 * float(al['aligned_frac']):.0f}% of the insertion"
                      if al and al["identity"] else f", k-mer coverage {lc['local_cov']}")
                site_fams = [x for x in m.get("site_in_ref", "").split(",") if x]
                if lc["local_call"] in ("tandem duplication", "inverted duplication"):
                    cls = f"tandem duplication (in {fam})"
                    detail = f"{lc['local_call']} of {lc['local_source']} ({lc['local_strand']}{ev}); {note}"
                elif (lc["local_call"] == "nearby duplication" and cov >= a.min_src_cov
                      and te_frac_of(lc["local_source"], fam) >= a.min_src_te):
                    cls = f"interspersed duplication (in {fam})"
                    detail = f"source {lc['local_source']} ({lc['local_strand']}{ev}), {lc['local_dist']} bp away; {note}"
                elif lc["local_call"] == "partly local":
                    tgt = fam if fam in site_fams else (site_fams[0] if site_fams else "")
                    cls = f"{fam} insertion into {tgt}" if tgt else f"{fam} insertion, no local source"
                    detail = f"{note}; TSD {m.get('tsd_bp', '')} bp, poly-A {m.get('polyA', '')}"
                else:
                    cls, detail = "non-repeat (next stage)", f"{ms}; {lc['local_call']}"
            elif sfrac >= a.min_sat_frac:
                cls, detail = f"satellite: {stype}", f"{stype} {sfrac:.2f} of insertion"
            elif st.startswith("STR ") or st.startswith("interspersed STR"):
                cls, detail = "STR", st
            elif st.startswith("VNTR ") or st.startswith("interspersed VNTR"):
                cls, detail = "VNTR", st
            elif ms == "Alu chain / nested Alus" or ms.endswith("without insertion hallmarks"):
                cls, detail = "MEI nested Alu", m["rm_chain"]
            elif st == "tandem-repeat sequence, uncatalogued":
                cls = (t.get("tr_type") or "VNTR") + "_novel"
                ev = novel_loc.get(iid, {}).get("local_call", "")
                detail = f"not in the TR catalog; period {t.get('period', '')} bp" + (f"; {ev}" if ev else "")
            elif st in ("in TR, partly matching", "in TR, different tandem repeat"):
                cls, detail = "other TR", st
            else:
                cls, detail = "non-repeat (next stage)", st
            row = dict(ID=iid, len=n, repeat_class=cls, repeat_detail=detail, sat_type=stype,
                       sat_frac=round(sfrac, 3),
                       mei_stage=ms, tr_stage=st, completeness=m["completeness"], te_subfamily=m["te_subfamily"],
                       kept_parts=m["kept_parts"], rm_chain=m["rm_chain"])
            o.write("\t".join(str(row[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
