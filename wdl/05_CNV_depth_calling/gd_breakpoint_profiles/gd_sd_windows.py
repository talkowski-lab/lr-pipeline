#!/usr/bin/env python3
"""For each genomic-disorder (GD) locus in a list of GD call IDs, find the segmental duplication (SD) at its left and right
breakpoints and the plotting window around each (SD +/- --flank-frac x SD length), one row per locus (DEL and DUP of the
same locus and coordinates share a row).

GD call IDs look like GD_<locus>_<DEL|DUP>_<chrom>_<start>_<end>[_<n>] (one per line; a header line such as GD_CALL_IDS
and any leading count column, e.g. from `uniq -c`, are ignored).

Breakpoint SDs, in order of preference:
  1. --bp-sd (GD_NAHR.breakpoint_SDs.bed.gz): left = the row whose end equals the GD start (or lists it in
     right_GD_candidates); right = the row whose start equals the GD end (or lists it in left_GD_candidates)
  2. --segdup (merged hg38 SD BED): the merged SD block touching the GD edge (within --touch bp), extended outward
     through neighbouring blocks separated by <= --max-gap bp
A side with no SD gets a window of --fallback-flank bp outside the GD edge, flagged in sd_source.

Output TSV columns: locus, chrom, gd_start, gd_end, types, n_calls, highlight_pattern,
  left_sd_name, left_sd_start, left_sd_end, left_sd_source, left_win_start, left_win_end,
  right_sd_name, right_sd_start, right_sd_end, right_sd_source, right_win_start, right_win_end
(windows 0-based half-open, clipped at 0; highlight_pattern is a regex matching the renamed sample IDs of the locus).

Usage:
  gd_sd_windows.py --calls gd_calls.txt --bp-sd references/GD_NAHR.breakpoint_SDs.bed.gz \\
      --segdup references/hg38.SegDup.sorted.merged.bed.gz --out gd_loci.windows.tsv [--flank-frac 0.5]
"""
import argparse
import re
import sys
from collections import OrderedDict

import pandas as pd

CALL = re.compile(r"(GD_(.+)_(DEL|DUP)_(chr[0-9XYM]+)_(\d+)_(\d+))(?:_\d+)?$")


def listed(v, x):
    return pd.notna(v) and str(x) in str(v).split(",")


def merged_side(sd, chrom, x, side, touch, max_gap):
    """merged SD block touching x, chained outward (side 'left' grows to lower coordinates)"""
    d = sd[sd.c == chrom].sort_values("s").reset_index(drop=True)
    hit = d[(d.s <= x + touch) & (d.e >= x - touch)]
    if hit.empty:
        return None
    if side == "left":
        i = hit.index.min()
        s, e = int(d.at[i, "s"]), int(min(d.at[i, "e"], x))
        while i > 0 and s - d.at[i - 1, "e"] <= max_gap:
            i -= 1
            s = int(min(s, d.at[i, "s"]))
    else:
        i = hit.index.max()
        s, e = int(max(d.at[i, "s"], x)), int(d.at[i, "e"])
        while i + 1 < len(d) and d.at[i + 1, "s"] - e <= max_gap:
            i += 1
            e = int(max(e, d.at[i, "e"]))
    return s, e


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calls", required=True, help="GD call IDs, one per line")
    p.add_argument("--bp-sd", required=True, help="GD NAHR breakpoint-SD table")
    p.add_argument("--segdup", required=True, help="merged SD BED")
    p.add_argument("--out", required=True)
    p.add_argument("--flank-frac", type=float, default=0.5)
    p.add_argument("--touch", type=int, default=1000, help="max distance of a merged SD block from the GD edge")
    p.add_argument("--max-gap", type=int, default=5000, help="max gap between chained merged SD blocks")
    p.add_argument("--fallback-flank", type=int, default=100000, help="window width when a side has no SD")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()

    loci = OrderedDict()
    for line in open(a.calls):
        f = line.split()
        if not f:
            continue
        m = CALL.match(f[-1])
        if not m:
            continue
        _, locus, t, chrom, s, e = m.groups()
        key = (locus, chrom, int(s), int(e))
        loci.setdefault(key, {"types": set(), "n": 0})
        loci[key]["types"].add(t)
        loci[key]["n"] += 1
    if not loci:
        sys.exit("no GD call IDs found")

    bp = pd.read_csv(a.bp_sd, sep="\t")
    bp = bp.rename(columns={bp.columns[0]: "chrom"})
    sd = pd.read_csv(a.segdup, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])

    rows = []
    for (locus, chrom, gs, ge), v in loci.items():
        d = bp[bp.chrom == chrom]
        L = d[(d.end == gs) | d.right_GD_candidates.apply(lambda x: listed(x, gs))]
        R = d[(d.start == ge) | d.left_GD_candidates.apply(lambda x: listed(x, ge))]
        row = dict(locus=locus, chrom=chrom, gd_start=gs, gd_end=ge, types="/".join(sorted(v["types"])), n_calls=v["n"],
                   highlight_pattern=f"GD_{re.escape(locus)}_(DEL|DUP)_{chrom}_{gs}_{ge}")
        for side, hit, edge in (("left", L, gs), ("right", R, ge)):
            if not hit.empty:
                h = hit.iloc[0]
                s, e, name, src = int(h.start), int(h.end), f"{h.cluster}:{h.breakpoint}", "breakpoint_SD"
            else:
                m = merged_side(sd, chrom, edge, side, a.touch, a.max_gap)
                if m:
                    (s, e), name, src = m, f"{side}_SD", "merged_SD"
                else:
                    s, e = (edge - a.fallback_flank, edge) if side == "left" else (edge, edge + a.fallback_flank)
                    name, src = "none", "no_SD_fixed_window"
            flank = round((e - s) * a.flank_frac) if src != "no_SD_fixed_window" else 0
            row.update({f"{side}_sd_name": name, f"{side}_sd_start": s, f"{side}_sd_end": e, f"{side}_sd_source": src,
                        f"{side}_win_start": max(0, s - flank), f"{side}_win_end": e + flank})
        rows.append(row)
    out = pd.DataFrame(rows).sort_values(["chrom", "gd_start"])
    out.to_csv(a.out, sep="\t", index=False)

    def kb(x, y):
        return f"{(y - x) / 1e3:.0f}kb"
    for r in out.itertuples(index=False):
        print(f"{r.locus:24s} {r.types:7s} {r.chrom}:{r.gd_start}-{r.gd_end}  "
              f"L {r.left_sd_name} {kb(r.left_sd_start, r.left_sd_end)} -> {kb(r.left_win_start, r.left_win_end)}  "
              f"R {r.right_sd_name} {kb(r.right_sd_start, r.right_sd_end)} -> {kb(r.right_win_start, r.right_win_end)}  "
              f"[{r.left_sd_source}/{r.right_sd_source}]", file=sys.stderr)
    print(f"{len(out)} loci written to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
