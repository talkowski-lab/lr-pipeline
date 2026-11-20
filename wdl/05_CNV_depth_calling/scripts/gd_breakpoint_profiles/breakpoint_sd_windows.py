#!/usr/bin/env python3
"""Windows around breakpoint segmental duplications (SDs) of a genomic disorder, from the GD NAHR breakpoint-SD table
(build_GD_breakpoint_SD_bed output: #chrom start end cluster breakpoint ...), each widened on both ends by a fraction of
its own length (default half: window = start - 0.5 x length .. end + 0.5 x length, clipped at 0).

Output BED: chrom, window_start, window_end, name (<cluster>_<breakpoint>), sd_start, sd_end, sd_length, flank (bp);
also printed with a ready-to-use region string (chrom:start-end, 1-based) for tabix.

Usage:
  breakpoint_sd_windows.py --bed references/GD_NAHR.breakpoint_SDs.bed.gz --cluster 15q11-13 --breakpoints BP1,BP2 \\
      --out aou_II/chr15_bp1_bp2.windows.bed [--flank-frac 0.5]
"""
import argparse
import sys

import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", required=True)
    p.add_argument("--cluster", required=True, help="cluster name, e.g. 15q11-13 or 16p11")
    p.add_argument("--breakpoints", required=True, help="comma-separated breakpoint names, e.g. BP1,BP2")
    p.add_argument("--flank-frac", type=float, default=0.5, help="flank on each end as a fraction of the SD length")
    p.add_argument("--out", required=True)
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()

    d = pd.read_csv(a.bed, sep="\t")
    d = d.rename(columns={d.columns[0]: "chrom"})
    want = a.breakpoints.split(",")
    sel = d[(d["cluster"] == a.cluster) & d["breakpoint"].isin(want)].copy()
    missing = set(want) - set(sel["breakpoint"])
    if missing:
        sys.exit(f"not found in cluster {a.cluster}: {', '.join(sorted(missing))}")
    sel["sd_length"] = sel["end"] - sel["start"]
    sel["flank"] = (sel["sd_length"] * a.flank_frac).round().astype(int)
    sel["window_start"] = (sel["start"] - sel["flank"]).clip(lower=0)
    sel["window_end"] = sel["end"] + sel["flank"]
    sel["name"] = sel["cluster"] + "_" + sel["breakpoint"]
    out = sel[["chrom", "window_start", "window_end", "name", "start", "end", "sd_length", "flank"]].rename(
        columns={"start": "sd_start", "end": "sd_end"}).sort_values(["chrom", "window_start"])
    out.to_csv(a.out, sep="\t", index=False, header=False)
    for r in out.itertuples(index=False):
        print(f"{r.name}\tSD {r.chrom}:{r.sd_start:,}-{r.sd_end:,} ({r.sd_length:,} bp)\tflank {r.flank:,} bp\t"
              f"window {r.chrom}:{r.window_start + 1}-{r.window_end}")


if __name__ == "__main__":
    main()
