#!/usr/bin/env python3
"""Reference segments for the read recurrence plots of one region: hg38 region +/- flank, and the region lifted to
T2T-CHM13 (UCSC chain file, hg38 -> CHM13; see read_t2t_recurrence.py) +/- flank. Each segment is written as a small
FASTA whose one sequence is named chrom:start-end (1-based), which read_ref_recurrence.py / read_t2t_recurrence.py
accept as --ref / --t2t-ref.

Outputs:
  <out>.hg38.fa(.fai)        hg38 segment
  <out>.t2t.fa(.fai)         T2T segment (or --t2t-region)
  <out>.t2t_region.txt       the T2T segment as chrom:start-end
  <out>.lift.tsv             every lifted point of the region (hg38_pos chain_id t2t_chrom t2t_pos t2t_strand)

Usage:
  prepare_region.py --region chr16:29324175-29638676 --flank 5000 --hg38-ref hg38.fa --t2t-ref chm13v2.0.fa.bgz \\
      --chain grch38-chm13v2.chain --out 16p_bp4 [--t2t-region chr16:...] [--lift-step 1000]
Needs read_ref_recurrence.py and read_t2t_recurrence.py in the same directory.
"""
import argparse
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_ref_recurrence as rr  # noqa: E402
import read_t2t_recurrence as rt  # noqa: E402


def write_segment(ref, chrom, start, end, path):
    seq = pysam.FastaFile(ref).fetch(chrom, start - 1, end).upper()
    with open(path, "w") as fh:
        fh.write(f">{chrom}:{start}-{start + len(seq) - 1}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i + 80] + "\n")
    pysam.faidx(path)
    print(f"{path}: {chrom}:{start:,}-{start + len(seq) - 1:,} ({len(seq):,} bp)", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--region", required=True, help="hg38 chrom:start-end")
    p.add_argument("--flank", type=int, default=5000)
    p.add_argument("--hg38-ref", required=True)
    p.add_argument("--t2t-ref", required=True)
    p.add_argument("--chain", required=True, help="UCSC chain file hg38 -> T2T-CHM13")
    p.add_argument("--t2t-region", help="T2T chrom:start-end to use instead of the lifted region (flank still added)")
    p.add_argument("--lift-step", type=int, default=1000)
    p.add_argument("--out", required=True, help="output prefix")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()

    chrom, rs, re_ = rr.parse_region(a.region)
    write_segment(a.hg38_ref, chrom, max(1, rs - a.flank), re_ + a.flank, a.out + ".hg38.fa")
    chains = rt.load_chains(a.chain, chrom)
    tchrom, t_lo, t_hi, frac = rt.lift_region(chains, a.region, a.lift_step, a.out + ".lift.tsv")
    if a.t2t_region:
        tchrom, ts, te = rr.parse_region(a.t2t_region)
    else:
        ts, te = t_lo + 1, t_hi + 1
    t2t_len = pysam.FastaFile(a.t2t_ref).get_reference_length(tchrom)
    ts, te = max(1, ts - a.flank), min(t2t_len, te + a.flank)
    write_segment(a.t2t_ref, tchrom, ts, te, a.out + ".t2t.fa")
    with open(a.out + ".t2t_region.txt", "w") as fh:
        fh.write(f"{tchrom}:{ts}-{te}\n")
    print(f"lifted fraction {frac}", file=sys.stderr)


if __name__ == "__main__":
    main()
