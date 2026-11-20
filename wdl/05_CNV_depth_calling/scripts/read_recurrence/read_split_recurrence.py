#!/usr/bin/env python3
"""Recurrence (dot) plot data for a read or assembly contig with a split alignment (primary + supplementary alignments,
e.g. the two sides of a large deletion), against all its aligned reference segments together.

1. The read is fetched by name in --region; its other alignments are found through the SA tag (as in
   read_ref_recurrence.py). Its primary and supplementary alignments on the chromosomes of --region (--all-chroms: on
   any chromosome) define the reference segments: each alignment's reference span +/- --pad, overlapping segments merged
   (or --segments). Secondary alignments are not used.
2. The read segment is the span of those alignments on the read +/- --pad (or --read-range), in original read
   coordinates.
3. Every --window bp window of the read segment is matched against every reference segment (edit distance <=
   --max-ed-frac x window, both strands; exact seed-and-verify as in read_ref_recurrence.py). --orient-to-ref
   reverse-complements the read when most single-copy hits are on the reverse strand.

Outputs (plot with plot_read_ref_recurrence.R, which draws the segments side by side):
  <out>.hits.tsv    read_start read_end strand ref_start ref_end edit_distance identity segment (1-based index)
  <out>.aln.tsv     the read's alignment blocks (alignment strand read_start read_end chrom ref_start ref_end)
  <out>.info.tsv    read, segments (chrom:start-end;...), read segment, settings, each alignment used
  <out>.alignments.tsv  the alignments used: alignment strand chrom ref_start ref_end read_start read_end mapq

Usage:
  read_split_recurrence.py --bam asm_h1.bam --read-name h1tg003138l --region chr16:29268226-29300000 \\
      --ref hg38.fa --out 16p_BP4-5.h1tg003138l --pad 20000 --threads 8 --orient-to-ref
      [--segments chr16:29240000-29620000,chr16:30160000-30600000] [--all-chroms] [--window 100] [--step 100]
Needs read_ref_recurrence.py in the same directory.
"""
import argparse
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_ref_recurrence as rr  # noqa: E402


def used_alignments(alns, chroms):
    """Primary/supplementary alignments on the given chromosomes (all when chroms is None), sorted by position."""
    return sorted((r for r in alns if chroms is None or r.reference_name in chroms),
                  key=lambda r: (r.reference_name, r.reference_start))


def segments_from(use, pad):
    """Reference segments: each alignment's span +/- pad (1-based, inclusive), overlapping segments merged."""
    segs = []
    for c, s, e in sorted((r.reference_name, max(1, r.reference_start + 1 - pad), r.reference_end + pad) for r in use):
        if segs and segs[-1][0] == c and s <= segs[-1][2]:
            segs[-1] = (c, segs[-1][1], max(segs[-1][2], e))
        else:
            segs.append((c, s, e))
    return segs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    rr.common_args(p)
    p.add_argument("--ref", required=True, help="reference FASTA (with .fai), whole genome or chrom:start-end segments")
    p.add_argument("--segments", help="comma-separated chrom:start-end reference segments instead of the alignments' spans")
    p.add_argument("--all-chroms", action="store_true", help="use alignments on every chromosome, not only --region's")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    a.step = a.step or a.window
    rr.gs_env([a.bam, a.ref], a.requester_pays_project)

    bam = pysam.AlignmentFile(a.bam, reference_filename=a.ref if a.bam.endswith(".cram") else None)
    alns, prim, seq = rr.fetch_read(bam, a.read_name, a.region)
    qlen = len(seq)
    chroms = None if a.all_chroms or not a.region else {rr.parse_region(a.region)[0]}
    use = used_alignments(alns, chroms)
    if not use:
        sys.exit("no alignment of the read on the chromosome of --region")
    spans = []
    with open(a.out + ".alignments.tsv", "w") as fh:
        fh.write("alignment\tstrand\tchrom\tref_start\tref_end\tread_start\tread_end\tmapq\n")
        for r in use:
            q0, q1 = rr.read_coords(r, r.query_alignment_start, r.query_alignment_end)
            spans.append((q0, q1))
            kind = "supplementary" if r.is_supplementary else "primary"
            strand = "-" if r.is_reverse else "+"
            fh.write(f"{kind}\t{strand}\t{r.reference_name}\t{r.reference_start + 1}\t{r.reference_end}\t{q0 + 1}\t{q1}\t"
                     f"{r.mapping_quality}\n")
            print(f"  {kind} {strand} {r.reference_name}:{r.reference_start + 1:,}-{r.reference_end:,} -> read "
                  f"{q0 + 1:,}-{q1:,}", file=sys.stderr)

    # Reference segments: alignment spans +/- pad, merged when they overlap
    if a.segments:
        segs = [rr.parse_region(x) for x in a.segments.split(",")]
    else:
        segs = segments_from(use, a.pad)
    total = sum(e - s + 1 for _, s, e in segs)
    if total > a.max_ref:
        sys.exit(f"reference segments total {total:,} bp, more than --max-ref {a.max_ref:,}; use --segments")

    if a.read_range:
        q_lo, q_hi = (int(x.replace(",", "")) for x in a.read_range.split("-"))
        q_lo -= 1
    else:
        q_lo, q_hi = min(s for s, _ in spans) - a.pad, max(e for _, e in spans) + a.pad
    q_lo, q_hi = max(0, q_lo), min(qlen, q_hi)
    print(f"read {a.read_name}: {qlen:,} bp; read segment {q_lo + 1:,}-{q_hi:,}; reference segments "
          + ", ".join(f"{c}:{s:,}-{e:,}" for c, s, e in segs), file=sys.stderr)

    rows = []
    for k, (c, s, e) in enumerate(segs, 1):
        ref = rr.encode(rr.fetch_ref(a.ref, c, s, e))
        part, max_ed, _ = rr.run_windows(seq, q_lo, q_hi, ref, s, a.window, a.step, a.max_ed_frac, a.min_sep, a.threads)
        rows += [r + (k,) for r in part]
        print(f"segment {k} {c}:{s:,}-{e:,}: {len(part)} hits", file=sys.stderr)
    blocks = rr.alignment_blocks(alns)
    seg_span, flipped = (q_lo + 1, q_hi), False
    if a.orient_to_ref:
        rows, blocks, seg_span, flipped = rr.orient_to_ref(rows, blocks, q_lo, q_hi, qlen)
    with open(a.out + ".hits.tsv", "w") as fh:
        fh.write("read_start\tread_end\tstrand\tref_start\tref_end\tedit_distance\tidentity\tsegment\n")
        for r in sorted(rows):
            fh.write("\t".join(map(str, r)) + "\n")
    with open(a.out + ".aln.tsv", "w") as fh:
        fh.write("alignment\tstrand\tread_start\tread_end\tchrom\tref_start\tref_end\n")
        for b in blocks:
            fh.write("\t".join(map(str, b)) + "\n")
    info = [("read_name", a.read_name), ("read_length", qlen), ("reference", os.path.basename(a.ref)),
            ("chrom", segs[0][0]), ("segments", ";".join(f"{c}:{s}-{e}" for c, s, e in segs)),
            ("window", a.window), ("step", a.step), ("read_seg_start", seg_span[0]), ("read_seg_end", seg_span[1]),
            ("max_edit_distance", int(round(a.max_ed_frac * a.window))),
            ("read_orientation", "reverse_complemented" if flipped else "original"),
            ("primary", f"{prim.reference_name}:{prim.reference_start + 1}-{prim.reference_end} "
                        f"({'-' if prim.is_reverse else '+'})")]
    with open(a.out + ".info.tsv", "w") as fh:
        fh.write("key\tvalue\n")
        for key, v in info:
            fh.write(f"{key}\t{v}\n")
    print(f"{len(rows)} hits -> {a.out}.hits.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
