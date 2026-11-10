#!/usr/bin/env python3
"""Find the split-aligned reads (or assembly contigs) of a region and write one self-contained bundle per read for
read_split_recurrence.py.

A read is split when it has at least --min-alignments primary/supplementary alignments (SA tag followed, secondary
alignments ignored) on the region's chromosome (--all-chroms: on any chromosome), at least one of them overlapping the
region by >= --min-overlap bp. For each split read (at most --max-reads, most alignments first), <outdir>/<out>.<file_id>.tar.gz
holds:
  read.bam (+ .bai)    the read's primary and supplementary alignments
  segments.fa (+ .fai) the reference segments (--segments, or each alignment's span +/- --pad, overlapping ones
                       merged); one sequence per segment, named chrom:start-end
  segments.txt         the segments, comma-separated (for read_split_recurrence.py --segments)
  read_name.txt        the read name
  segdup.bed, genes.gtf  --segdup / --genes features overlapping the segments (empty files when not given)
file_id = read name with characters other than letters, digits, '.', '-', '_' replaced by '_'.

Outputs: the bundles and <out>.split_reads.tsv (read_name file_id n_alignments alignments segments; alignments as
chrom:start-end(strand,primary|supplementary), 1-based).

Usage:
  list_split_reads.py --bam asm_h1.bam --region chr16:28000000-30500000 --ref hg38.fa --out hap1 --outdir bundles \\
      [--pad 20000] [--min-overlap 1000] [--min-alignments 2] [--all-chroms] [--max-reads 50] \\
      [--segments chr16:29274176-29688676,chr16:30138532-30393458] \\
      [--segdup segdups.bed.gz] [--genes genes.gtf.gz]
Needs read_ref_recurrence.py and read_split_recurrence.py in the same directory.
"""
import argparse
import gzip
import os
import re
import shutil
import sys
import tarfile

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_ref_recurrence as rr  # noqa: E402
import read_split_recurrence as rs  # noqa: E402


def read_features(path, kind, chroms):
    """BED (chrom, start0, end) or GTF gene lines (chrom, start1, end) on the given chromosomes, as (chrom, s, e, line)."""
    if not path:
        return []
    out = []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0].startswith("#") or (chroms is not None and f[0] not in chroms):
                continue
            if kind == "bed":
                out.append((f[0], int(f[1]) + 1, int(f[2]), line))
            elif len(f) >= 9 and f[2] == "gene":
                out.append((f[0], int(f[3]), int(f[4]), line))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="BAM (with index), local or gs://")
    p.add_argument("--region", required=True, help="chrom:start-end")
    p.add_argument("--ref", required=True, help="reference FASTA (with .fai)")
    p.add_argument("--out", required=True, help="output prefix (also the bundle name prefix)")
    p.add_argument("--outdir", default=".")
    p.add_argument("--pad", type=int, default=20000)
    p.add_argument("--min-overlap", type=int, default=1000)
    p.add_argument("--min-alignments", type=int, default=2)
    p.add_argument("--all-chroms", action="store_true")
    p.add_argument("--max-reads", type=int, default=50)
    p.add_argument("--max-ref", type=int, default=5000000, help="skip reads whose segments total more than this (bp)")
    p.add_argument("--segments", help="comma-separated chrom:start-end: fixed reference segments for every read (e.g. the "
                   "two SDs of a GD +/- flank) instead of each read's alignment spans +/- --pad")
    p.add_argument("--segdup")
    p.add_argument("--genes")
    p.add_argument("--requester-pays-project")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    rr.gs_env([a.bam, a.ref], a.requester_pays_project)
    os.makedirs(a.outdir, exist_ok=True)

    chrom, rs_, re_ = rr.parse_region(a.region)
    fixed = [rr.parse_region(x) for x in a.segments.split(",")] if a.segments else None
    chroms = None if a.all_chroms else {chrom}
    bam = pysam.AlignmentFile(a.bam)
    # candidates: reads overlapping the region by >= min_overlap with a supplementary alignment somewhere (SA tag)
    overlap, has_sa = {}, set()
    for r in bam.fetch(chrom, rs_ - 1, re_):
        if r.is_secondary or r.is_unmapped:
            continue
        overlap[r.query_name] = overlap.get(r.query_name, 0) + max(0, min(r.reference_end, re_) - max(r.reference_start, rs_ - 1))
        if r.has_tag("SA") or r.is_supplementary:
            has_sa.add(r.query_name)
    split = []
    for name in sorted(n for n in has_sa if overlap[n] >= a.min_overlap):
        alns, _, _ = rr.fetch_read(bam, name, a.region)
        use = rs.used_alignments(alns, chroms)
        if len(use) >= a.min_alignments:
            split.append((name, use))
    split.sort(key=lambda x: -len(x[1]))
    print(f"{len(overlap)} reads overlap {a.region}; {len(split)} split (>= {a.min_alignments} alignments"
          f"{'' if a.all_chroms else ' on ' + chrom})", file=sys.stderr)

    fa = pysam.FastaFile(a.ref)
    segdups = read_features(a.segdup, "bed", chroms)
    genes = read_features(a.genes, "gtf", chroms)
    with open(a.out + ".split_reads.tsv", "w") as summary:
        summary.write("read_name\tfile_id\tn_alignments\talignments\tsegments\tbundled\n")
        for k, (name, use) in enumerate(split):
            fid = re.sub(r"[^A-Za-z0-9._-]", "_", name)
            segs = fixed if fixed else rs.segments_from(use, a.pad)
            total = sum(e - s + 1 for _, s, e in segs)
            desc = ",".join(f"{r.reference_name}:{r.reference_start + 1}-{r.reference_end}"
                            f"({'-' if r.is_reverse else '+'},{'supplementary' if r.is_supplementary else 'primary'})"
                            for r in use)
            seg_txt = ",".join(f"{c}:{s}-{e}" for c, s, e in segs)
            ok = k < a.max_reads and total <= a.max_ref
            summary.write(f"{name}\t{fid}\t{len(use)}\t{desc}\t{seg_txt}\t{ok}\n")
            if not ok:
                continue
            d = os.path.join(a.outdir, f"{a.out}.{fid}")
            os.makedirs(d, exist_ok=True)
            alns, _, _ = rr.fetch_read(bam, name, a.region)
            with pysam.AlignmentFile(os.path.join(d, "unsorted.bam"), "wb", header=bam.header) as o:
                for r in alns:
                    o.write(r)
            pysam.sort("-o", os.path.join(d, "read.bam"), os.path.join(d, "unsorted.bam"))
            os.remove(os.path.join(d, "unsorted.bam"))
            pysam.index(os.path.join(d, "read.bam"))
            with open(os.path.join(d, "segments.fa"), "w") as fh:
                for c, s, e in segs:
                    e = min(e, fa.get_reference_length(c))
                    seq = fa.fetch(c, s - 1, e).upper()
                    fh.write(f">{c}:{s}-{e}\n")
                    for i in range(0, len(seq), 80):
                        fh.write(seq[i:i + 80] + "\n")
            pysam.faidx(os.path.join(d, "segments.fa"))
            names = [n for n in pysam.FastaFile(os.path.join(d, "segments.fa")).references]
            with open(os.path.join(d, "segments.txt"), "w") as fh:
                fh.write(",".join(names) + "\n")
            with open(os.path.join(d, "read_name.txt"), "w") as fh:
                fh.write(name + "\n")
            for fname, feats in (("segdup.bed", segdups), ("genes.gtf", genes)):
                with open(os.path.join(d, fname), "w") as fh:
                    for c, s, e, line in feats:
                        if any(c == sc and e >= ss and s <= se for sc, ss, se in segs):
                            fh.write(line)
            with tarfile.open(d + ".tar.gz", "w:gz") as t:
                t.add(d, arcname=os.path.basename(d))
            shutil.rmtree(d)
            print(f"  {name}: {len(use)} alignments; segments {seg_txt} -> {d}.tar.gz", file=sys.stderr)


if __name__ == "__main__":
    main()
