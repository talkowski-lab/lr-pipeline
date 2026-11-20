#!/usr/bin/env python3
"""List the reads (or assembly contigs) aligned to a region and write one small BAM per read holding all its primary and
supplementary alignments (followed through the SA tag, also outside the region), so each read can be analysed in its
own job with read_ref_recurrence.py / read_t2t_recurrence.py.

A read is kept when its non-secondary alignments overlap the region by >= --min-overlap bp in total; with
--spanning-only it must also cover the region's first and last base. At most --max-reads reads are kept (largest overlap
first).

Outputs:
  <out>.reads.tsv   read_name file_id n_alignments strands overlap_bp aligned_start aligned_end spans_region kept
                    (aligned_start/end: leftmost/rightmost aligned base within the region)
  <outdir>/<out>.<file_id>.bam (+ .bai)   one per kept read; file_id = read name with characters other than
                    letters, digits, '.', '-', '_' replaced by '_'

Usage:
  list_region_reads.py --bam asm_h1.bam --region chr16:29324175-29638676 --out hap1 --outdir reads \\
      [--min-overlap 1000] [--spanning-only] [--max-reads 100]
"""
import argparse
import os
import re
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_ref_recurrence as rr  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True, help="BAM (with index), local or gs://")
    p.add_argument("--region", required=True, help="chrom:start-end")
    p.add_argument("--out", required=True, help="output prefix (also the per-read BAM name prefix)")
    p.add_argument("--outdir", default=".", help="directory for the per-read BAMs")
    p.add_argument("--min-overlap", type=int, default=1000)
    p.add_argument("--spanning-only", action="store_true")
    p.add_argument("--max-reads", type=int, default=100)
    p.add_argument("--requester-pays-project")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    rr.gs_env([a.bam], a.requester_pays_project)
    os.makedirs(a.outdir, exist_ok=True)

    chrom, rs, re_ = rr.parse_region(a.region)
    bam = pysam.AlignmentFile(a.bam)
    reads = {}
    for r in bam.fetch(chrom, rs - 1, re_):
        if r.is_secondary or r.is_unmapped:
            continue
        s, e = max(r.reference_start, rs - 1), min(r.reference_end, re_)
        d = reads.setdefault(r.query_name, {"n": 0, "strands": set(), "ov": 0, "lo": e, "hi": s})
        d["n"] += 1
        d["strands"].add("-" if r.is_reverse else "+")
        d["ov"] += max(0, e - s)
        d["lo"], d["hi"] = min(d["lo"], s), max(d["hi"], e)
    rows = []
    for name, d in reads.items():
        spans = d["lo"] <= rs - 1 and d["hi"] >= re_
        ok = d["ov"] >= a.min_overlap and (spans or not a.spanning_only)
        rows.append([name, re.sub(r"[^A-Za-z0-9._-]", "_", name), d["n"], ",".join(sorted(d["strands"])), d["ov"],
                     d["lo"] + 1, d["hi"], spans, ok])
    rows.sort(key=lambda x: -x[4])
    n_kept = 0
    for row in rows:
        if row[-1] and n_kept >= a.max_reads:
            row[-1] = False
        if not row[-1]:
            continue
        n_kept += 1
        alns, _, _ = rr.fetch_read(bam, row[0], a.region)
        path = os.path.join(a.outdir, f"{a.out}.{row[1]}.bam")
        with pysam.AlignmentFile(path + ".unsorted", "wb", header=bam.header) as o:
            for r in alns:
                o.write(r)
        pysam.sort("-o", path, path + ".unsorted")
        os.remove(path + ".unsorted")
        pysam.index(path)
    with open(a.out + ".reads.tsv", "w") as fh:
        fh.write("read_name\tfile_id\tn_alignments\tstrands\toverlap_bp\taligned_start\taligned_end\tspans_region\tkept\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
    print(f"{len(rows)} reads overlap {a.region}; {n_kept} kept -> {a.outdir}/{a.out}.*.bam", file=sys.stderr)


if __name__ == "__main__":
    main()
