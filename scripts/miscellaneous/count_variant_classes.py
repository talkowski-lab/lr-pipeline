#!/usr/bin/env python3
"""Count SNVs, indels (1-49 bp) and SVs (>= 50 bp) per contig across a list of VCFs, and their sum over all contigs.

Each ALT allele of a record is classified (multi-allelic records count once per ALT allele; '*' and '.' ALTs are skipped):
  snv    REF and ALT both 1 bp
  indel  1 <= |len(ALT) - len(REF)| <= --max-indel (49)
  sv     |len(ALT) - len(REF)| > --max-indel, or a symbolic / breakend ALT (<DEL>, <INS>, N[chr1:123[, ...)
  other  same-length REF / ALT longer than 1 bp (MNPs, length-neutral complex alleles)
Records are read with `bcftools query` (genotypes are not parsed). VCFs are read in full; no index is needed.

Inputs: --vcfs (one or more VCF / BCF paths) and / or --counts (per-contig tables written by this script, e.g. one per scatter
shard, to add up). Counts of a contig seen in several inputs are added.
Outputs:
  --out-per-contig  contig, records, snv, indel, sv, other, alleles  (one row per contig, in first-seen order)
  --out-total       the same columns summed over all contigs (contig = "total")

Example:
  count_variant_classes.py --vcfs chr1.vcf.gz chr2.vcf.gz --out-per-contig counts.per_contig.tsv --out-total counts.total.tsv
  count_variant_classes.py --counts shard*.per_contig.tsv --out-per-contig all.per_contig.tsv --out-total all.total.tsv
"""
import argparse
import csv
import subprocess
from collections import OrderedDict

COLS = ["records", "snv", "indel", "sv", "other", "alleles"]


def classify(ref, alt, max_indel):
    if alt in ("*", ".") or not alt:
        return None
    if alt.startswith("<") or "[" in alt or "]" in alt:
        return "sv"
    d = abs(len(alt) - len(ref))
    if d == 0:
        return "snv" if len(ref) == 1 else "other"
    return "indel" if d <= max_indel else "sv"


def count_vcf(path, max_indel, counts):
    p = subprocess.Popen(["bcftools", "query", "-f", "%CHROM\t%REF\t%ALT\n", path], stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 20)
    for line in p.stdout:
        chrom, ref, alts = line.rstrip("\n").split("\t")
        c = counts.setdefault(chrom, dict.fromkeys(COLS, 0))
        c["records"] += 1
        for alt in alts.split(","):
            k = classify(ref, alt, max_indel)
            if k:
                c[k] += 1
                c["alleles"] += 1
    if p.wait() != 0:
        raise SystemExit(f"bcftools query failed on {path}")


def add_table(path, counts):
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            c = counts.setdefault(r["contig"], dict.fromkeys(COLS, 0))
            for k in COLS:
                c[k] += int(r[k])


def write(path, rows):
    with open(path, "w") as o:
        o.write("\t".join(["contig"] + COLS) + "\n")
        for name, c in rows:
            o.write("\t".join([name] + [str(c[k]) for k in COLS]) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcfs", nargs="*", default=[], help="VCF / BCF files to count")
    ap.add_argument("--counts", nargs="*", default=[], help="per-contig tables from this script to add")
    ap.add_argument("--max-indel", type=int, default=49, help="largest indel length (bp); longer = SV [49]")
    ap.add_argument("--out-per-contig", required=True)
    ap.add_argument("--out-total", required=True)
    a = ap.parse_args()
    if not a.vcfs and not a.counts:
        ap.error("give --vcfs and / or --counts")
    counts = OrderedDict()
    for path in a.vcfs:
        count_vcf(path, a.max_indel, counts)
        print(f"counted {path}", flush=True)
    for path in a.counts:
        add_table(path, counts)
    total = dict.fromkeys(COLS, 0)
    for c in counts.values():
        for k in COLS:
            total[k] += c[k]
    write(a.out_per_contig, counts.items())
    write(a.out_total, [("total", total)])
    print("\t".join(f"{k}={total[k]}" for k in COLS))


if __name__ == "__main__":
    main()
