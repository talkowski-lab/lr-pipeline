#!/usr/bin/env python3
"""Genotypes (GT) of selected samples from a VCF/BCF as a table: CHROM POS ID REF ALT + one GT column per sample, in the
order of the sample list (pysam; no bcftools needed). Missing calls are '.', phased calls use '|'.
Samples listed twice are kept once; samples not in the VCF are skipped with a warning.

Usage:
  extract_gt.py --vcf region.vcf.gz --samples samples.txt --out region.GT.tsv.gz [--region chr22:18601988-19032295]
"""
import argparse
import gzip
import sys

import pysam


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True)
    p.add_argument("--samples", required=True, help="one sample ID per line (first column used)")
    p.add_argument("--out", required=True, help="output table (.gz allowed)")
    p.add_argument("--region", help="chrom:start-end (1-based; needs an index)")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()

    wanted = list(dict.fromkeys(ln.split()[0] for ln in open(a.samples) if ln.strip()))
    vcf = pysam.VariantFile(a.vcf)
    in_vcf = set(vcf.header.samples)
    missing = [s for s in wanted if s not in in_vcf]
    if missing:
        print(f"warning: {len(missing)} samples not in the VCF, skipped: {', '.join(missing[:10])}", file=sys.stderr)
    samples = [s for s in wanted if s in in_vcf]
    vcf.subset_samples(samples)
    if a.region:
        chrom, rng = a.region.split(":")
        start, end = (int(x.replace(",", "")) for x in rng.split("-"))
        records = vcf.fetch(chrom, start - 1, end)
    else:
        records = vcf.fetch()
    n = 0
    with (gzip.open(a.out, "wt") if a.out.endswith(".gz") else open(a.out, "w")) as out:
        out.write("\t".join(["CHROM", "POS", "ID", "REF", "ALT"] + samples) + "\n")
        for rec in records:
            gts = []
            for s in samples:
                c = rec.samples[s]
                gt = c.get("GT")
                if gt is None:
                    gts.append(".")
                    continue
                sep = "|" if c.phased else "/"
                gts.append(sep.join("." if x is None else str(x) for x in gt))
            out.write("\t".join([rec.chrom, str(rec.pos), rec.id or ".", rec.ref, ",".join(rec.alts or ["."])] + gts) + "\n")
            n += 1
    print(f"wrote {a.out}: {n} records x {len(samples)} samples", file=sys.stderr)


if __name__ == "__main__":
    main()
