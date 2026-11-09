#!/usr/bin/env python3
"""Fetch one or more regions of a VCF (local or gs://, requester pays) with pysam and write each as a bgzipped, indexed VCF.

For gs:// paths the OAuth token comes from $GCS_OAUTH_TOKEN, else `gcloud auth print-access-token`; the requester-pays
project from --requester-pays-project, else $GCS_REQUESTER_PAYS_PROJECT, else $GOOGLE_PROJECT.

Usage:
  fetch_vcf_regions.py --vcf gs://.../chr22.g.vcf.bgz --out-prefix 22q11/22q11 \\
      --region A=chr22:18601988-19032295 --region B=chr22:20260680-20416702 [--requester-pays-project PROJECT]
  -> 22q11/22q11_A.vcf.gz (+ .tbi), 22q11/22q11_B.vcf.gz (+ .tbi); a region without NAME= is named chrom_start_end
"""
import argparse
import os
import subprocess
import sys

import pysam


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="VCF/BCF path (local or gs://), with its index next to it")
    p.add_argument("--region", action="append", required=True, help="[NAME=]chrom:start-end (1-based), repeatable")
    p.add_argument("--out-prefix", required=True, help="outputs are <prefix>_<NAME>.vcf.gz")
    p.add_argument("--requester-pays-project")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()

    if a.vcf.startswith("gs://"):
        if not os.environ.get("GCS_OAUTH_TOKEN"):
            os.environ["GCS_OAUTH_TOKEN"] = subprocess.check_output(["gcloud", "auth", "print-access-token"], text=True).strip()
        rp = a.requester_pays_project or os.environ.get("GCS_REQUESTER_PAYS_PROJECT") or os.environ.get("GOOGLE_PROJECT")
        if rp:
            os.environ["GCS_REQUESTER_PAYS_PROJECT"] = rp
    d = os.path.dirname(a.out_prefix)
    if d:
        os.makedirs(d, exist_ok=True)

    vcf = pysam.VariantFile(a.vcf)
    for spec in a.region:
        name, region = spec.split("=", 1) if "=" in spec else (None, spec)
        chrom, rng = region.split(":")
        start, end = (int(x.replace(",", "")) for x in rng.split("-"))
        name = name or f"{chrom}_{start}_{end}"
        if chrom not in vcf.header.contigs:
            sys.exit(f"{chrom} is not a contig of {a.vcf}")
        out_path = f"{a.out_prefix}_{name}.vcf.gz"
        out = pysam.VariantFile(out_path, "wz", header=vcf.header)
        n = 0
        for rec in vcf.fetch(chrom, start - 1, end):
            out.write(rec)
            n += 1
        out.close()
        pysam.tabix_index(out_path, preset="vcf", force=True)
        print(f"{name} {chrom}:{start}-{end}: {n} records -> {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
