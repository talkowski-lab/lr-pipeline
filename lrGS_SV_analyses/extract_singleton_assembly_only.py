#!/usr/bin/env python3
"""Extract singleton (AC=1) variants whose sole carrier's genotype is
supported ONLY by assembly-based callers (dipcall, hapdiff -- never kanpig
or any read-based caller), from a single chromosome VCF.

A variant qualifies when:
  - INFO/AC == 1 for some ALT allele (site-level singleton)
  - the one sample carrying that ALT allele has FORMAT/EV present, and every
    caller listed in EV is in {dipcall, hapdiff}

Output columns: ID, allele_type, allele_length, dbSNP_ID, gnomAD_V4_match_ID, FILTER

Usage: python3 extract_singleton_assembly_only.py <in.vcf.gz> <out.tsv>
"""
import argparse
import re
import sys

import pysam

ASSEMBLY_CALLERS = {"dipcall", "hapdiff"}
CALLER_RE = re.compile(r"^([A-Za-z0-9]+)")


def parse_callers(ev_val):
    if ev_val is None:
        return None
    if isinstance(ev_val, tuple):
        parts = []
        for v in ev_val:
            if v is None:
                continue
            parts.extend(str(v).split(","))
    else:
        s = str(ev_val)
        if s in (".", ""):
            return None
        parts = s.split(",")
    callers = set()
    for p in parts:
        p = p.strip()
        if not p:
            continue
        m = CALLER_RE.match(p)
        if m:
            callers.add(m.group(1).lower())
    return callers if callers else None


def info_str(v):
    if v is None:
        return "."
    if isinstance(v, tuple):
        vals = [str(x) for x in v if x is not None]
        return ",".join(vals) if vals else "."
    return str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("vcf", help="Input VCF/BCF path (local or remote; index alongside if remote)")
    ap.add_argument("out_tsv", help="Output TSV path")
    ap.add_argument("--progress-every", type=int, default=500000,
                     help="Log progress to stderr every N records (0 to disable)")
    args = ap.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    n_total = 0
    n_hits = 0
    with open(args.out_tsv, "w") as out:
        out.write("ID\tallele_type\tallele_length\tdbSNP_ID\tgnomAD_V4_match_ID\tFILTER\n")
        for rec in vcf:
            n_total += 1
            if args.progress_every and n_total % args.progress_every == 0:
                print(f"...{n_total} records, {n_hits} hits", file=sys.stderr)

            ac = rec.info.get("AC")
            if ac is None:
                continue
            if not isinstance(ac, tuple):
                ac = (ac,)
            target_alt_indices = [i for i, v in enumerate(ac) if v == 1]
            if not target_alt_indices:
                continue

            for alt_idx in target_alt_indices:
                allele_num = alt_idx + 1
                carrier = None
                for sample_id, sample_val in rec.samples.items():
                    gt = sample_val["GT"]
                    if gt is None:
                        continue
                    if any(a == allele_num for a in gt):
                        carrier = sample_val
                        break
                if carrier is None:
                    continue
                callers = parse_callers(carrier.get("EV"))
                if not callers or not callers.issubset(ASSEMBLY_CALLERS):
                    continue

                vid = rec.id if rec.id else f"{rec.chrom}-{rec.pos}-{rec.ref}-{rec.alts[alt_idx]}"
                allele_type = info_str(rec.info.get("allele_type"))
                allele_length = info_str(rec.info.get("allele_length"))
                dbsnp = info_str(rec.info.get("dbSNP_ID"))
                gnomad = info_str(rec.info.get("gnomAD_V4_match_ID"))
                filt = ",".join(rec.filter.keys()) if rec.filter.keys() else "."
                out.write(f"{vid}\t{allele_type}\t{allele_length}\t{dbsnp}\t{gnomad}\t{filt}\n")
                n_hits += 1

    print(f"DONE {args.vcf}: total_records={n_total} hits={n_hits}", file=sys.stderr)


if __name__ == "__main__":
    main()
