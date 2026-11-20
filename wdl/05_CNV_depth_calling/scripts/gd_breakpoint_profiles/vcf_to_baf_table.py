#!/usr/bin/env python3
"""Per-sample B-allele frequency (BAF) table from a VCF with FORMAT/AD.

BAF = AD[alt] / sum(AD) for each ALT allele of a record (multiallelic records give one row per ALT allele). A genotype
is reported as NA when AD is missing, sum(AD) is 0, or sum(AD) < --min-depth. DP is reported in the long format only.

Output (tab-separated, .gz allowed):
  wide (default)  CHROM POS ID REF ALT <one BAF column per sample>
  long            CHROM POS ID REF ALT sample GT DP AD_ref AD_alt BAF
Sample names are written as is, or with --anonymize as random names (sample_0001, ... in random order); the
name -> random-name key is then written to --key-out (required with --anonymize) so the table can be decoded later.

Example:
  python3 vcf_to_baf_table.py in.vcf.gz samples.txt baf.tsv.gz --anonymize --key-out baf.sample_key.tsv --seed 1
"""
import argparse
import gzip
import random

import pysam


def fmt(x):
    return "NA" if x is None else f"{x:.4f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vcf", help="input VCF / VCF.gz / BCF with FORMAT/GT, DP and AD")
    ap.add_argument("samples", help="file with one sample ID per line (samples of interest)")
    ap.add_argument("output", help="output table (.gz allowed)")
    ap.add_argument("--format", choices=["wide", "long"], default="wide", help="table layout (default: wide)")
    ap.add_argument("--region", help="restrict to a region, e.g. chr22 or chr22:1-1000000 (needs an index)")
    ap.add_argument("--min-depth", type=int, default=0, help="report NA when sum(AD) is below this (default: 0)")
    ap.add_argument("--pass-only", action="store_true", help="keep only records with FILTER PASS (or .)")
    ap.add_argument("--anonymize", action="store_true", help="replace sample IDs with random names")
    ap.add_argument("--key-out", help="with --anonymize: write the sample -> random-name key here")
    ap.add_argument("--prefix", default="sample", help="random-name prefix (default: sample)")
    ap.add_argument("--seed", type=int, help="random seed for --anonymize (default: unseeded)")
    args = ap.parse_args()
    if args.anonymize and not args.key_out:
        ap.error("--anonymize requires --key-out")

    wanted = [s.strip() for s in open(args.samples) if s.strip()]
    vcf = pysam.VariantFile(args.vcf)
    in_vcf = set(vcf.header.samples)
    missing = [s for s in wanted if s not in in_vcf]
    if missing:
        print(f"warning: {len(missing)} samples not in the VCF, skipped: {', '.join(missing[:10])}"
              + (" ..." if len(missing) > 10 else ""))
    samples = list(dict.fromkeys(s for s in wanted if s in in_vcf))
    if not samples:
        raise SystemExit("none of the requested samples are in the VCF")
    vcf.subset_samples(samples)

    names = {s: s for s in samples}
    if args.anonymize:
        order = samples[:]
        random.Random(args.seed).shuffle(order)
        width = max(4, len(str(len(order))))
        names = {s: f"{args.prefix}_{i + 1:0{width}d}" for i, s in enumerate(order)}
        with open(args.key_out, "w") as fh:
            fh.write("sample\trandom_name\n")
            for s in samples:
                fh.write(f"{s}\t{names[s]}\n")

    out = gzip.open(args.output, "wt") if args.output.endswith(".gz") else open(args.output, "w")
    base = ["CHROM", "POS", "ID", "REF", "ALT"]
    if args.format == "wide":
        out.write("\t".join(base + [names[s] for s in samples]) + "\n")
    else:
        out.write("\t".join(base + ["sample", "GT", "DP", "AD_ref", "AD_alt", "BAF"]) + "\n")

    n = 0
    records = vcf.fetch(region=args.region) if args.region else vcf.fetch()
    for rec in records:
        if args.pass_only and list(rec.filter.keys()) not in ([], ["PASS"]):
            continue
        for ai, alt in enumerate(rec.alts or [], start=1):
            lead = [rec.chrom, str(rec.pos), rec.id or ".", rec.ref, alt]
            vals = []
            for s in samples:
                call = rec.samples[s]
                ad = call.get("AD")
                baf, ref_ad, alt_ad = None, None, None
                if ad is not None and len(ad) > ai and ad[0] is not None and ad[ai] is not None:
                    tot = sum(x for x in ad if x is not None)
                    ref_ad, alt_ad = ad[0], ad[ai]
                    if tot > 0 and tot >= args.min_depth:
                        baf = alt_ad / tot
                if args.format == "wide":
                    vals.append(fmt(baf))
                else:
                    gt = call.get("GT")
                    sep = "|" if call.phased else "/"
                    gt_s = sep.join("." if a is None else str(a) for a in gt) if gt else "."
                    dp = call.get("DP")
                    out.write("\t".join(lead + [names[s], gt_s, "NA" if dp is None else str(dp),
                                                "NA" if ref_ad is None else str(ref_ad),
                                                "NA" if alt_ad is None else str(alt_ad), fmt(baf)]) + "\n")
            if args.format == "wide":
                out.write("\t".join(lead + vals) + "\n")
            n += 1
    out.close()
    print(f"{n} variant rows x {len(samples)} samples written to {args.output}")


if __name__ == "__main__":
    main()
