#!/usr/bin/env python3
# Derived from scripts/reformat_variants_for_visualization.py of yuliamostovoy/variant-interpretation
# https://github.com/yuliamostovoy/variant-interpretation/blob/main/scripts/reformat_variants_for_visualization.py
"""
Turn a candidate-variant list (generated from a VCF, typically via `bcftools query`) into
the canonical bgzipped 6-column BED the long-read visualization workflow consumes.

REQUIRED INPUT: a tab-separated file with these 6 columns, in this order:

    1. chrom    e.g. chr1
    2. start    0-based start (use %POS0 in bcftools query)
    3. end      end coordinate (use %END in bcftools query)
    4. ID       variant ID (use %ID; "." is allowed and auto-filled)
    5. svtype   SV type (use %INFO/SVTYPE; "." allowed, becomes VAR — fine for SNV/indel)
    6. samples  comma-separated IDs of the samples to plot this variant for

A header line is optional (auto-detected and skipped). Example generator:

    bcftools view -i 'GT="alt"' calls.vcf.gz \\
      | bcftools query -f '%CHROM\\t%POS0\\t%END\\t%ID\\t%INFO/SVTYPE\\t[%SAMPLE,]\\n' \\
      > candidates.tsv

(the trailing comma left by `[%SAMPLE,]` is stripped automatically.)

The `samples` IDs must match the individual IDs in the ped file and the sample->BAM map.
Output is plain TSV (chrom,start,end,ID,svtype,samples + header); the WDL task bgzips it.
"""

import argparse
import sys


def is_header(fields):
    # header if the "start" column isn't an integer
    if len(fields) < 2:
        return True
    try:
        int(fields[1])
        return False
    except ValueError:
        return True


def clean_samples(raw):
    # accept comma/semicolon/space separators; drop empties; strip any "=GT" suffix
    for sep in (";", " "):
        raw = raw.replace(sep, ",")
    out = []
    for tok in raw.split(","):
        tok = tok.strip().split("=")[0]
        if tok:
            out.append(tok)
    return ",".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="6-column TSV (chrom,start,end,ID,svtype,samples); header optional")
    ap.add_argument("--output", default="variants_for_visualization.bed", help="output BED (plain TSV; bgzipped by the WDL task)")
    args = ap.parse_args()

    rows = []
    with open(args.input) as fh:
        for lineno, line in enumerate(fh):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if lineno == 0 and is_header(f):
                continue
            if len(f) < 6:
                sys.exit(
                    f"ERROR: line {lineno + 1} has {len(f)} columns, expected 6 "
                    f"(chrom,start,end,ID,svtype,samples): {line.rstrip()}"
                )
            chrom, start, end, vid, svtype, samples = f[:6]
            start, end = int(start), int(end)
            samples = clean_samples(samples)
            svtype = svtype.strip().upper()
            if svtype in ("", "."):
                svtype = "VAR"
            if vid.strip() in ("", "."):
                vid = f"{chrom}_{start}_{end}_{svtype}"
            rows.append((chrom, start, end, vid, svtype, samples))

    rows.sort(key=lambda r: (r[0], r[1]))

    with open(args.output, "w") as out:
        out.write("#chrom\tstart\tend\tID\tsvtype\tsamples\n")
        for chrom, start, end, vid, svtype, samples in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{vid}\t{svtype}\t{samples}\n")

    sys.stderr.write(f"Wrote {len(rows)} variants to {args.output}\n")


if __name__ == "__main__":
    main()
