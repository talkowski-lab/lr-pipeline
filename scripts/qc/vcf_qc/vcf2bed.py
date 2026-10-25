#!/usr/bin/env python3

# Copyright (c) 2009-2026, Broad Institute, Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the conditions in GATK-SV's
# BSD 3-Clause license are met. See:
# https://github.com/broadinstitute/gatk-sv/blob/main/LICENSE.TXT
#
# Derived from src/svtk/svtk/cli/utils.py of broadinstitute/gatk-sv (branch kj_project_gnomad_lr)
# https://github.com/broadinstitute/gatk-sv/blob/main/src/svtk/svtk/cli/utils.py
# Derived from src/svtk/svtk/utils/utils.py of broadinstitute/gatk-sv (branch kj_project_gnomad_lr)
# https://github.com/broadinstitute/gatk-sv/blob/main/src/svtk/svtk/utils/utils.py

"""
Convert a VCF to a BED, optionally with called samples, INFO fields and FILTER as columns.

A pysam-only port of `svtk vcf2bed` covering the options this repository uses.
"""

import argparse
import functools

import pysam

NULL_GT = [(0, 0), (None, None), (0,), (None,)]
CPX_INS_CLASSES = ["dDUP", "dDUP_iDEL", "INS_iDEL"]


def get_called_samples(record):
    samples = [sample for sample in record.samples.keys() if record.samples[sample]["GT"] not in NULL_GT]
    if record.info.get("SVTYPE", None) == "CNV":
        for sample, gt in record.samples.items():
            cn = gt.get("CN", gt.get("RD_CN", None))
            if cn is not None and cn != 2:
                samples.append(sample)
    return sorted(samples)


def format_info(info):
    if info is None:
        return "NA"
    elif isinstance(info, (tuple, list)):
        return ",".join([str(x) for x in info])
    else:
        return str(info)


def parse_bnd_pos(alt):
    alt = alt.strip("ATCGNRYSWKMBDHV")
    # Strip brackets separately, otherwise GL contigs would be altered
    alt = alt.strip("[]")
    # HLA contigs include colons
    data = alt.split(":")
    return ":".join(data[:-1]).replace("CHR", "chr"), int(data[-1])


def format_entry(chrom, start, end, svtype, name, suffix):
    return f"{chrom}\t{start}\t{end}\t{name}\t{svtype}{suffix}\n"


def bed_entries(vcf, include_samples, info_keys, include_filters, split_bnd, split_cpx):
    for record in vcf:
        columns = []
        if include_samples:
            columns.append(",".join(get_called_samples(record)))
        if info_keys:
            columns.append("\t".join(format_info(record.stop if key == "END" else record.info.get(key)) for key in info_keys))
        if include_filters:
            columns.append(",".join(record.filter))
        suffix = "".join(f"\t{column}" for column in columns)
        entry = functools.partial(format_entry, name=record.id, suffix=suffix)

        start, end = sorted([int(record.pos), int(record.stop)])
        # VCF is 1-based, BED is 0-based
        start = max([0, start - 1])
        svtype = record.alts[0].strip("<>")
        record_svtype = record.info.get("SVTYPE", None)

        if record_svtype == "BND":
            yield entry(record.chrom, max([0, int(record.pos) - 1]), record.pos, svtype)
            if split_bnd:
                if "[" in record.alts[0] or "]" in record.alts[0]:
                    chrom, end = parse_bnd_pos(record.alts[0])
                else:
                    chrom, end = record.info.get("CHR2", None), record.info.get("END2", None)
                yield entry(chrom, max([0, int(end) - 1]), end, svtype)
        elif record_svtype == "INS":
            yield entry(record.chrom, start, end, svtype)
        elif record_svtype == "CTX":
            start = max([0, int(record.pos) - 1])
            yield entry(record.chrom, start, start + 1, svtype)
            if split_bnd:
                yield entry(record.info["CHR2"], max([0, int(record.stop) - 1]), record.stop, svtype)
        elif "CPX_INTERVALS" in record.info and split_cpx:
            for interval in record.info["CPX_INTERVALS"]:
                interval_svtype, region = interval.split("_")
                chrom, coords = region.split(":")
                interval_start, interval_end = coords.split("-")
                yield entry(chrom, max([0, int(interval_start) - 1]), interval_end, interval_svtype)
            # A complex insertion also reports its insertion point as a 1 bp DEL
            if record.info.get("CPX_TYPE", None) in CPX_INS_CLASSES:
                yield entry(record.chrom, max([0, int(record.pos) - 1]), record.pos, "DEL")
        else:
            if start == end:
                end += 1
            yield entry(record.chrom, start, end, svtype)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", help="VCF to convert.")
    parser.add_argument("bed", help="Converted BED.")
    parser.add_argument(
        "--no-samples",
        dest="include_samples",
        action="store_false",
        help="Don't include the comma-delimited list of called samples for each variant.",
    )
    parser.add_argument(
        "-i",
        "--info",
        action="append",
        help="INFO field to include as a column, in the order requested. May be given more than once; `ALL` includes "
        "every INFO field in header order.",
    )
    parser.add_argument("--include-filters", action="store_true", help="Include FILTER as the last column.")
    parser.add_argument("--split-bnd", action="store_true", help="Report two entries for each BND, one per breakend.")
    parser.add_argument("--split-cpx", action="store_true", help="Report an entry for each CPX rearrangement interval.")
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    info_keys = list(vcf.header.info.keys()) if args.info and "ALL" in args.info else args.info or []
    header = "#chrom start end name svtype".split()
    if args.include_samples:
        header.append("samples")
    header += info_keys
    if args.include_filters:
        header.append("FILTER")
    with open(args.bed, "w") as out:
        out.write("\t".join(header) + "\n")
        for line in bed_entries(vcf, args.include_samples, info_keys, args.include_filters, args.split_bnd, args.split_cpx):
            out.write(line)


if __name__ == "__main__":
    main()
