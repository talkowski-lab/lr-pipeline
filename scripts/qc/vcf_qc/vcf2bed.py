#!/usr/bin/env python3

# Copyright (c) 2009-2026, Broad Institute, Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the conditions in GATK-SV's
# BSD 3-Clause license are met. See:
# https://github.com/broadinstitute/gatk-sv/blob/main/LICENSE.TXT
#
# Derived from broadinstitute/gatk-sv (branch kj_project_gnomad_lr) src/svtk/svtk/cli/utils.py and
# src/svtk/svtk/utils/utils.py, reduced to the `svtk vcf2bed -i ALL` behavior that collectQC.vcf_wide.sh uses

"""
Convert a VCF to a BED with called samples and every INFO field as columns.
"""

import argparse

import pysam

NULL_GT = [(0, 0), (None, None), (0,), (None,)]


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


def bed_entries(vcf, info_keys):
    for record in vcf:
        start, end = sorted([int(record.pos), int(record.stop)])
        # VCF is 1-based, BED is 0-based
        start = max([0, start - 1])
        svtype = record.alts[0].strip("<>")
        samples = ",".join(get_called_samples(record))
        infos = "\t".join(format_info(record.stop if key == "END" else record.info.get(key)) for key in info_keys)

        record_svtype = record.info.get("SVTYPE", None)
        if record_svtype == "BND":
            start = max([0, int(record.pos) - 1])
            end = record.pos
        elif record_svtype == "CTX":
            start = max([0, int(record.pos) - 1])
            end = start + 1
        elif record_svtype != "INS" and start == end:
            end += 1
        yield f"{record.chrom}\t{start}\t{end}\t{record.id}\t{svtype}\t{samples}\t{infos}\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", help="VCF to convert.")
    parser.add_argument("bed", help="Converted BED.")
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    info_keys = list(vcf.header.info.keys())
    with open(args.bed, "w") as out:
        out.write("\t".join("#chrom start end name svtype samples".split() + info_keys) + "\n")
        for entry in bed_entries(vcf, info_keys):
            out.write(entry)


if __name__ == "__main__":
    main()
