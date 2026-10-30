#!/usr/bin/env python3
"""Give every record of a biallelic VCF an allele type, a signed allele length and a unique ID, for EvaluateCallsetAccuracy.

INFO/allele_type and INFO/allele_length are kept when a record carries them. Otherwise the type comes from SVTYPE,
then from REF and ALT (an SNV or MNV when they are the same length, else an insertion or deletion), and the length
from SVLEN, then from the ALT minus REF length difference; deletion lengths are negative. A duplication written as
an insertion, with a sequence ALT longer than its REF, is typed as an insertion and its SVTYPE rewritten to INS, so
it counts and matches as one. With --replace_ids every record is renamed <id_prefix>_<n>; otherwise only records
whose ID is missing or repeats an earlier one are.
"""

import argparse
from collections import Counter

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", required=True, help="Biallelic input VCF")
    parser.add_argument("--output", required=True, help="Output bgzipped VCF, indexed alongside")
    parser.add_argument("--id_prefix", required=True, help="Prefix of the IDs given to renamed records")
    parser.add_argument("--replace_ids", action="store_true", help="Rename every record rather than only missing or repeated IDs")
    return parser.parse_args()


def info_value(record, key):
    """A record's INFO value, its first when a list, or None when the record or the header lacks the field."""
    if key not in record.header.info or key not in record.info:
        return None
    value = record.info[key]
    return value[0] if isinstance(value, tuple) else value


def allele(record):
    """The record's allele type and signed length."""
    ref, alt = record.ref, record.alts[0]
    existing_type, svtype = info_value(record, "allele_type"), info_value(record, "SVTYPE")
    if existing_type is not None:
        allele_type = str(existing_type).lower()
    elif svtype is not None:
        allele_type = str(svtype).lower()
    elif len(alt) == len(ref):
        allele_type = "snv" if len(ref) == 1 else "mnv"
    else:
        allele_type = "ins" if len(alt) > len(ref) else "del"
    if allele_type == "dup" and not alt.startswith("<") and len(alt) > len(ref):
        allele_type = "ins"

    existing_length = info_value(record, "allele_length")
    if existing_length is not None:
        return allele_type, int(existing_length)
    if allele_type in ("snv", "mnv"):
        return allele_type, 0
    svlen = info_value(record, "SVLEN")
    length = abs(svlen) if svlen is not None else abs(len(alt) - len(ref))
    return allele_type, -length if allele_type == "del" else length


def main():
    args = parse_args()
    types = Counter()
    renamed = 0
    seen = set()
    with pysam.VariantFile(args.vcf) as vcf_in:
        header = vcf_in.header.copy()
        if "allele_type" not in header.info:
            header.add_line('##INFO=<ID=allele_type,Number=1,Type=String,Description="Allele type">')
        if "allele_length" not in header.info:
            header.add_line('##INFO=<ID=allele_length,Number=1,Type=Integer,Description="Allele length">')
        with pysam.VariantFile(args.output, "wz", header=header) as vcf_out:
            for n, record in enumerate(vcf_in, 1):
                record.translate(header)
                allele_type, length = allele(record)
                if allele_type == "ins" and info_value(record, "SVTYPE") == "DUP":
                    record.info["SVTYPE"] = "INS"
                record.info["allele_type"] = allele_type
                record.info["allele_length"] = length
                if args.replace_ids or record.id is None or record.id in seen:
                    record.id = f"{args.id_prefix}_{n}"
                    renamed += 1
                seen.add(record.id)
                types[allele_type] += 1
                vcf_out.write(record)
    pysam.tabix_index(args.output, preset="vcf", force=True)
    print(f"Typed {sum(types.values())} records {dict(types)}, renaming {renamed}")


if __name__ == "__main__":
    main()
