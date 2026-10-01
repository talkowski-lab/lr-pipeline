#!/usr/bin/env python3

import argparse
import re
import sys
from typing import Tuple

from pysam import VariantFile

TARGET_CLASSES = {"tandem", "tandem_inverted", "tandem_complex"}
REMAP_COORDS_RE = re.compile(r"^(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)$")


def get_single_value(value):
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def parse_remap_coords(coords: str) -> Tuple[str, int, int]:
    match = REMAP_COORDS_RE.match(coords)
    if not match:
        raise ValueError(f"Invalid remap_coords value: {coords}")

    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if end < start:
        raise ValueError(f"Invalid remap_coords interval with end < start: {coords}")
    return chrom, start, end


def build_header(header):
    # Mutate the header in place (do NOT copy) so that records iterated from the
    # input file inherit these INFO definitions. If we copied the header instead,
    # the iterated records would still reference the original header, and looking
    # up / setting an undeclared field such as ORIG_POS raises
    # "ValueError: Invalid header".
    if "ORIG_POS" not in header.info:
        header.add_line(
            '##INFO=<ID=ORIG_POS,Number=1,Type=Integer,Description="Original insertion position before remap tandem-to-dup post-processing">'
        )

    if "ORIG_ALT" not in header.info:
        header.add_line(
            '##INFO=<ID=ORIG_ALT,Number=1,Type=String,Description="Original ALT allele (including any leading anchor base) before remap tandem-to-dup post-processing rewrote ALT to <DUP>">'
        )

    if "SVLEN" not in header.info:
        header.add_line(
            '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
        )

    if "END" not in header.info:
        header.add_line(
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">'
        )

    return header


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert remap.py tandem-style insertion calls into DUP records by "
            "rewriting POS/END/SVLEN from remap_coords and storing the original "
            "insertion position in ORIG_POS."
        )
    )
    parser.add_argument("--input", required=True, help="Input VCF/BCF annotated by truvari anno remap")
    parser.add_argument("--output", required=True, help="Output VCF/BCF path")
    args = parser.parse_args()

    converted = 0
    skipped_missing_coords = 0
    skipped_contig_mismatch = 0

    with VariantFile(args.input) as vcf_in:
        header = build_header(vcf_in.header)
        with VariantFile(args.output, "w", header=header) as vcf_out:
            for record in vcf_in:
                remap_classification = get_single_value(record.info.get("remap_classification"))

                if remap_classification in TARGET_CLASSES:
                    remap_coords = get_single_value(record.info.get("remap_coords"))
                    if remap_coords is None:
                        skipped_missing_coords += 1
                        sys.stderr.write(
                            f"Skipping {record.contig}:{record.pos}:{record.id or '.'} - missing remap_coords\n"
                        )
                        vcf_out.write(record)
                        continue

                    remap_chrom, remap_start, remap_end = parse_remap_coords(remap_coords)
                    if remap_chrom != record.contig:
                        skipped_contig_mismatch += 1
                        sys.stderr.write(
                            f"Skipping {record.contig}:{record.pos}:{record.id or '.'} - remap_coords chromosome {remap_chrom} does not match record chromosome\n"
                        )
                        vcf_out.write(record)
                        continue

                    if get_single_value(record.info.get("ORIG_POS")) is None:
                        record.info["ORIG_POS"] = record.pos

                    # Preserve the original ALT allele (including any leading
                    # anchor base) before it is overwritten with the symbolic
                    # <DUP> allele below. Guarded so a re-run (where ALT is
                    # already <DUP>) does not clobber it.
                    orig_alt = get_single_value(record.alts)
                    if (
                        orig_alt is not None
                        and get_single_value(record.info.get("ORIG_ALT")) is None
                    ):
                        record.info["ORIG_ALT"] = orig_alt

                    record.info["SVTYPE"] = "DUP"
                    record.pos = remap_start
                    # pysam treats END as a reserved attribute: it is derived from
                    # record.stop and cannot be assigned via record.info["END"]
                    # (doing so raises KeyError). Setting record.stop updates the
                    # END INFO field on write.
                    record.stop = remap_end
                    # remap_coords is 1-based inclusive on both ends (truvari sets
                    # start = SAM POS and end = POS + ref_advance - 1), so the
                    # duplicated span is (end - start + 1) bases.
                    record.info["SVLEN"] = remap_end - remap_start + 1
                    record.alts = ("<DUP>",)
                    converted += 1

                vcf_out.write(record)

    sys.stderr.write(
        "Converted "
        f"{converted} records; skipped {skipped_missing_coords} without remap_coords and "
        f"{skipped_contig_mismatch} with chromosome mismatches.\n"
    )


if __name__ == "__main__":
    main()
