#!/usr/bin/env python3

"""Losslessly reverse remap_tandem_to_cpx.py, restoring DUP/CPX records to INS.

remap_tandem_to_cpx.py rewrites truvari-remap tandem-style insertions into either
a plain <DUP> (insertion site inside the duplicated region) or a GATK-SV
<CPX>/dDUP record (dispersed duplication), while snapshotting every field it
overwrites into ORIG_* INFO tags. This script consumes those snapshots and
restores the original insertion record exactly, for downstream tools that expect
INS rather than DUP/CPX.

A record is reverted iff it carries ORIG_ALT (written only by the forward
conversion). For each such record it:
  * restores ALT   <- ORIG_ALT (the original inserted sequence, incl. anchor base)
  * restores POS   <- ORIG_POS
  * restores END   <- ORIG_END (via record.stop)
  * restores SVTYPE/SVLEN <- ORIG_SVTYPE/ORIG_SVLEN if those were present in the
    original; otherwise DELETES SVTYPE/SVLEN (they were absent originally)
  * removes the fields the forward step ADDED: CPX_TYPE, CPX_INTERVALS, SOURCE
  * removes the ORIG_* bookkeeping tags

All other columns and INFO fields (ID, REF, QUAL, FILTER, remap_* tags, sample
genotypes, etc.) are never touched by the forward conversion and so are already
identical to the original.

The ORIG_*/CPX_* INFO *header declarations* added by the forward step are, by
default, also stripped from the output header (use --keep-header-lines to leave
them). Header lines for SVTYPE/SVLEN/END are always kept, since those commonly
predate the conversion.
"""

import argparse
import sys

import pysam
from pysam import VariantFile

# Bookkeeping INFO fields written by remap_tandem_to_cpx.py.
ORIG_FIELDS = ("ORIG_POS", "ORIG_ALT", "ORIG_SVTYPE", "ORIG_SVLEN", "ORIG_END")
# Complex-SV fields the forward step adds to dispersed (dDUP) records. A raw
# truvari-remap insertion never carries these, so they are safe to drop on revert.
ADDED_CPX_FIELDS = ("CPX_TYPE", "CPX_INTERVALS", "SOURCE")
# INFO header declarations to strip from the output header (bookkeeping only).
HEADER_STRIP_IDS = set(ORIG_FIELDS) | set(ADDED_CPX_FIELDS)


def get_single_value(value):
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def build_output_header(in_header, strip_lines: bool):
    """Return the header to write. Optionally drop the bookkeeping INFO lines."""
    if not strip_lines:
        return in_header

    new_header = pysam.VariantHeader()
    for rec in in_header.records:
        # rec.type is e.g. 'INFO', 'FORMAT', 'CONTIG', 'GENERIC', 'STRUCTURED'.
        if rec.type == "INFO" and rec.get("ID") in HEADER_STRIP_IDS:
            continue
        new_header.add_line(str(rec).rstrip())
    for sample in in_header.samples:
        new_header.add_sample(sample)
    return new_header


def revert_record(record):
    orig_alt = get_single_value(record.info.get("ORIG_ALT"))
    orig_pos = get_single_value(record.info.get("ORIG_POS"))
    orig_svtype = get_single_value(record.info.get("ORIG_SVTYPE"))
    orig_svlen = get_single_value(record.info.get("ORIG_SVLEN"))
    orig_end = get_single_value(record.info.get("ORIG_END"))

    # Restore POS first, then ALT (a resolved sequence), then END last so the
    # record end is not transiently derived from the symbolic allele.
    if orig_pos is not None:
        record.pos = orig_pos
    record.alts = (orig_alt,)

    if orig_svtype is not None:
        record.info["SVTYPE"] = orig_svtype
    elif "SVTYPE" in record.info:
        del record.info["SVTYPE"]

    if orig_svlen is not None:
        record.info["SVLEN"] = orig_svlen
    elif "SVLEN" in record.info:
        del record.info["SVLEN"]

    if orig_end is not None:
        record.stop = orig_end

    for key in ADDED_CPX_FIELDS:
        if key in record.info:
            del record.info[key]
    for key in ORIG_FIELDS:
        if key in record.info:
            del record.info[key]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Losslessly reverse remap_tandem_to_cpx.py: restore <DUP>/<CPX> records "
            "back to their original INS representation using the ORIG_* snapshots."
        )
    )
    parser.add_argument("--input", required=True, help="VCF/BCF produced by remap_tandem_to_cpx.py")
    parser.add_argument("--output", required=True, help="Output VCF/BCF path")
    parser.add_argument(
        "--keep-header-lines",
        action="store_true",
        help="Keep the ORIG_*/CPX_TYPE/CPX_INTERVALS/SOURCE INFO header declarations "
        "in the output header (default: strip them).",
    )
    args = parser.parse_args()

    reverted = 0
    passed = 0
    strip = not args.keep_header_lines

    with VariantFile(args.input) as vcf_in:
        out_header = build_output_header(vcf_in.header, strip_lines=strip)
        with VariantFile(args.output, "w", header=out_header) as vcf_out:
            for record in vcf_in:
                if get_single_value(record.info.get("ORIG_ALT")) is None:
                    passed += 1
                else:
                    revert_record(record)
                    reverted += 1
                # When the header was rebuilt, reindex the record onto it so any
                # stripped INFO declarations are dropped cleanly on write.
                if strip:
                    record.translate(out_header)
                vcf_out.write(record)

    sys.stderr.write(
        f"Reverted {reverted} DUP/CPX records to INS; passed {passed} records through unchanged.\n"
    )


if __name__ == "__main__":
    main()
