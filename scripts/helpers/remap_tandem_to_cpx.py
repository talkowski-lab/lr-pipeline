#!/usr/bin/env python3

"""Convert truvari-remap tandem-style insertions into DUP or GATK-SV CPX/dDUP records.

`truvari anno remap` labels literal-sequence insertions whose ALT sequence maps
back to the reference near the call site as ``tandem`` / ``tandem_inverted``
(a single local duplicated segment) or ``tandem_complex`` (the inserted sequence
reconstructed from multiple stitched local segments). It does NOT rewrite the
records; it only writes the ``remap_*`` INFO tags onto the original INS record.

This script rewrites those tagged insertions into a representation that GATK's
SVAnnotate can functionally annotate. The choice is driven by GEOMETRY, not by
the remap classification label, because the ``tandem``/``tandem_inverted`` test
in remap is a proximity test (nearest edge within one allele length of POS), not
a containment test -- so even a single-segment "tandem" can have its insertion
site OUTSIDE the duplicated segment:

  * Insertion site INSIDE (or within --flank bp of) the duplicated region
      -> plain SVTYPE=DUP / <DUP> over the duplicated span. This matches the
         behavior of remap_tandem_to_dup.py. The copy is adjacent to its source,
         so a single duplication interval captures the consequence.

  * Insertion site OUTSIDE the duplicated region (dispersed duplication)
      -> SVTYPE=CPX / <CPX> with CPX_TYPE=dDUP. The duplicated source region(s)
         go into CPX_INTERVALS as one DUP_ entry per remapped segment, and POS is
         kept at the insertion site. SVAnnotate then annotates each duplicated
         segment with duplication rules AND synthesizes a 1bp insertion segment
         at POS (triggered by CPX_TYPE containing "dDUP"), annotating the
         insertion site independently and unioning the consequences.

The original INS position and ALT are preserved in ORIG_POS / ORIG_ALT.

Caveat: GATK SVAnnotate splits each CPX_INTERVALS entry on "_", so a contig name
containing an underscore (e.g. GRCh38 alt/decoy/*_random contigs) cannot be
represented in CPX_INTERVALS. Such records fall back to the plain <DUP>
conversion with a warning.
"""

import argparse
import re
import sys
from typing import List, Optional, Tuple

from pysam import VariantFile

TARGET_CLASSES = {"tandem", "tandem_inverted", "tandem_complex"}
COORDS_RE = re.compile(r"^(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)$")
# remap_segments entries look like "chr1:10120-10250(+):q1-131", joined by "|".
SEGMENT_RE = re.compile(
    r"^(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)\((?P<ori>[+-])\):q(?P<qs>\d+)-(?P<qe>\d+)$"
)


def get_single_value(value):
    if isinstance(value, (list, tuple)):
        return value[0] if value else None
    return value


def parse_coords(coords: str) -> Tuple[str, int, int]:
    match = COORDS_RE.match(coords)
    if not match:
        raise ValueError(f"Invalid coords value: {coords}")
    chrom = match.group("chrom")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if end < start:
        raise ValueError(f"Invalid interval with end < start: {coords}")
    return chrom, start, end


def parse_segments(segments: str) -> List[Tuple[str, int, int]]:
    """Parse a pipe-joined remap_segments string into [(chrom, start, end), ...].

    Segments are returned sorted by start coordinate. Orientation and query
    coordinates are discarded: for annotating the *source* locus of a dispersed
    duplication, every stitched segment is a duplicated region (DUP), regardless
    of whether the inserted copy is inverted.
    """
    parsed = []
    for token in segments.split("|"):
        token = token.strip()
        if not token:
            continue
        match = SEGMENT_RE.match(token)
        if not match:
            raise ValueError(f"Unparseable remap_segments entry: {token}")
        parsed.append(
            (match.group("chrom"), int(match.group("start")), int(match.group("end")))
        )
    parsed.sort(key=lambda s: (s[1], s[2]))
    return parsed


def build_header(header):
    # Mutate the header in place (do NOT copy) so records iterated from the input
    # inherit these INFO definitions; otherwise setting an undeclared field raises
    # "ValueError: Invalid header".
    additions = [
        ('SVTYPE', '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'),
        ('END', '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">'),
        ('SVLEN', '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'),
        ('CPX_TYPE', '##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Class of complex variant.">'),
        ('CPX_INTERVALS', '##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description="Genomic intervals constituting complex variant.">'),
        ('SOURCE', '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Source of inserted/duplicated sequence.">'),
        ('ORIG_POS', '##INFO=<ID=ORIG_POS,Number=1,Type=Integer,Description="Original insertion POS before remap tandem-to-dup/cpx post-processing (for lossless reversal)">'),
        ('ORIG_ALT', '##INFO=<ID=ORIG_ALT,Number=1,Type=String,Description="Original ALT allele (inserted sequence, incl. anchor base) before ALT was rewritten to a symbolic allele (for lossless reversal)">'),
        ('ORIG_SVTYPE', '##INFO=<ID=ORIG_SVTYPE,Number=1,Type=String,Description="Original SVTYPE INFO value before conversion, if present (for lossless reversal)">'),
        ('ORIG_SVLEN', '##INFO=<ID=ORIG_SVLEN,Number=1,Type=Integer,Description="Original SVLEN INFO value before conversion, if present (for lossless reversal)">'),
        ('ORIG_END', '##INFO=<ID=ORIG_END,Number=1,Type=Integer,Description="Original END (record end) before conversion (for lossless reversal)">'),
    ]
    for key, line in additions:
        if key not in header.info:
            header.add_line(line)

    for alt_id, desc in (("DUP", "Duplication"), ("CPX", "Complex SV")):
        try:
            header.add_line(f'##ALT=<ID={alt_id},Description="{desc}">')
        except Exception:
            pass  # already declared

    return header


def record_key(record) -> str:
    return f"{record.contig}:{record.pos}:{record.id or '.'}"


def preserve_originals(record):
    """Snapshot every field the conversion overwrites, so it can be reversed
    losslessly by remap_cpx_to_ins.py. Each snapshot is guarded so a re-run does
    not clobber a genuine original with an already-converted value.

    SVTYPE/SVLEN are snapshotted only if the input record actually carried them;
    their absence in the ORIG_* tags tells the reverse script to delete (rather
    than restore) those fields. END is always defined (record.stop) so it is
    always snapshotted.
    """
    if get_single_value(record.info.get("ORIG_POS")) is None:
        record.info["ORIG_POS"] = record.pos
    orig_alt = get_single_value(record.alts)
    if orig_alt is not None and get_single_value(record.info.get("ORIG_ALT")) is None:
        record.info["ORIG_ALT"] = orig_alt
    if get_single_value(record.info.get("ORIG_SVTYPE")) is None:
        orig_svtype = get_single_value(record.info.get("SVTYPE"))
        if orig_svtype is not None:
            record.info["ORIG_SVTYPE"] = orig_svtype
    if get_single_value(record.info.get("ORIG_SVLEN")) is None:
        orig_svlen = get_single_value(record.info.get("SVLEN"))
        if orig_svlen is not None:
            record.info["ORIG_SVLEN"] = orig_svlen
    if get_single_value(record.info.get("ORIG_END")) is None:
        record.info["ORIG_END"] = record.stop


def to_dup(record, span_chrom: str, span_start: int, span_end: int):
    """Plain <DUP> over the duplicated span (matches remap_tandem_to_dup.py)."""
    preserve_originals(record)
    record.info["SVTYPE"] = "DUP"
    record.pos = span_start
    # pysam derives END from record.stop; assigning record.info["END"] raises.
    record.stop = span_end
    # remap coords are 1-based inclusive on both ends -> span is end-start+1 bases.
    record.info["SVLEN"] = span_end - span_start + 1
    record.alts = ("<DUP>",)


def to_cpx_ddup(record, segments: List[Tuple[str, int, int]]):
    """<CPX> / dDUP: duplicated segments in CPX_INTERVALS, POS kept at the sink.

    POS is left at the original insertion site so SVAnnotate synthesizes the
    insertion-point segment there. END is normalized to POS (a point sink).
    """
    preserve_originals(record)

    intervals = [f"DUP_{chrom}:{start}-{end}" for chrom, start, end in segments]
    total_dup_bases = sum(end - start + 1 for _chrom, start, end in segments)

    # SOURCE is Number=1: use the overall min-start..max-end span of the source.
    src_chrom = segments[0][0]
    src_start = min(s[1] for s in segments)
    src_end = max(s[2] for s in segments)

    record.info["SVTYPE"] = "CPX"
    record.info["CPX_TYPE"] = "dDUP"
    record.info["CPX_INTERVALS"] = tuple(intervals)
    record.info["SOURCE"] = f"DUP_{src_chrom}:{src_start}-{src_end}"
    record.info["SVLEN"] = total_dup_bases
    # Keep POS at the insertion site; normalize END to a 1bp point sink at POS.
    record.stop = record.pos
    record.alts = ("<CPX>",)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert truvari-remap tandem-style insertions into DUP or CPX/dDUP "
            "records depending on whether the insertion site falls inside the "
            "duplicated region (DUP) or outside it (dispersed duplication -> CPX/dDUP)."
        )
    )
    parser.add_argument("--input", required=True, help="Input VCF/BCF annotated by truvari anno remap")
    parser.add_argument("--output", required=True, help="Output VCF/BCF path")
    parser.add_argument(
        "--flank",
        type=int,
        default=10,
        help=(
            "Tolerance (bp) around the duplicated span within which the insertion "
            "site still counts as 'contained' (plain DUP). Insertion sites farther "
            "outside the span are treated as dispersed duplications (CPX/dDUP). "
            "Default: 10."
        ),
    )
    args = parser.parse_args()

    n_dup = 0
    n_cpx = 0
    n_skip_missing = 0
    n_skip_contig_mismatch = 0
    n_fallback_underscore = 0

    with VariantFile(args.input) as vcf_in:
        header = build_header(vcf_in.header)
        with VariantFile(args.output, "w", header=header) as vcf_out:
            for record in vcf_in:
                cls = get_single_value(record.info.get("remap_classification"))
                if cls not in TARGET_CLASSES:
                    vcf_out.write(record)
                    continue

                # Already converted on a prior run (ORIG_ALT set): pass through
                # untouched so snapshots and coordinates are not double-applied.
                if get_single_value(record.info.get("ORIG_ALT")) is not None:
                    vcf_out.write(record)
                    continue

                coords = get_single_value(record.info.get("remap_coords"))
                if coords is None:
                    n_skip_missing += 1
                    sys.stderr.write(f"Skipping {record_key(record)} - missing remap_coords\n")
                    vcf_out.write(record)
                    continue

                span_chrom, span_start, span_end = parse_coords(coords)
                if span_chrom != record.contig:
                    n_skip_contig_mismatch += 1
                    sys.stderr.write(
                        f"Skipping {record_key(record)} - remap_coords chromosome "
                        f"{span_chrom} does not match record chromosome {record.contig}\n"
                    )
                    vcf_out.write(record)
                    continue

                # Build the list of duplicated source segments. Prefer the
                # per-segment breakdown (present for tandem_complex, and for the
                # refined-misplaced-dup edge case); fall back to the overall span.
                segments_field = get_single_value(record.info.get("remap_segments"))
                segments: Optional[List[Tuple[str, int, int]]] = None
                if segments_field:
                    try:
                        segments = parse_segments(segments_field)
                    except ValueError as exc:
                        sys.stderr.write(
                            f"Warning: {record_key(record)} - {exc}; using remap_coords span\n"
                        )
                        segments = None
                if not segments:
                    segments = [(span_chrom, span_start, span_end)]

                # Geometry: is the insertion site (original POS) inside/adjacent to
                # the duplicated span, or dispersed (outside)?
                ins_pos = get_single_value(record.info.get("ORIG_POS")) or record.pos
                contained = (span_start - args.flank) <= ins_pos <= (span_end + args.flank)

                if contained:
                    to_dup(record, span_chrom, span_start, span_end)
                    n_dup += 1
                else:
                    # CPX_INTERVALS entries are split on "_" by SVAnnotate, so a
                    # contig with "_" cannot be represented. Fall back to <DUP>.
                    if any("_" in chrom for chrom, _s, _e in segments):
                        n_fallback_underscore += 1
                        sys.stderr.write(
                            f"Warning: {record_key(record)} - contig name contains '_'; "
                            f"cannot encode CPX_INTERVALS, falling back to <DUP>\n"
                        )
                        to_dup(record, span_chrom, span_start, span_end)
                        n_dup += 1
                    else:
                        to_cpx_ddup(record, segments)
                        n_cpx += 1

                vcf_out.write(record)

    sys.stderr.write(
        f"Converted {n_dup} records to <DUP> (contained) and {n_cpx} to <CPX>/dDUP "
        f"(dispersed); {n_fallback_underscore} dispersed records fell back to <DUP> "
        f"due to underscore in contig; skipped {n_skip_missing} without remap_coords "
        f"and {n_skip_contig_mismatch} with chromosome mismatches.\n"
    )


if __name__ == "__main__":
    main()
