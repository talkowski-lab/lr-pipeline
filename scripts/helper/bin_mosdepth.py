#!/usr/bin/env python3
"""Bin run-length encoded mosdepth coverage with base-weighted medians."""

import argparse
import gzip
from collections import Counter
from pathlib import Path
from statistics import median


def _open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def _open_output(path):
    return gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "w")


def _parse_depth(value, truncate_depth):
    depth = float(value)
    return int(depth) if truncate_depth else depth


def _parse_record(line, line_number, truncate_depth):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 4:
        raise ValueError(f"Expected at least 4 columns at line {line_number}")
    chrom, start_text, end_text, depth_text = fields[:4]
    start = int(start_text)
    end = int(end_text)
    depth = _parse_depth(depth_text, truncate_depth)
    if start < 0 or end <= start or depth < 0:
        raise ValueError(f"Invalid interval at line {line_number}: {line.rstrip()}")
    return chrom, start, end, depth


def _write_row(output, chrom, start, end, depth, coordinate_system):
    if coordinate_system == "one-based":
        start += 1
    output.write(f"{chrom}\t{start}\t{end}\t{_format_depth(depth)}\n")


def _format_depth(depth):
    if isinstance(depth, float) and depth.is_integer():
        return str(int(depth))
    return f"{depth:.10g}" if isinstance(depth, float) else str(depth)


def _write_histogram(histogram, path):
    with open(path, "w") as output:
        for depth in sorted(histogram):
            output.write(f"{depth:g}\t{histogram[depth]}\n")


def bin_mosdepth(
    input_path,
    output_path,
    bin_size,
    coordinate_system="zero-based",
    contig=None,
    output_contig=None,
    allow_partial=False,
    require_contiguous=False,
    skip_comments=False,
    exclude_contigs=(),
    truncate_depth=False,
    histogram_path=None,
):
    """Write fixed-width median bins and return histogram of emitted medians.

    In non-contiguous mode, bins are emitted only when they contain exactly
    ``bin_size`` covered bases. In contiguous mode, bins are emitted as soon as
    crossed and a final partial bin is emitted only when ``allow_partial`` is
    true. The latter matches IdentifyLowCoverageRegions' chromosome streaming
    behavior.
    """
    if bin_size <= 0:
        raise ValueError("bin_size must be greater than zero")
    if coordinate_system not in {"zero-based", "one-based"}:
        raise ValueError("coordinate_system must be zero-based or one-based")

    excluded = set(exclude_contigs)
    histogram = Counter()
    contigs = []
    bins = {}
    previous_end = None
    current_chrom = None
    coverage_counts = Counter()
    covered_bases = 0
    current_start = 0
    current_end = bin_size

    def finish_contiguous(output, end):
        nonlocal coverage_counts, covered_bases
        if covered_bases == 0:
            return
        if covered_bases == bin_size or allow_partial:
            depth = median(
                value for value, count in coverage_counts.items() for _ in range(count)
            )
            _write_row(
                output,
                current_chrom,
                current_start,
                end,
                depth,
                coordinate_system,
            )
            histogram[depth] += 1
        coverage_counts = Counter()
        covered_bases = 0

    with _open_text(input_path) as source, _open_output(output_path) as output:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or (skip_comments and line.startswith("#")):
                continue
            input_chrom = line.split("\t", 1)[0]
            if input_chrom in excluded or (
                contig is not None and input_chrom != contig
            ):
                continue
            chrom, start, end, depth = _parse_record(
                line, line_number, truncate_depth
            )

            if require_contiguous:
                if chrom != current_chrom:
                    if current_chrom is not None:
                        finish_contiguous(output, previous_end)
                    if start != 0:
                        raise ValueError(f"Chromosome {chrom} starts at {start}, not 0")
                    current_chrom = chrom
                    current_start = 0
                    current_end = bin_size
                    previous_end = 0
                elif start != previous_end:
                    raise ValueError(
                        f"Coverage intervals are not contiguous at line {line_number}: "
                        f"expected {chrom}:{previous_end}, found {chrom}:{start}"
                    )

                position = start
                while position < end:
                    overlap_end = min(end, current_end)
                    overlap = overlap_end - position
                    coverage_counts[depth] += overlap
                    covered_bases += overlap
                    position = overlap_end
                    if position == current_end:
                        finish_contiguous(output, current_end)
                        current_start = current_end
                        current_end += bin_size
                previous_end = end
                continue

            if chrom not in contigs:
                contigs.append(chrom)
            start_bin = (start // bin_size) * bin_size
            end_bin = ((end - 1) // bin_size) * bin_size
            for bin_start in range(start_bin, end_bin + 1, bin_size):
                bin_end = bin_start + bin_size
                overlap = min(end, bin_end) - max(start, bin_start)
                key = (chrom, bin_start)
                if key not in bins:
                    bins[key] = Counter()
                bins[key][depth] += overlap

        if require_contiguous:
            if current_chrom is None:
                raise ValueError("Coverage input contains no records")
            finish_contiguous(output, previous_end)
        else:
            for chrom in contigs:
                for bin_start in sorted(
                    start
                    for current_chrom_name, start in bins
                    if current_chrom_name == chrom
                ):
                    depth_counts = bins[(chrom, bin_start)]
                    covered = sum(depth_counts.values())
                    if covered != bin_size and not (allow_partial and covered > 0):
                        continue
                    depth = median(
                        value
                        for value, count in depth_counts.items()
                        for _ in range(count)
                    )
                    emitted_end = (
                        bin_start + covered
                        if allow_partial and covered < bin_size
                        else bin_start + bin_size
                    )
                    _write_row(
                        output,
                        output_contig or chrom,
                        bin_start,
                        emitted_end,
                        depth,
                        coordinate_system,
                    )
                    histogram[depth] += 1

    if histogram_path is not None:
        _write_histogram(histogram, histogram_path)
    return histogram


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bin-size", required=True, type=int)
    parser.add_argument(
        "--coordinate-system", choices=("zero-based", "one-based"), default="zero-based"
    )
    parser.add_argument("--contig")
    parser.add_argument("--output-contig")
    parser.add_argument("--allow-partial", action="store_true")
    parser.add_argument("--require-contiguous", action="store_true")
    parser.add_argument("--skip-comments", action="store_true")
    parser.add_argument("--exclude-contig", action="append", default=[])
    parser.add_argument("--truncate-depth", action="store_true")
    parser.add_argument("--histogram-output", type=Path)
    args = parser.parse_args()
    bin_mosdepth(
        args.input,
        args.output,
        args.bin_size,
        coordinate_system=args.coordinate_system,
        contig=args.contig,
        output_contig=args.output_contig,
        allow_partial=args.allow_partial,
        require_contiguous=args.require_contiguous,
        skip_comments=args.skip_comments,
        exclude_contigs=args.exclude_contig,
        truncate_depth=args.truncate_depth,
        histogram_path=args.histogram_output,
    )


if __name__ == "__main__":
    main()
