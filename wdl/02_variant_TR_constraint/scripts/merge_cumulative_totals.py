#!/usr/bin/env python3
"""Merge cumulative union tables (ranked_sample_burden_table outputs of
PerSampleVariantBurden.wdl, e.g. one per contig) by summing every count column
(3rd column onward: cumulative_n_total_nonref, cumulative_n_snv, ...) at each
cumulative_n_samples.

Every file in the input folder ending with --suffix is read. The output file and any
earlier merged output (a file with an n_files column) are skipped, so rerunning
never double counts. Rows are matched on cumulative_n_samples.
- Count columns: those present in every file are summed, in the first file's order.
  Columns missing from some files (e.g. 3-column files mixed with 10-column ones)
  are reported and dropped, since a partial sum would be too low.
- Cutoffs: only cutoffs present in every file are written; the others are reported
  and dropped, for the same reason.
- sample_id is kept when all files agree at that cutoff (same sample ranking), else
  written as "." with a warning.

Usage: merge_cumulative_totals.py <input_folder> <out.tsv> [--suffix .ranked_sample_burden.tsv]
Writes <out.tsv>: sample_id, cumulative_n_samples, <summed count columns...>, n_files
"""
import argparse
import os
import sys

KEY = "cumulative_n_samples"


def read_table(path):
    """(count column names, {cumulative_n_samples: (sample_id, {column: count})}) from one file."""
    rows = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:2] != ["sample_id", KEY] or len(header) < 3:
            sys.exit(f"{path}: header must start with sample_id, {KEY}, then count columns; got {header}")
        count_cols = header[2:]
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != len(header):
                sys.exit(f"{path}: row has {len(cols)} columns, header has {len(header)}: {line.rstrip()}")
            n = int(cols[1])
            if n in rows:
                sys.exit(f"{path}: {KEY}={n} appears more than once")
            rows[n] = (cols[0], dict(zip(count_cols, map(int, cols[2:]))))
    return count_cols, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_folder", help="folder containing the tables to merge")
    parser.add_argument("out", help="output TSV path")
    parser.add_argument("--suffix", default=".ranked_sample_burden.tsv",
                        help="merge files ending with this (default: .ranked_sample_burden.tsv)")
    args = parser.parse_args()

    out_abs = os.path.abspath(args.out)
    paths = sorted(os.path.join(args.input_folder, f) for f in os.listdir(args.input_folder)
                   if f.endswith(args.suffix) and os.path.abspath(os.path.join(args.input_folder, f)) != out_abs)
    if not paths:
        sys.exit(f"no files ending with {args.suffix} in {args.input_folder}")
    tables = {}
    for p in paths:
        cols, rows = read_table(p)
        if "n_files" in cols:
            print(f"Skipping {os.path.basename(p)}: already a merged output (has n_files)", file=sys.stderr)
            continue
        tables[p] = (cols, rows)
    if not tables:
        sys.exit(f"no unmerged files ending with {args.suffix} in {args.input_folder}")
    paths = list(tables)
    print(f"Read {len(paths)} files ending with {args.suffix}", file=sys.stderr)

    first_cols = tables[paths[0]][0]
    shared_cols = [c for c in first_cols if all(c in cols for cols, _ in tables.values())]
    all_cols = {c for cols, _ in tables.values() for c in cols}
    dropped_cols = sorted(all_cols - set(shared_cols))
    if dropped_cols:
        print(f"WARNING: columns missing from some files are dropped: {dropped_cols}", file=sys.stderr)
        for p, (cols, _) in tables.items():
            missing = [c for c in dropped_cols if c not in cols]
            if missing:
                print(f"  {os.path.basename(p)} lacks {missing}", file=sys.stderr)

    all_cutoffs = {n for _, rows in tables.values() for n in rows}
    shared = set.intersection(*(set(rows) for _, rows in tables.values()))
    dropped = sorted(all_cutoffs - shared)
    if dropped:
        print(f"WARNING: {len(dropped)} cutoffs are missing from some files and are dropped: "
              f"{dropped[:10]}{' ...' if len(dropped) > 10 else ''}", file=sys.stderr)
        for p, (_, rows) in tables.items():
            missing = [n for n in dropped if n not in rows]
            if missing:
                print(f"  {os.path.basename(p)} lacks {len(missing)} of them, e.g. {missing[:5]}", file=sys.stderr)

    n_sample_mismatch = 0
    with open(args.out, "w") as out:
        out.write("\t".join(["sample_id", KEY] + shared_cols + ["n_files"]) + "\n")
        for n in sorted(shared):
            sample_ids = {rows[n][0] for _, rows in tables.values()}
            if len(sample_ids) == 1:
                sample_id = sample_ids.pop()
            else:
                sample_id = "."
                n_sample_mismatch += 1
            sums = [sum(rows[n][1][c] for _, rows in tables.values()) for c in shared_cols]
            out.write("\t".join([sample_id, str(n)] + list(map(str, sums)) + [str(len(tables))]) + "\n")
    if n_sample_mismatch:
        print(f"WARNING: at {n_sample_mismatch} cutoffs the files name different samples (different "
              f"rankings); sample_id written as '.'", file=sys.stderr)
    print(f"Saved: {args.out} ({len(shared)} cutoffs x {len(shared_cols)} count columns, "
          f"summed over {len(paths)} files)", file=sys.stderr)


if __name__ == "__main__":
    main()
