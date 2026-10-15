#!/usr/bin/env python3
"""Merge cumulative union tables (ranked_sample_burden_table outputs of
PerSampleVariantBurden.wdl, e.g. one per contig or contig group) by summing
cumulative_n_total_nonref (3rd column) at each cumulative_n_samples.

Every *.tsv in the input folder is read (the output file itself is skipped). Files
may have all 10 columns or only the first 3; only the first 3 are used. Rows are
matched on cumulative_n_samples. Only cutoffs present in every file are written;
cutoffs missing from some files are reported and dropped, since a partial sum
would be too low. sample_id is kept when all files agree at that cutoff (same
sample ranking), else written as "." with a warning.

Usage: merge_cumulative_totals.py <input_folder> <out.tsv>
Writes <out.tsv>: sample_id, cumulative_n_samples, cumulative_n_total_nonref, n_files
"""
import argparse
import glob
import os
import sys

KEY = "cumulative_n_samples"
TOTAL = "cumulative_n_total_nonref"


def read_table(path):
    """{cumulative_n_samples: (sample_id, cumulative_n_total_nonref)} from one file."""
    rows = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header[:3] != ["sample_id", KEY, TOTAL]:
            sys.exit(f"{path}: first 3 columns must be sample_id, {KEY}, {TOTAL}; got {header[:3]}")
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            n = int(cols[1])
            if n in rows:
                sys.exit(f"{path}: {KEY}={n} appears more than once")
            rows[n] = (cols[0], int(cols[2]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input_folder", help="folder containing the *.tsv tables to merge")
    parser.add_argument("out", help="output TSV path")
    args = parser.parse_args()

    out_abs = os.path.abspath(args.out)
    paths = sorted(p for p in glob.glob(os.path.join(args.input_folder, "*.tsv"))
                   if os.path.abspath(p) != out_abs)
    if not paths:
        sys.exit(f"no *.tsv files in {args.input_folder}")
    tables = {p: read_table(p) for p in paths}
    print(f"Read {len(paths)} files", file=sys.stderr)

    all_cutoffs = set().union(*(t.keys() for t in tables.values()))
    shared = set.intersection(*(set(t) for t in tables.values()))
    dropped = sorted(all_cutoffs - shared)
    if dropped:
        print(f"WARNING: {len(dropped)} cutoffs are missing from some files and are dropped: "
              f"{dropped[:10]}{' ...' if len(dropped) > 10 else ''}", file=sys.stderr)
        for p, t in tables.items():
            missing = [n for n in dropped if n not in t]
            if missing:
                print(f"  {os.path.basename(p)} lacks {len(missing)} of them, e.g. {missing[:5]}", file=sys.stderr)

    n_sample_mismatch = 0
    with open(args.out, "w") as out:
        out.write("\t".join(["sample_id", KEY, TOTAL, "n_files"]) + "\n")
        for n in sorted(shared):
            sample_ids = {t[n][0] for t in tables.values()}
            if len(sample_ids) == 1:
                sample_id = sample_ids.pop()
            else:
                sample_id = "."
                n_sample_mismatch += 1
            total = sum(t[n][1] for t in tables.values())
            out.write(f"{sample_id}\t{n}\t{total}\t{len(tables)}\n")
    if n_sample_mismatch:
        print(f"WARNING: at {n_sample_mismatch} cutoffs the files name different samples (different "
              f"rankings); sample_id written as '.'", file=sys.stderr)
    print(f"Saved: {args.out} ({len(shared)} cutoffs, summed over {len(paths)} files)", file=sys.stderr)


if __name__ == "__main__":
    main()
