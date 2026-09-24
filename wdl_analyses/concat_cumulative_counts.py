#!/usr/bin/env python3
"""Concatenate the per-i one-row outputs of cumulative_union_counts.sh (in
scatter order, i.e. already rank-ordered) into the final cumulative table,
prepending the header.

Usage: concat_cumulative_counts.py <out.tsv> <row1.tsv> [<row2.tsv> ...]
"""
import sys

HEADER = ["sample_id", "cumulative_n_samples", "cumulative_n_total_nonref",
          "cumulative_n_snv", "cumulative_n_del_1_49", "cumulative_n_ins_1_49",
          "cumulative_n_del_50_499", "cumulative_n_ins_50_499",
          "cumulative_n_del_gt499", "cumulative_n_ins_gt499"]


def main():
    out_path, row_paths = sys.argv[1], sys.argv[2:]
    with open(out_path, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        for row_path in row_paths:
            with open(row_path) as f:
                out.write(f.readline())
    print(f"Saved: {out_path} ({len(row_paths)} rows)")


if __name__ == "__main__":
    main()
