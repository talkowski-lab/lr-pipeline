#!/usr/bin/env python3
"""Concatenate the per-i one-row outputs of cumulative_union_counts.sh (in
scatter order, i.e. already rank-ordered) into the final cumulative table,
prepending the header. Rows have 10 columns when split by type/size
(BY_ALLELE_TYPE=true) or 3 (sample_id, i, total) when not; the header is
chosen to match, and rows of mixed widths are an error.

Usage: concat_cumulative_counts.py <out.tsv> <row1.tsv> [<row2.tsv> ...]
"""
import sys

HEADER = ["sample_id", "cumulative_n_samples", "cumulative_n_total_nonref",
          "cumulative_n_snv", "cumulative_n_del_1_49", "cumulative_n_ins_1_49",
          "cumulative_n_del_50_499", "cumulative_n_ins_50_499",
          "cumulative_n_del_gt499", "cumulative_n_ins_gt499"]


def main():
    out_path, row_paths = sys.argv[1], sys.argv[2:]
    rows = []
    for row_path in row_paths:
        with open(row_path) as f:
            rows.append(f.readline())
    widths = sorted({len(r.rstrip("\n").split("\t")) for r in rows})
    if widths not in ([3], [len(HEADER)]):
        sys.exit(f"rows must all have 3 or {len(HEADER)} columns; found widths {widths}")
    header = HEADER if widths == [len(HEADER)] else HEADER[:3]
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        out.writelines(rows)
    print(f"Saved: {out_path} ({len(row_paths)} rows)")


if __name__ == "__main__":
    main()
