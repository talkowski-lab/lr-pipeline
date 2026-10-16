#!/usr/bin/env python3
"""Sum per-sample variant counts (from per_sample_type_counts.sh) across
chromosome-shard TSVs into one genome-wide flat table: one row per sample,
one column per type/size/REGION category. No ranking, no cumulative sums,
no population-list filtering -- every sample seen in any shard is included.
"""
import argparse
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shards", nargs="+", required=True,
                         help="per_sample_counts.tsv files, one per chromosome shard")
    parser.add_argument("--out", required=True, help="output TSV path")
    return parser.parse_args()


def main():
    args = parse_args()

    cat_order = []
    cat_set = set()
    totals = defaultdict(lambda: defaultdict(int))

    for shard_path in args.shards:
        with open(shard_path) as f:
            header = f.readline().rstrip("\n").split("\t")
            cats = header[1:]
            for c in cats:
                if c not in cat_set:
                    cat_set.add(c)
                    cat_order.append(c)
            idx = {c: i for i, c in enumerate(header)}
            for line in f:
                cols = line.rstrip("\n").split("\t")
                sample = cols[idx["sample"]]
                for c in cats:
                    totals[sample][c] += int(cols[idx[c]])

    samples = sorted(totals)
    with open(args.out, "w") as out:
        out.write("\t".join(["sample"] + cat_order) + "\n")
        for sample in samples:
            row = [sample] + [str(totals[sample][c]) for c in cat_order]
            out.write("\t".join(row) + "\n")

    print(f"Summed {len(args.shards)} chromosome shards; "
          f"wrote {args.out} ({len(samples)} samples, {len(cat_order)} count columns)")


if __name__ == "__main__":
    main()
