#!/usr/bin/env python3
"""Build the cumulative union-site table from per-contig min-rank histograms.

Each union_min_rank_histogram.sh output counts, per rank r, the PASS sites whose
best-ranked non-ref carrier has rank r. Summing those over contigs and taking the
cumulative sum over r gives, for each i, the number of distinct sites non-ref in at
least one of the top-i samples (a union, not a sum), split by type/size.

Rows are written only at checkpoint sample sizes (rank cutoffs), not at every i:
  1..10 by 1; 12..20 by 2; 25..100 by 5; 110..1,000 by 10; 1,100..10,000 by 100;
  then by 1,000; plus always the last rank N (e.g. N = 1,319: ..., 1,200, 1,300, 1,319).
All ranks are still accumulated, so each checkpoint row is exact.

Output columns (same as the earlier per-i implementation): sample_id (the i-th
ranked sample), cumulative_n_samples (i), cumulative_n_total_nonref,
cumulative_n_snv, cumulative_n_del_1_49, cumulative_n_ins_1_49,
cumulative_n_del_50_499, cumulative_n_ins_50_499, cumulative_n_del_gt499,
cumulative_n_ins_gt499.

Usage: cumulative_union_from_min_rank.py --ordered-samples ordered.txt \
           --histograms h1.tsv [h2.tsv ...] --out out.tsv
"""
import argparse

# (upper bound of the range, step within it); beyond the last bound the last step applies.
CHECKPOINT_STEPS = [(10, 1), (20, 2), (100, 5), (1000, 10), (10000, 100)]
LARGE_STEP = 1000

CATS = ["n_total_nonref", "n_snv", "n_del_1_49", "n_ins_1_49", "n_del_50_499",
        "n_ins_50_499", "n_del_gt499", "n_ins_gt499"]


def checkpoints(n):
    """Rank cutoffs at which to report the cumulative union, ending with n itself."""
    out, i, lower = [], 0, 0
    for upper, step in CHECKPOINT_STEPS + [(float("inf"), LARGE_STEP)]:
        i = lower + step
        while i <= min(upper, n):
            out.append(i)
            i += step
        lower = upper
        if upper >= n:
            break
    if not out or out[-1] != n:
        out.append(n)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ordered-samples", required=True, help="one sample ID per line, rank 1 first")
    ap.add_argument("--histograms", nargs="+", required=True, help="per-contig *.min_rank_histogram.tsv")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.ordered_samples) as f:
        samples = [line.strip() for line in f if line.strip()]
    per_rank = [[0] * len(CATS) for _ in range(len(samples) + 1)]

    for path in args.histograms:
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            if header != ["rank"] + CATS:
                raise SystemExit(f"{path}: unexpected header {header}")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                r = int(cols[0])
                if not 1 <= r <= len(samples):
                    raise SystemExit(f"{path}: rank {r} outside 1..{len(samples)}")
                for k, v in enumerate(cols[1:]):
                    per_rank[r][k] += int(v)

    report = set(checkpoints(len(samples)))
    running = [0] * len(CATS)
    with open(args.out, "w") as out:
        out.write("\t".join(["sample_id", "cumulative_n_samples"] + [f"cumulative_{c}" for c in CATS]) + "\n")
        for i, sample in enumerate(samples, start=1):
            running = [a + b for a, b in zip(running, per_rank[i])]
            if i in report:
                out.write("\t".join([sample, str(i)] + [str(v) for v in running]) + "\n")
    print(f"Saved: {args.out} ({len(report)} checkpoint rows of {len(samples)} ranks, from "
          f"{len(args.histograms)} contig histograms; union of all ranked samples = {running[0]} PASS sites)")


if __name__ == "__main__":
    main()
