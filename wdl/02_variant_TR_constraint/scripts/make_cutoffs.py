#!/usr/bin/env python3
"""Checkpoint sample sizes (rank cutoffs) at which to compute the cumulative union.

Cutoffs: 1..10 by 1; 12..20 by 2; 25..100 by 5; 110..1,000 by 10; 1,100..10,000 by
100; then by 1,000; plus always the last rank N (e.g. N = 1,319: ..., 1,200, 1,300,
1,319). Each cutoff becomes one CumulativeUnionCounts task.

Usage: make_cutoffs.py <ordered_samples.txt> <cutoffs.txt>
  ordered_samples.txt: one sample ID per line, rank 1 first (RankSamples output)
  cutoffs.txt: one cutoff per line, ascending
"""
import sys

# (upper bound of the range, step within it); beyond the last bound LARGE_STEP applies.
CHECKPOINT_STEPS = [(10, 1), (20, 2), (100, 5), (1000, 10), (10000, 100)]
LARGE_STEP = 1000


def checkpoints(n):
    """Rank cutoffs at which to report the cumulative union, ending with n itself."""
    out, lower = [], 0
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
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    with open(sys.argv[1]) as f:
        n = sum(1 for line in f if line.strip())
    if n == 0:
        sys.exit(f"{sys.argv[1]} has no samples")
    cutoffs = checkpoints(n)
    with open(sys.argv[2], "w") as out:
        out.write("\n".join(map(str, cutoffs)) + "\n")
    print(f"{len(cutoffs)} cutoffs for {n} ranked samples: {cutoffs[:12]} ... {cutoffs[-4:]}")


if __name__ == "__main__":
    main()
