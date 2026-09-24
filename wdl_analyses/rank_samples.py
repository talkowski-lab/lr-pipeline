#!/usr/bin/env python3
"""
1. Sum per-sample n_total_nonref (from per_sample_type_counts.sh, already
   FILTER=PASS-restricted) across all chromosome-shard TSVs to get each
   sample's genome-wide individual non-ref site count.
2. Restrict to samples present in the sample/population list; report any
   mismatches. Accepts either:
     - the real ancestry file, e.g. ancestry_label.HGSVC_HPRC.tsv.gz
       (gzip, header with "sample" + "ethnicity" columns among others --
       "ethnicity" holds the SAS/EAS/EUR/AMR/AFR-style population label;
       picked by header name so extra columns like Population/color are
       ignored), or
     - a plain generic 2-column (sample_id, population) file, with or
       without a header, for backward compatibility.
3. Sort samples ascending by that individual total (ties broken by sample
   ID, for a reproducible order), and write:
     - {out_prefix}.ordered_samples.txt -- one sample ID per line, in rank
       order (input to the cumulative union-site counting step)
     - {out_prefix}.sample_rank.tsv -- sample_id, rank, n_total_nonref
       (diagnostic/audit trail for the ranking itself)

Usage:
    rank_samples.py <sample_population_list.tsv[.gz]> <out_prefix> <shard1.tsv> [<shard2.tsv> ...]
"""
import gzip
import sys
from collections import defaultdict

VALID_POPS = {"SAS", "EAS", "EUR", "AMR", "AFR"}


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    with open(path, "rb") as f:
        is_gzip = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if is_gzip else open(path)


def load_sample_pop(pop_list_path):
    """Column-name-based lookup for "sample" + "ethnicity" (as in
    ancestry_label.HGSVC_HPRC.tsv.gz); falls back to positional columns
    1 and 2 for a plain generic sample/population list."""
    sample_pop = {}
    with open_maybe_gzip(pop_list_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        lower = [c.lower() for c in header]
        if "sample" in lower and "ethnicity" in lower:
            sample_i, pop_i = lower.index("sample"), lower.index("ethnicity")
        elif "sample" in lower and "population" in lower:
            sample_i, pop_i = lower.index("sample"), lower.index("population")
        else:
            sample_i, pop_i = None, None  # no recognized header -- treat first line as data

        def parse_row(parts):
            if len(parts) <= max(sample_i, pop_i):
                return
            sample_pop[parts[sample_i]] = parts[pop_i]

        if sample_i is not None:
            for line in f:
                parse_row(line.rstrip("\n").split("\t"))
        else:
            for line in ["\t".join(header)] + [l.rstrip("\n") for l in f]:
                parts = line.split("\t")
                if len(parts) < 2 or parts[0] in ("sample", "sample_id"):
                    continue
                sample_pop[parts[0]] = parts[1]
    return sample_pop


def main():
    pop_list_path, out_prefix = sys.argv[1], sys.argv[2]
    shard_paths = sys.argv[3:]

    sample_pop = load_sample_pop(pop_list_path)
    print(f"Population list: {len(sample_pop)} samples", file=sys.stderr)
    n_unknown_pop = sum(1 for p in sample_pop.values() if p not in VALID_POPS)
    if n_unknown_pop:
        print(f"WARNING: {n_unknown_pop} samples have a population label outside "
              f"{sorted(VALID_POPS)}", file=sys.stderr)

    totals = defaultdict(int)
    for shard_path in shard_paths:
        with open(shard_path) as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = {c: i for i, c in enumerate(header)}
            for line in f:
                cols = line.rstrip("\n").split("\t")
                sample = cols[idx["sample"]]
                totals[sample] += int(cols[idx["n_total_nonref"]])
    print(f"Summed {len(shard_paths)} chromosome shards; {len(totals)} distinct samples seen in VCFs",
          file=sys.stderr)

    vcf_samples = set(totals)
    pop_samples = set(sample_pop)
    in_both = vcf_samples & pop_samples
    only_in_vcf = vcf_samples - pop_samples
    only_in_pop_list = pop_samples - vcf_samples
    print(f"Samples in both VCFs and population list (ranked below): {len(in_both)}", file=sys.stderr)
    if only_in_vcf:
        print(f"Samples in VCFs but NOT in population list (excluded): {len(only_in_vcf)}", file=sys.stderr)
    if only_in_pop_list:
        print(f"Samples in population list but NOT seen in any VCF (ignored): {len(only_in_pop_list)}",
              file=sys.stderr)

    ordered = sorted(in_both, key=lambda s: (totals[s], s))

    with open(f"{out_prefix}.ordered_samples.txt", "w") as f:
        for sample in ordered:
            f.write(sample + "\n")

    with open(f"{out_prefix}.sample_rank.tsv", "w") as f:
        f.write("sample_id\trank\tn_total_nonref\n")
        for i, sample in enumerate(ordered, start=1):
            f.write(f"{sample}\t{i}\t{totals[sample]}\n")

    print(f"Saved: {out_prefix}.ordered_samples.txt and {out_prefix}.sample_rank.tsv "
          f"({len(ordered)} samples)", file=sys.stderr)


if __name__ == "__main__":
    main()
