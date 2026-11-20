#!/usr/bin/env python3
"""Combine per-contig per_sample_category_counts.py outputs into one genome-wide
table: sum the count columns per sample across contigs, and take the UNION
(not concatenation) of each sample's gene-set columns across contigs before
recomputing the unique-gene-count columns, since the same sample's per-contig
gene sets must be merged as sets to avoid inflating counts if a symbol were
ever repeated (not expected across distinct chromosomes, but still the
correct operation).

Usage: concat_sample_category_counts.py <out.tsv> <shard1.tsv> <shard2.tsv> ...
"""
import sys
from collections import defaultdict

COUNT_CATEGORIES = [
    "n_plof_snv", "n_plof_indel_del", "n_plof_indel_ins", "n_plof_sv_del", "n_plof_sv_ins",
    "n_missense", "n_synonymous", "n_intronic", "n_intergenic",
]
GENE_CATEGORIES = ["plof_snv_genes", "plof_indel_genes", "plof_sv_genes"]


def main():
    out_path, shard_paths = sys.argv[1], sys.argv[2:]

    counts = defaultdict(lambda: defaultdict(int))
    genes = defaultdict(lambda: {g: set() for g in GENE_CATEGORIES})
    sample_order = []
    seen = set()

    for path in shard_paths:
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            col = {name: i for i, name in enumerate(header)}
            for line in f:
                fields = line.rstrip("\n").split("\t")
                sample = fields[col["sample"]]
                if sample not in seen:
                    seen.add(sample)
                    sample_order.append(sample)
                for c in COUNT_CATEGORIES:
                    counts[sample][c] += int(fields[col[c]])
                for g in GENE_CATEGORIES:
                    val = fields[col[g]]
                    if val:
                        genes[sample][g].update(val.split(","))

    header = (
        ["sample"] + COUNT_CATEGORIES + GENE_CATEGORIES
        + ["n_unique_genes_plof_snv", "n_unique_genes_plof_indel", "n_unique_genes_plof_sv", "n_unique_genes_plof_any"]
    )
    with open(out_path, "w") as out:
        out.write("\t".join(header) + "\n")
        for s in sample_order:
            row = [s] + [str(counts[s][c]) for c in COUNT_CATEGORIES]
            row += [",".join(sorted(genes[s][g])) for g in GENE_CATEGORIES]
            row += [str(len(genes[s][g])) for g in GENE_CATEGORIES]
            row.append(str(len(set.union(*genes[s].values()))))
            out.write("\t".join(row) + "\n")

    print(f"Combined {len(shard_paths)} shards across {len(sample_order)} samples; wrote {out_path}")


if __name__ == "__main__":
    main()
