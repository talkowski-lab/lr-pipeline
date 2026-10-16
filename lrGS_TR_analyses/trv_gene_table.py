#!/usr/bin/env python3
"""Gene table of STR counts by genic context from the merged TRV table.

Input: merged TRV BED (merge_trv_tables.sh output; columns used: FILTER, motifs, AC, genic_context, genes, dist_closest_5UTR,
gene_closest_5UTR, dist_closest_3UTR, gene_closest_3UTR, LOEUF_v4, LOEUF_v4_percentile).
Sites: FILTER PASS, summed AC > 0, shortest motif length in [--min-motif, --max-motif] (default 1-6 bp = STR).
Per gene:
  n_coding / n_UTR / n_intronic   STRs whose genes column lists the gene (multi-gene STRs count toward every listed gene)
  n_intergenic_nearest            intergenic STRs whose nearest 5'UTR or 3'UTR belongs to the gene (smaller of the two
                                  distances; both genes on a tie), split into _5UTR_side and _3UTR_side
  n_total                         sum of the above
  median_dist_intergenic          median distance (bp) of the gene's assigned intergenic STRs to that UTR
  LOEUF_v4, LOEUF_v4_percentile   from --gene-info when given (gene_name match)
Output: TSV sorted by n_total, and a printed summary.

Usage:
  trv_gene_table.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --gene-info gene_information_loeufv4...txt.gz --out genes.tsv
"""
import argparse
import collections
import gzip

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--gene-info", help="gene table with gene_name, LOEUF_v4, LOEUF_v4_percentile")
    p.add_argument("--min-motif", type=int, default=1)
    p.add_argument("--max-motif", type=int, default=6)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    counts = collections.defaultdict(collections.Counter)
    dists = collections.defaultdict(list)
    n_used = collections.Counter()
    with gzip.open(args.trv, "rt") as f:
        h = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
        for line in f:
            t = line.rstrip("\n").split("\t")
            if t[h["FILTER"]] != "PASS":
                continue
            if sum(int(x) for x in t[h["AC"]].split(",") if x not in (".", "")) <= 0:
                continue
            m = min(len(x) for x in t[h["motifs"]].split(","))
            if not args.min_motif <= m <= args.max_motif:
                continue
            ctx = t[h["genic_context"]]
            n_used[ctx] += 1
            if ctx in ("coding", "UTR", "intronic"):
                for g in t[h["genes"]].split(","):
                    counts[g][f"n_{ctx}"] += 1
            elif ctx == "intergenic":
                d5, g5 = t[h["dist_closest_5UTR"]], t[h["gene_closest_5UTR"]]
                d3, g3 = t[h["dist_closest_3UTR"]], t[h["gene_closest_3UTR"]]
                d5 = float(d5) if d5 not in (".", "") else np.inf
                d3 = float(d3) if d3 not in (".", "") else np.inf
                if not np.isfinite(min(d5, d3)):
                    n_used["intergenic_unassigned"] += 1
                    continue
                if d5 < d3:
                    sides = [("5UTR", d5, g5)]
                elif d3 < d5:
                    sides = [("3UTR", d3, g3)]
                else:
                    sides = [("5UTR", d5, g5), ("3UTR", d3, g3)]
                for side, d, gs in sides:
                    for g in gs.split(","):
                        counts[g]["n_intergenic_nearest"] += 1
                        counts[g][f"n_intergenic_nearest_{side}_side"] += 1
                        dists[g].append(d)

    cols = ["n_coding", "n_UTR", "n_intronic", "n_intergenic_nearest", "n_intergenic_nearest_5UTR_side",
            "n_intergenic_nearest_3UTR_side"]
    rows = []
    for g, c in counts.items():
        row = {"gene": g, **{k: c.get(k, 0) for k in cols}}
        row["n_total"] = row["n_coding"] + row["n_UTR"] + row["n_intronic"] + row["n_intergenic_nearest"]
        row["median_dist_intergenic"] = float(np.median(dists[g])) if dists[g] else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    if args.gene_info:
        gi = pd.read_csv(args.gene_info, sep="\t", usecols=["gene_name", "LOEUF_v4", "LOEUF_v4_percentile"], na_values=["NA"])
        gi = gi.drop_duplicates("gene_name").rename(columns={"gene_name": "gene"})
        out = out.merge(gi, on="gene", how="left")
    out = out.sort_values(["n_total", "n_coding"], ascending=False)
    out.to_csv(args.out, sep="\t", index=False, float_format="%.6g", na_rep="NA")

    print(f"STRs used (motif {args.min_motif}-{args.max_motif} bp, PASS, AC > 0): " +
          ", ".join(f"{k} {v:,}" for k, v in n_used.items()))
    has = {k: int((out[k] > 0).sum()) for k in ["n_coding", "n_UTR", "n_intronic", "n_intergenic_nearest"]}
    print(f"genes: {len(out):,}; with >= 1 STR by context: " + ", ".join(f"{k[2:]} {v:,}" for k, v in has.items()))
    if "LOEUF_v4" in out:
        print(f"genes with LOEUF_v4: {int(out.LOEUF_v4.notna().sum()):,}")
    print(out.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
