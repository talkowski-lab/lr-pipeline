#!/usr/bin/env python3
"""Per-gene counts of coding STR sites, split by LoF (frameshift) vs in-frame.

Input: <cohort>.coding_TRV_alleles.alleles.tsv.gz from classify_coding_TRV_alleles.py (one row per ALT allele with AC > 0
at PASS coding TRV sites; columns chrom, start, end, ID, min_motif_len, genes, AC, cds_contained, delta, class).

Site filters: min_motif_len <= --max-motif (STR; default 6) and locus entirely inside one CDS block (cds_contained).
Site class (from its ALT alleles with AC >= --min-ac):
  LoF            >= 1 frameshift allele (len change not a multiple of 3; putative, no NMD / last-exon rule)
  inframe        no frameshift, >= 1 in-frame deletion or insertion allele
  no_len_change  only alleles with no length change
Sites with no allele passing --min-ac are dropped. Multi-gene sites are attributed to every listed gene.

Outputs (prefix --out-prefix):
  .sites.tsv     one row per site: ID, coords, genes, motif length, site class, per-class allele counts and summed AC
  .per_gene.tsv  one row per gene: site counts by class, allele counts, summed AC, motif-3 sites
  .summary.tsv   totals: sites and genes by class

Usage:
  summarize_coding_STR_per_gene.py --alleles hprc_hgsvc.coding_TRV_alleles.alleles.tsv.gz --out-prefix out/hprc_hgsvc.coding_STR \
      [--max-motif 6] [--min-ac 1]
"""
import argparse

import pandas as pd

ALLELE_CLASSES = ["frameshift", "inframe_del", "inframe_ins", "no_len_change"]
SITE_CLASSES = ["LoF", "inframe", "no_len_change"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--alleles", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--max-motif", type=int, default=6, help="Max shortest-motif length for an STR (default 6)")
    p.add_argument("--min-ac", type=int, default=1, help="Min AC for an allele to count (default 1)")
    args = p.parse_args()

    a = pd.read_csv(args.alleles, sep="\t", dtype={"genes": str})
    n_all_sites = a["ID"].nunique()
    a = a[a["cds_contained"] & (a["AC"] >= args.min_ac)]
    n_vntr_sites = a.loc[a["min_motif_len"] > args.max_motif, "ID"].nunique()
    a = a[a["min_motif_len"] <= args.max_motif]

    for c in ALLELE_CLASSES:
        a[f"n_{c}_alleles"] = (a["class"] == c).astype(int)
        a[f"AC_{c}"] = a["AC"].where(a["class"] == c, 0)
    agg = {f"n_{c}_alleles": "sum" for c in ALLELE_CLASSES}
    agg.update({f"AC_{c}": "sum" for c in ALLELE_CLASSES})
    sites = a.groupby(["ID", "chrom", "start", "end", "genes", "min_motif_len"], as_index=False).agg(agg)
    sites["site_class"] = "no_len_change"
    sites.loc[sites["n_inframe_del_alleles"] + sites["n_inframe_ins_alleles"] > 0, "site_class"] = "inframe"
    sites.loc[sites["n_frameshift_alleles"] > 0, "site_class"] = "LoF"
    sites = sites.sort_values(["chrom", "start"])
    sites.to_csv(f"{args.out_prefix}.sites.tsv", sep="\t", index=False)

    g = sites.assign(gene=sites["genes"].str.split(",")).explode("gene")
    per_gene = g.groupby("gene").agg(n_coding_STR_sites=("ID", "size"),
                                     n_motif3_sites=("min_motif_len", lambda x: int((x == 3).sum())))
    for c in SITE_CLASSES:
        per_gene[f"n_{c}_sites"] = g[g["site_class"] == c].groupby("gene").size()
    for c in ALLELE_CLASSES:
        per_gene[f"n_{c}_alleles"] = g.groupby("gene")[f"n_{c}_alleles"].sum()
        per_gene[f"AC_{c}"] = g.groupby("gene")[f"AC_{c}"].sum()
    per_gene = per_gene.fillna(0).astype(int).sort_values(["n_coding_STR_sites", "n_LoF_sites"], ascending=False)
    per_gene.to_csv(f"{args.out_prefix}.per_gene.tsv", sep="\t")

    rows = [("coding TRV sites in input (any motif, incl. partial CDS)", n_all_sites),
            (f"CDS-contained sites with motif > {args.max_motif} bp (excluded, VNTR)", n_vntr_sites),
            ("coding STR sites", len(sites)),
            ("genes with coding STR", per_gene.shape[0])]
    for c in SITE_CLASSES:
        rows.append((f"{c} sites", int((sites["site_class"] == c).sum())))
        rows.append((f"genes with >=1 {c} site", int((per_gene[f"n_{c}_sites"] > 0).sum())))
    both = (per_gene["n_LoF_sites"] > 0) & (per_gene["n_inframe_sites"] > 0)
    rows.append(("genes with LoF and in-frame sites", int(both.sum())))
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(f"{args.out_prefix}.summary.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows, columns=["metric", "value"]).to_string(index=False))


if __name__ == "__main__":
    main()
