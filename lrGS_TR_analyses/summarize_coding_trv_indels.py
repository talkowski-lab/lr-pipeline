#!/usr/bin/env python3
"""Coding TRVs and the indels inside them: share of indels not captured by short-read resources (srGS).

Inputs:
  --indels  from extract_indels_in_trvs.sh (TRV_ID, variant_ID, allele_type, abs_length, srGS_captured, ...) run on the
            coding TRV IDs (non-ref PASS TRVs whose TR span lies inside one CDS block)
  --loci    per-locus table from build_trv_locus_table.py (ID, size, genic_context, FILTER, n_alt)
Not captured by srGS = "." in both dbGaP_ID and gnomAD_V4_match_ID.
Two tables, both with bins 1-49, 50-149, 150-249, >=250 bp and an "all" row, split by all / del / ins:
  by_indel_size  bin = |allele_length| of the indel
  by_TRV_size    bin = TRV reference size (TRID span); also counts coding TRVs in the bin and those with >= 1 indel
Output: <out-prefix>.tsv (both tables, long format) and the tables printed.

Usage:
  summarize_coding_trv_indels.py --indels hprc_hgsvc.coding_TRV_indels.tsv.gz --loci hprc_hgsvc.TRV_loci.tsv.gz \
      --label hprc_hgsvc --out-prefix out/hprc_hgsvc.coding_TRV_indels
"""
import argparse

import numpy as np
import pandas as pd

EDGES = [1, 50, 150, 250, np.inf]
LABELS = ["1–49 bp", "50–149 bp", "150–249 bp", "≥250 bp"]


def summarize(ind, bin_col):
    rows = []
    for typ, s_t in [("all", ind), ("del", ind[ind["allele_type"] == "del"]), ("ins", ind[ind["allele_type"] == "ins"])]:
        for lab in ["all"] + LABELS:
            s = s_t if lab == "all" else s_t[s_t[bin_col] == lab]
            n = len(s)
            nc = int((s["srGS_captured"] == 0).sum())
            rows.append({"indel_type": typ, "bin": lab, "n_indels": n, "n_not_captured_srGS": nc,
                         "pct_not_captured_srGS": round(100 * nc / n, 1) if n else np.nan,
                         "n_TRVs_with_indel": int(s["TRV_ID"].nunique())})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--indels", required=True)
    p.add_argument("--loci", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    loci = pd.read_csv(args.loci, sep="\t", usecols=["ID", "FILTER", "n_alt", "genic_context", "size"])
    coding = loci[(loci["genic_context"] == "coding") & (loci["FILTER"] == "PASS") & (loci["n_alt"] > 0)].copy()
    coding["trv_bin"] = pd.cut(coding["size"], EDGES, right=False, labels=LABELS).astype(str)
    ind = pd.read_csv(args.indels, sep="\t", dtype={"dbGaP_ID": str, "gnomAD_V4_match_ID": str})
    assert ind["TRV_ID"].isin(coding["ID"]).all()
    ind = ind.merge(coding[["ID", "size", "trv_bin"]].rename(columns={"ID": "TRV_ID", "size": "trv_size"}), on="TRV_ID")
    ind["indel_bin"] = pd.cut(ind["abs_length"].clip(lower=1), EDGES, right=False, labels=LABELS).astype(str)

    print(f"{args.label}: coding TRVs (non-ref PASS, inside one CDS block): {len(coding):,}; "
          f"with >= 1 PASS indel: {ind['TRV_ID'].nunique():,}; indels: {len(ind):,} "
          f"(not captured by srGS: {int((ind['srGS_captured'] == 0).sum()):,}, "
          f"{100 * (ind['srGS_captured'] == 0).mean():.1f}%)")

    t1 = summarize(ind, "indel_bin").assign(table="by_indel_size")
    t2 = summarize(ind, "trv_bin").assign(table="by_TRV_size")
    n_trv = coding["trv_bin"].value_counts()
    t2["n_coding_TRVs_in_bin"] = t2["bin"].map(lambda b: len(coding) if b == "all" else int(n_trv.get(b, 0)))
    for name, t in [("1. Binned by indel size (|allele_length|)", t1), ("2. Binned by TRV reference size (TRID span)", t2)]:
        print(f"\n{name}")
        cols = ["indel_type", "bin"] + (["n_coding_TRVs_in_bin"] if "n_coding_TRVs_in_bin" in t else []) + \
               ["n_TRVs_with_indel", "n_indels", "n_not_captured_srGS", "pct_not_captured_srGS"]
        print(t[cols].to_string(index=False))
    pd.concat([t1, t2]).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
