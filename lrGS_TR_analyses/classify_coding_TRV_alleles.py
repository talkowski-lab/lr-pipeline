#!/usr/bin/env python3
"""Classify ALT alleles of coding TRV sites by their length change vs REF (delta = len(ALT) - len(REF)).

Only PASS sites with genic_context == coding (from annotate_TRV_genic_context.sh) are used, and only alleles with AC > 0.
A site is "CDS-contained" if its locus [start, end) lies entirely inside one CDS block of the GTF; otherwise it
straddles a CDS boundary and its alleles are labelled partial_CDS (the length change may not be in coding sequence).

Allele classes (CDS-contained sites):
  frameshift      delta % 3 != 0 (putative pLoF; no NMD / last-exon rules applied)
  inframe_del     delta < 0 and delta % 3 == 0
  inframe_ins     delta > 0 and delta % 3 == 0
  no_len_change   delta == 0 (sequence-only change; not classified further)

Outputs:
  <out-prefix>.alleles.tsv.gz  one row per allele: ID, coords, genes, min_motif_len, cds_contained, delta, AC, class
  <out-prefix>.summary.tsv     per AC filter x class: n alleles, n sites (>=1 allele of class), n genes (multi-gene: all)

Usage:
  classify_coding_TRV_alleles.py --genic-tsv genic.tsv.gz --trv-bed trv.bed.gz --gtf gencode.gtf.gz \
      --label AoU_I --out-prefix out/x
"""
import argparse
import gzip

import numpy as np
import pandas as pd

CLASSES = ["frameshift", "inframe_del", "inframe_ins", "no_len_change", "partial_CDS"]


def cds_blocks(gtf):
    """Per-chrom sorted, merged CDS blocks (0-based half-open) as (starts, ends) arrays."""
    rows = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split("\t", 5)
            if t[2] == "CDS":
                rows.append((t[0], int(t[3]) - 1, int(t[4])))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end"]).sort_values(["chrom", "start"])
    blocks = {}
    for chrom, sub in df.groupby("chrom"):
        merged = []
        for s, e in zip(sub["start"], sub["end"]):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        m = np.array(merged)
        blocks[chrom] = (m[:, 0], m[:, 1])
    return blocks


def contained(chroms, starts, ends, blocks):
    out = np.zeros(len(starts), dtype=bool)
    for chrom in np.unique(chroms):
        idx = np.where(chroms == chrom)[0]
        if chrom not in blocks:
            continue
        bs, be = blocks[chrom]
        k = np.searchsorted(bs, starts[idx], side="right") - 1
        ok = k >= 0
        out[idx[ok]] = ends[idx[ok]] <= be[k[ok]]
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genic-tsv", required=True)
    p.add_argument("--trv-bed", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    g = pd.read_csv(args.genic_tsv, sep="\t", dtype=str,
                    usecols=["chrom", "start", "end", "ID", "FILTER", "min_motif_len", "alt_len_diffs", "genic_context", "genes"])
    g = g[(g["FILTER"] == "PASS") & (g["genic_context"] == "coding")]
    ac = pd.concat(ch[ch["ID"].isin(set(g["ID"]))] for ch in pd.read_csv(
        args.trv_bed, sep="\t", usecols=["ID", "AC"], dtype=str, chunksize=500000))
    g = g.merge(ac, on="ID", how="left")
    assert g["AC"].notna().all()

    starts = g["start"].astype(int).to_numpy()
    ends = g["end"].astype(int).to_numpy()
    g["cds_contained"] = contained(g["chrom"].to_numpy(), starts, ends, cds_blocks(args.gtf))

    g["delta"] = g["alt_len_diffs"].str.split(",")
    g["AC"] = g["AC"].str.split(",")
    assert (g["delta"].str.len() == g["AC"].str.len()).all()
    a = g.drop(columns=["alt_len_diffs", "FILTER", "genic_context"]).explode(["delta", "AC"])
    a["delta"] = a["delta"].astype(int)
    a["AC"] = a["AC"].astype(int)
    a = a[a["AC"] > 0]

    d = a["delta"].to_numpy()
    cls = np.select([d == 0, d % 3 != 0, d < 0], ["no_len_change", "frameshift", "inframe_del"], "inframe_ins")
    a["class"] = np.where(a["cds_contained"], cls, "partial_CDS")
    a.to_csv(f"{args.out_prefix}.alleles.tsv.gz", sep="\t", index=False, compression="gzip")

    rows = []
    for filt, m in [("AC>=1", a["AC"] >= 1), ("AC>=2", a["AC"] >= 2), ("AC>=5", a["AC"] >= 5)]:
        sub = a[m]
        for c in CLASSES:
            s = sub[sub["class"] == c]
            genes = {x for gs in s["genes"].unique() for x in gs.split(",")}
            rows.append({"cohort": args.label, "allele_filter": filt, "class": c, "n_alleles": len(s),
                         "n_sites": s["ID"].nunique(), "n_genes": len(genes)})
        genes_any = {x for gs in sub.loc[sub["class"] != "partial_CDS", "genes"].unique() for x in gs.split(",")}
        rows.append({"cohort": args.label, "allele_filter": filt, "class": "any_CDS_contained",
                     "n_alleles": int((sub["class"] != "partial_CDS").sum()),
                     "n_sites": sub.loc[sub["class"] != "partial_CDS", "ID"].nunique(), "n_genes": len(genes_any)})
    out = pd.DataFrame(rows)
    out.to_csv(f"{args.out_prefix}.summary.tsv", sep="\t", index=False)
    n_sites = g["ID"].nunique()
    print(f"{args.label}: {n_sites:,} PASS coding sites; CDS-contained {int(g['cds_contained'].sum()):,}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
