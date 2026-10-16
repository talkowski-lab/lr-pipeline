#!/usr/bin/env python3
"""Coding STR variability vs copy number for one motif length, split by gene LOEUF tertile.

Inputs: TRV bed (ID, FILTER, SOURCE, TRID, REF, ALT, AC, AN), the genic-context table from annotate_TRV_genic_context.sh
(ID, genic_context, genes, alt_len_diffs), a GENCODE GTF (CDS blocks) and the gnomAD constraint table.
Loci: FILTER PASS, single-component TRID, motif length == --motif, SOURCE == --source, copy number >= --min-copies.
Coding loci must lie entirely inside one merged CDS block. Intergenic loci of the same motif are kept as a reference.

LOEUF: --loeuf-column (default lof.oe_ci.upper) from rows with mane_select == true; genes matched by symbol to the
GTF gene_name in the genes column. Tertiles are computed over all MANE genes with LOEUF (genome-wide cutoffs);
T1 = lowest LOEUF (most constrained). Multi-gene loci are attributed to each gene's tertile.

Per locus metrics (ALT alleles with AC > 0):
  nonref_AF       sum(AC) / AN                                    ("mutation rate")
  n_alt           number of distinct ALT alleles                  ("overall ALT count")
  max_abs_delta   max |len(ALT) - len(REF)| in bp, 0 if no ALT    ("max difference ALT vs REF")

Outputs (prefix --out-prefix):
  .pdf   three panels vs copy number; lines = LOEUF tertiles, plus intergenic reference. nonref_AF and n_alt: mean +-95% CI;
         max_abs_delta: median with IQR band (heavy tail: a few loci carry kb-scale SV-sized ALT alleles; mean is in the TSV)
  .tsv   per group x copy-number bin: n_loci, mean/median/sd of each metric
  .genes.tsv  coding loci with gene, LOEUF and tertile

Usage:
  plot_STR_constraint_by_loeuf.py --trv-bed TRV.bed.gz --genic-tsv genic.tsv.gz --gtf gencode.gtf.gz \
      --constraint gnomad.v4.1.1.constraint_metrics.tsv.bgz --motif 3 --label hprc_hgsvc --out-prefix out/x
"""
import argparse

import matplotlib
import numpy as np
import pandas as pd

from plot_STR_constraint_allele_class import cds_blocks, cds_contained, style

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#0b0b0b"
GROUPS = ["T1 (most constrained)", "T2", "T3 (least constrained)", "intergenic"]
GROUP_COLORS = {"T1 (most constrained)": "#0d366b", "T2": "#2a78d6", "T3 (least constrained)": "#8fbcf0",
                "intergenic": "#8a8984"}
METRICS = [("nonref_AF", "Mean non-ref allele frequency (mutation rate)"),
           ("n_alt", "Mean distinct ALT alleles per locus"),
           ("max_abs_delta", "Median max |ALT − REF| length (bp; IQR band)")]
COPY_EDGES = [2, 3, 4, 6, 9, np.inf]


def load_loci(args):
    parts = []
    for ch in pd.read_csv(args.trv_bed, sep="\t", usecols=["#CHROM", "ID", "FILTER", "SOURCE", "TRID", "AC", "AN"], dtype=str,
                          chunksize=500000):
        ch = ch[(ch["FILTER"] == "PASS") & (ch["SOURCE"] == args.source) & ~ch["TRID"].str.contains(",")]
        t = ch["TRID"].str.split("-", expand=True)
        keep = t[3].str.len() == args.motif
        ch, t = ch[keep], t[keep]
        parts.append(pd.DataFrame({"chrom": ch["#CHROM"].to_numpy(), "ID": ch["ID"].to_numpy(), "AC": ch["AC"].to_numpy(),
                                   "AN": ch["AN"].astype(int).to_numpy(), "start": t[1].astype(int).to_numpy(),
                                   "end": t[2].astype(int).to_numpy()}))
    d = pd.concat(parts, ignore_index=True)
    d["copies"] = (d["end"] - d["start"]) / args.motif
    d = d[d["copies"] >= args.min_copies]
    ctx = pd.read_csv(args.genic_tsv, sep="\t", usecols=["ID", "genic_context", "genes", "alt_len_diffs"], dtype=str)
    d = d.merge(ctx, on="ID", how="left")
    d = d[d["genic_context"].isin(["coding", "intergenic"])].reset_index(drop=True)
    contained = cds_contained(d["chrom"].to_numpy(), d["start"].to_numpy(), d["end"].to_numpy(), cds_blocks(args.gtf))
    d = d[(d["genic_context"] == "intergenic") | contained].reset_index(drop=True)

    ac = d["AC"].str.split(",").map(lambda x: np.array(x, dtype=int))
    delta = d["alt_len_diffs"].str.split(",").map(lambda x: np.abs(np.array(x, dtype=int)))
    d["n_alt"] = [int((a > 0).sum()) for a in ac]
    d["nonref_AF"] = np.array([a.sum() for a in ac]) / d["AN"].to_numpy()
    d["max_abs_delta"] = [int(dl[a > 0].max()) if (a > 0).any() else 0 for a, dl in zip(ac, delta)]
    return d


def load_loeuf(path, column):
    c = pd.read_csv(path, sep="\t", usecols=["gene", "gene_id", "mane_select", column],
                    dtype={"gene": str, "gene_id": str, "mane_select": str}, compression="gzip")
    c = c[(c["mane_select"].str.lower() == "true") & c["gene_id"].str.startswith("ENSG") & c[column].notna()]
    c = c.drop_duplicates("gene")[["gene", column]].rename(columns={column: "LOEUF"})
    q1, q2 = c["LOEUF"].quantile([1 / 3, 2 / 3])
    c["tertile"] = np.where(c["LOEUF"] <= q1, GROUPS[0], np.where(c["LOEUF"] <= q2, GROUPS[1], GROUPS[2]))
    return c, q1, q2


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv-bed", required=True)
    p.add_argument("--genic-tsv", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--constraint", required=True, help="gnomAD constraint metrics table (bgzipped TSV)")
    p.add_argument("--loeuf-column", default="lof.oe_ci.upper")
    p.add_argument("--motif", type=int, required=True)
    p.add_argument("--source", default="TRExplorer")
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-copies", type=float, default=2)
    p.add_argument("--min-loci", type=int, default=20)
    args = p.parse_args()

    d = load_loci(args)
    loeuf, q1, q2 = load_loeuf(args.constraint, args.loeuf_column)
    coding = d[d["genic_context"] == "coding"].assign(gene=lambda x: x["genes"].str.split(",")).explode("gene")
    coding = coding.merge(loeuf, on="gene", how="left")
    n_unmatched = coding.loc[coding["LOEUF"].isna(), "gene"].nunique()
    coding[["ID", "chrom", "start", "end", "copies", "gene", "LOEUF", "tertile", "n_alt", "nonref_AF", "max_abs_delta"]].to_csv(
        f"{args.out_prefix}.genes.tsv", sep="\t", index=False, float_format="%.4g")
    coding = coding[coding["tertile"].notna()].rename(columns={"tertile": "group"})
    inter = d[d["genic_context"] == "intergenic"].assign(group="intergenic")
    allg = pd.concat([coding, inter], ignore_index=True)
    allg["bin"] = pd.cut(allg["copies"], COPY_EDGES, right=False)

    agg = {"n_loci": ("ID", "size"), "copies_median": ("copies", "median")}
    for m, _ in METRICS:
        agg.update({f"{m}_mean": (m, "mean"), f"{m}_median": (m, "median"), f"{m}_sd": (m, "std"),
                    f"{m}_q25": (m, lambda x: x.quantile(0.25)), f"{m}_q75": (m, lambda x: x.quantile(0.75))})
    tab = allg.groupby(["group", "bin"], observed=True).agg(**agg).reset_index()
    tab["bin"] = tab["bin"].astype(str)
    tab.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")

    print(f"LOEUF tertile cutoffs ({args.loeuf_column}, MANE genes n={len(loeuf):,}): T1 <= {q1:.3f} < T2 <= {q2:.3f} < T3")
    print(f"coding loci x gene rows: {len(coding):,} matched; genes without LOEUF: {n_unmatched}")
    print(allg.groupby("group").agg(n_loci=("ID", "size"), n_genes=("gene", "nunique"), nonref_AF=("nonref_AF", "mean"),
                                    n_alt=("n_alt", "mean"), max_abs_delta=("max_abs_delta", "mean"),
                                    copies_median=("copies", "median")).reindex(GROUPS).round(4).to_string())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, (m, ylabel) in zip(axes, METRICS):
        for g in GROUPS:
            t = tab[(tab["group"] == g) & (tab["n_loci"] >= args.min_loci)]
            if m == "max_abs_delta":
                y = t[f"{m}_median"].to_numpy()
                lo, hi = t[f"{m}_q25"].to_numpy(), t[f"{m}_q75"].to_numpy()
            else:
                y = t[f"{m}_mean"].to_numpy()
                ci = 1.96 * t[f"{m}_sd"].fillna(0).to_numpy() / np.sqrt(t["n_loci"].to_numpy())
                lo, hi = y - ci, y + ci
            n_g = int((allg["group"] == g).sum())
            ls = "--" if g == "intergenic" else "-"
            ax.fill_between(t["copies_median"], lo, hi, color=GROUP_COLORS[g], alpha=0.15, linewidth=0)
            ax.plot(t["copies_median"], y, color=GROUP_COLORS[g], linewidth=2, linestyle=ls, marker="o", markersize=3.5,
                    label=f"{g} (n = {n_g:,})")
        ax.set_xscale("log")
        ax.set_xticks([2, 3, 4, 5, 6, 8, 10, 15, 20])
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel(f"Copy number ({args.motif}-bp motif)", fontsize=9, color=INK)
        ax.set_ylabel(ylabel, fontsize=9, color=INK)
        ax.set_ylim(bottom=0)
        style(ax)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left", title="Coding loci by gene LOEUF tertile",
                   title_fontsize=8)
    fig.suptitle(f"{args.label}: coding {args.motif}-bp STRs ({args.source}, PASS, inside CDS) by gene LOEUF tertile "
                 f"(T1 ≤ {q1:.2f} < T2 ≤ {q2:.2f} < T3); bins with ≥{args.min_loci} loci", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
