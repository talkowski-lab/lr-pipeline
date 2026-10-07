#!/usr/bin/env python3
"""lrGS-unique variation inside TRVs: per-TRV BED and proportion of novel TRVs by size, motif and genic context.

Inputs:
  --loci     per-locus table from build_trv_locus_table.py (ID, chrom, start, end, FILTER, source, motif_len, size, n_alt,
             genic_context, genes)
  --counts   per-TRV counts from count_variants_in_trv.sh (TRV_ID, n_variants, n_novel, n_snv, n_novel_snv): PASS routine
             (non-TRV) variants with AC > 0 whose TRID column is the TRV, and how many have "." in both dbGaP_ID and
             gnomAD_V4_match_ID (lrGS-unique)
Per TRV:
  frac_novel  n_novel / n_variants
  novel_all   1 if n_variants > 0 and every routine variant inside it is lrGS-unique (TR locus unseen by short-read resources)
  novel_any   1 if >= 1 routine variant inside it is lrGS-unique
Plot: TRVs that are PASS, non-ref (n_alt > 0) and have >= 1 routine variant; proportion novel_all with Wilson 95% CI vs TRV
size and vs motif length, one line per genic context (coding = TR span inside one CDS block, intronic, intergenic).
Bins with fewer than --min-n TRVs are not drawn.

Outputs (prefix --out-prefix):
  .bed.gz   one row per TRV (sorted, bgzip-compatible gzip): chrom start end ID FILTER source motif_len size genic_context genes
            n_alt n_variants n_novel frac_novel novel_all novel_any n_snv n_novel_snv
  .pdf      2 panels (size, motif) x lines = context
  .tsv      plotted proportions per context x bin (n_TRV, n_novel_all, prop_novel_all, CI, prop_novel_any, mean frac_novel)

Usage:
  trv_novelty.py --loci hprc_hgsvc.TRV_loci.tsv.gz --counts hprc_hgsvc.variants_in_TRV.tsv.gz --label hprc_hgsvc \
      --out-prefix trv_novelty/hprc_hgsvc.TRV_novelty
"""
import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
CONTEXTS = ["coding", "intronic", "intergenic"]
COLORS = {"coding": "#2a78d6", "intronic": "#1baf7a", "intergenic": "#eda100"}
SIZE_EDGES = [1, 10, 15, 20, 30, 50, 100, 250, 500, 1000, np.inf]
MOTIF_EDGES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 51, 101, np.inf]


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return c - h, c + h


def bin_label(lo, hi):
    if not np.isfinite(hi):
        return f"≥{lo:g}"
    return f"{lo:g}" if hi - lo == 1 else f"{lo:g}–{hi - 1:g}"


def proportions(d, col, edges, min_n):
    rows = []
    labels = [bin_label(lo, hi) for lo, hi in zip(edges[:-1], edges[1:])]
    d = d.assign(bin=pd.cut(d[col], edges, right=False, labels=labels))
    for ctx in CONTEXTS:
        for i, lab in enumerate(labels):
            s = d[(d["genic_context"] == ctx) & (d["bin"] == lab)]
            n = len(s)
            k = int(s["novel_all"].sum())
            lo, hi = wilson(k, n) if n else (np.nan, np.nan)
            rows.append({"x_variable": col, "genic_context": ctx, "bin": lab, "bin_index": i, "n_TRV": n, "n_novel_all": k,
                         "prop_novel_all": k / n if n else np.nan, "ci_lo": lo, "ci_hi": hi,
                         "prop_novel_any": float(s["novel_any"].mean()) if n else np.nan,
                         "mean_frac_novel": float(s["frac_novel"].mean()) if n else np.nan, "plotted": n >= min_n})
    return pd.DataFrame(rows), labels


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", required=True)
    p.add_argument("--counts", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-n", type=int, default=20)
    args = p.parse_args()

    loci = pd.read_csv(args.loci, sep="\t", dtype={"genes": str, "chrom": str},
                       usecols=["ID", "chrom", "start", "end", "FILTER", "source", "motif_len", "size", "n_alt",
                                "genic_context", "genes"])
    cnt = pd.read_csv(args.counts, sep="\t").rename(columns={"TRV_ID": "ID"})
    unmatched = int((~cnt["ID"].isin(loci["ID"])).sum())
    d = loci.merge(cnt, on="ID", how="left")
    for c in ["n_variants", "n_novel", "n_snv", "n_novel_snv"]:
        d[c] = d[c].fillna(0).astype(int)
    d["frac_novel"] = np.where(d["n_variants"] > 0, d["n_novel"] / d["n_variants"].clip(lower=1), np.nan)
    d["novel_all"] = ((d["n_variants"] > 0) & (d["n_novel"] == d["n_variants"])).astype(int)
    d["novel_any"] = (d["n_novel"] > 0).astype(int)

    chrom_order = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23, "chrY": 24}
    bed = d.assign(o=d["chrom"].map(chrom_order).fillna(99)).sort_values(["o", "start"]).drop(columns="o")
    bed = bed[["chrom", "start", "end", "ID", "FILTER", "source", "motif_len", "size", "genic_context", "genes", "n_alt",
               "n_variants", "n_novel", "frac_novel", "novel_all", "novel_any", "n_snv", "n_novel_snv"]]
    bed.rename(columns={"chrom": "#chrom"}).to_csv(f"{args.out_prefix}.bed.gz", sep="\t", index=False, float_format="%.4g",
                                                   compression="gzip")

    pl = d[(d["FILTER"] == "PASS") & (d["n_alt"] > 0) & (d["n_variants"] > 0) & d["genic_context"].isin(CONTEXTS)]
    print(f"{args.label}: {len(d):,} TRVs; counts rows not matching a TRV: {unmatched}; "
          f"PASS non-ref TRVs with >=1 routine variant: {len(pl):,} of {int(((d.FILTER == 'PASS') & (d.n_alt > 0)).sum()):,}")
    print(pl.groupby("genic_context").agg(n_TRV=("ID", "size"), prop_novel_all=("novel_all", "mean"),
                                          prop_novel_any=("novel_any", "mean"), mean_frac_novel=("frac_novel", "mean"),
                                          variants_per_TRV=("n_variants", "mean")).reindex(CONTEXTS).round(4).to_string())

    tabs = []
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    for ax, (col, edges, xlabel) in zip(axes, [("size", SIZE_EDGES, "TRV size (bp; TRID span)"),
                                               ("motif_len", MOTIF_EDGES, "Motif length (bp; shortest TRID motif)")]):
        t, labels = proportions(pl, col, edges, args.min_n)
        tabs.append(t)
        for ctx in CONTEXTS:
            s = t[(t["genic_context"] == ctx) & t["plotted"]]
            ax.fill_between(s["bin_index"], s["ci_lo"], s["ci_hi"], color=COLORS[ctx], alpha=0.18, linewidth=0)
            ax.plot(s["bin_index"], s["prop_novel_all"], color=COLORS[ctx], linewidth=2, marker="o", markersize=4,
                    label=f"{ctx} (n = {int((pl['genic_context'] == ctx).sum()):,})")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
        ax.set_xlabel(xlabel, fontsize=9, color=INK)
        ax.set_ylabel("Proportion of TRVs novel to gnomAD-LR", fontsize=9, color=INK)
        ax.set_ylim(0, None)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left")
    axes[0].set_title("By TRV size", fontsize=10, color=INK)
    axes[1].set_title("By motif length", fontsize=10, color=INK)
    fig.suptitle(f"{args.label}: TRVs whose in-repeat variants all lack a dbGaP / gnomAD v4 match "
                 f"(PASS non-ref TRVs with ≥1 routine variant; 95% Wilson CI; bins with ≥{args.min_n} TRVs)",
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)
    pd.concat(tabs).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
