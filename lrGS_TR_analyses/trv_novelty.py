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
Plot: TRVs that are PASS, non-ref (n_alt > 0) and have >= 1 routine variant; proportions novel_all and novel_any with
Wilson 95% CI vs TRV
size and vs motif length, one line per genic context (coding = TR span inside one CDS block, intronic, intergenic).
Bins with fewer than --min-n TRVs are not drawn.

Outputs (prefix --out-prefix):
  .bed.gz   one row per TRV (sorted, bgzip-compatible gzip): chrom start end ID FILTER source motif_len size genic_context genes
            n_alt n_variants n_novel frac_novel novel_all novel_any n_snv n_novel_snv
  .pdf      2 x 2 panels: rows = novel_all / novel_any proportion, columns = size / motif; lines = context
  .variant_fraction.pdf  fraction of in-repeat variants that are lrGS-unique: mean of per-TRV frac_novel (+-95% CI, solid)
            and pooled sum(n_novel) / sum(n_variants) (dashed) vs size and motif, lines = context
  .tsv      plotted values per context x bin (n_TRV, novel_all / novel_any proportions with CI, mean / pooled frac_novel)

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
            k2 = int(s["novel_any"].sum())
            lo, hi = wilson(k, n) if n else (np.nan, np.nan)
            lo2, hi2 = wilson(k2, n) if n else (np.nan, np.nan)
            rows.append({"x_variable": col, "genic_context": ctx, "bin": lab, "bin_index": i, "n_TRV": n, "n_novel_all": k,
                         "prop_novel_all": k / n if n else np.nan, "ci_lo": lo, "ci_hi": hi,
                         "prop_novel_any": k2 / n if n else np.nan, "ci_lo_any": lo2, "ci_hi_any": hi2,
                         "mean_frac_novel": float(s["frac_novel"].mean()) if n else np.nan,
                         "se_frac_novel": float(s["frac_novel"].std() / np.sqrt(n)) if n > 1 else np.nan,
                         "n_variants": int(s["n_variants"].sum()), "n_novel_variants": int(s["n_novel"].sum()),
                         "pooled_frac_novel": float(s["n_novel"].sum() / s["n_variants"].sum()) if n else np.nan,
                         "plotted": n >= min_n})
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
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    rows_def = [("prop_novel_all", "ci_lo", "ci_hi", "All in-repeat variants lrGS-unique"),
                ("prop_novel_any", "ci_lo_any", "ci_hi_any", "≥1 in-repeat variant lrGS-unique")]
    for j, (col, edges, xlabel) in enumerate([("size", SIZE_EDGES, "TRV size (bp; TRID span)"),
                                              ("motif_len", MOTIF_EDGES, "Motif length (bp; shortest TRID motif)")]):
        t, labels = proportions(pl, col, edges, args.min_n)
        tabs.append(t)
        for i, (ycol, lo_col, hi_col, ylab) in enumerate(rows_def):
            ax = axes[i][j]
            for ctx in CONTEXTS:
                s = t[(t["genic_context"] == ctx) & t["plotted"]]
                ax.fill_between(s["bin_index"], s[lo_col], s[hi_col], color=COLORS[ctx], alpha=0.18, linewidth=0)
                ax.plot(s["bin_index"], s[ycol], color=COLORS[ctx], linewidth=2, marker="o", markersize=4,
                        label=f"{ctx} (n = {int((pl['genic_context'] == ctx).sum()):,})")
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
            ax.set_xlabel(xlabel, fontsize=9, color=INK)
            ax.set_ylabel(f"Proportion of TRVs\n{ylab}", fontsize=8.5, color=INK)
            ax.set_ylim(0, None)
            ax.set_title(f"{ylab} · by {'TRV size' if col == 'size' else 'motif length'}", fontsize=10, color=INK)
            for sp in ["top", "right"]:
                ax.spines[sp].set_visible(False)
            ax.tick_params(colors=INK2, labelsize=8)
            ax.grid(color=GRID, linewidth=0.6)
            ax.set_axisbelow(True)
            ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left")
    fig.suptitle(f"{args.label}: TRVs novel to gnomAD-LR (in-repeat routine variants lacking a dbGaP and gnomAD v4 match); "
                 f"PASS non-ref TRVs with ≥1 routine variant; 95% Wilson CI; bins with ≥{args.min_n} TRVs",
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)
    pd.concat(tabs).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    plot_variant_fraction(tabs, pl, args)


def plot_variant_fraction(tabs, pl, args):
    """Fraction of in-repeat variants that are lrGS-unique: mean of per-TRV fractions (+-95% CI) and pooled (dashed)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    xdefs = [("size", "TRV size (bp; TRID span)"), ("motif_len", "Motif length (bp; shortest TRID motif)")]
    for ax, t, (col, xlabel) in zip(axes, tabs, xdefs):
        labels = list(dict.fromkeys(t.sort_values("bin_index")["bin"]))
        for ctx in CONTEXTS:
            s = t[(t["genic_context"] == ctx) & t["plotted"]]
            y = s["mean_frac_novel"].to_numpy()
            ci = 1.96 * s["se_frac_novel"].fillna(0).to_numpy()
            ax.fill_between(s["bin_index"], y - ci, y + ci, color=COLORS[ctx], alpha=0.18, linewidth=0)
            ax.plot(s["bin_index"], y, color=COLORS[ctx], linewidth=2, marker="o", markersize=4,
                    label=f"{ctx}: mean per TRV (n = {int((pl['genic_context'] == ctx).sum()):,} TRVs)")
            ax.plot(s["bin_index"], s["pooled_frac_novel"], color=COLORS[ctx], linewidth=1.2, linestyle="--",
                    label=f"{ctx}: pooled over variants")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=8)
        ax.set_xlabel(xlabel, fontsize=9, color=INK)
        ax.set_ylabel("Fraction of in-repeat variants lrGS-unique", fontsize=9, color=INK)
        ax.set_ylim(0, None)
        ax.set_title(f"By {'TRV size' if col == 'size' else 'motif length'}", fontsize=10, color=INK)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper left", ncol=1)
    fig.suptitle(f"{args.label}: fraction of routine variants inside each TRV lacking a dbGaP and gnomAD v4 match "
                 f"(solid = mean of per-TRV fractions ± 95% CI; dashed = pooled; bins with ≥{args.min_n} TRVs)",
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.variant_fraction.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
