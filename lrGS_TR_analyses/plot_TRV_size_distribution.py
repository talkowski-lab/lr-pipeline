#!/usr/bin/env python3
"""Repeat-size distribution of TRV sites in one genic context (default coding), cohorts overlaid.

Input: per-locus tables from build_trv_locus_table.py as LABEL=PATH, colours as LABEL=HEX.
Sites: FILTER PASS, >= 1 ALT allele with AC > 0, genic_context == --context (coding = TR span inside one CDS block).
Size = TRID span (end - start, bp; REF repeat length without VCF padding).
Panels (equal-width bins within each panel, log10 y):
  1-100 bp     sites per exact size (1 bp)
  101-1000 bp  sites per 10-bp bin
  all sizes    cumulative fraction of sites vs size (log x)
Outputs: <out-prefix>.pdf and <out-prefix>.tsv (sites per exact size per cohort).

Usage:
  plot_TRV_size_distribution.py --loci hprc_hgsvc=a.tsv.gz AoU_I=b.tsv.gz --colors hprc_hgsvc=#2a78d6 AoU_I=#eb6834 \
      [--context coding] --out-prefix out/coding_TRV_size_distribution
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


def style(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", nargs="+", required=True, help="LABEL=PATH per-locus tables")
    p.add_argument("--colors", nargs="+", required=True, help="LABEL=HEX")
    p.add_argument("--context", default="coding", choices=["coding", "coding_partial", "UTR", "intronic", "intergenic"])
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()
    colors = dict(c.split("=", 1) for c in args.colors)

    sizes, rows = {}, []
    for item in args.loci:
        label, path = item.split("=", 1)
        d = pd.read_csv(path, sep="\t", usecols=["FILTER", "n_alt", "genic_context", "size"])
        v = d.loc[(d["FILTER"] == "PASS") & (d["n_alt"] > 0) & (d["genic_context"] == args.context), "size"].to_numpy()
        sizes[label] = v
        vc = pd.Series(v).value_counts().sort_index()
        rows += [{"cohort": label, "genic_context": args.context, "size_bp": int(s), "n_sites": int(n)} for s, n in vc.items()]
        q = np.percentile(v, [10, 25, 50, 75, 90, 99])
        print(f"{label}: {len(v):,} sites; size p10/p25/median/p75/p90/p99 = " + "/".join(f"{x:g}" for x in q)
              + f"; max {v.max():,}; >100 bp {int((v > 100).sum()):,}; >1000 bp {int((v > 1000).sum()):,}")
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1.5, 1.2, 1]})
    for label, v in sizes.items():
        c = colors[label]
        # 1-100 bp exact
        xs = np.arange(1, 101)
        y = pd.Series(v).value_counts().reindex(xs).to_numpy(dtype=float)
        y[np.isnan(y) | (y == 0)] = np.nan
        axes[0].plot(xs, y, color=c, linewidth=1.2, marker="o", markersize=2.5,
                     label=f"{label} (n = {int(((v >= 1) & (v <= 100)).sum()):,})")
        # 101-1000 bp, 10-bp bins
        edges = np.arange(101, 1011, 10)
        h = np.histogram(v, edges)[0].astype(float)
        h[h == 0] = np.nan
        axes[1].plot(edges[:-1] + 5, h, color=c, linestyle="none", marker="o", markersize=3.5,
                     label=f"{label} (n = {int(((v > 100) & (v <= 1000)).sum()):,})")
        # ECDF all
        sv = np.sort(v)
        axes[2].step(sv, np.arange(1, len(sv) + 1) / len(sv), where="post", color=c, linewidth=1.8,
                     label=f"{label} (median {np.median(v):g} bp, n = {len(v):,})")
    for ax in axes[:2]:
        ax.set_yscale("log", base=10)
        ax.set_ylabel("TRV sites (log10)", fontsize=9, color=INK)
    axes[0].set_xlabel("Repeat size (bp), exact, 1–100", fontsize=9, color=INK)
    axes[0].set_title("Size 1–100 bp (per bp)", fontsize=10, color=INK)
    axes[1].set_xlabel("Repeat size (bp), 10-bp bins, 101–1000", fontsize=9, color=INK)
    axes[1].set_title("Size 101–1000 bp (per 10 bp)", fontsize=10, color=INK)
    axes[2].set_xscale("log")
    axes[2].set_ylim(0, 1.01)
    axes[2].set_xlabel("Repeat size (bp, log)", fontsize=9, color=INK)
    axes[2].set_ylabel("Cumulative fraction of sites", fontsize=9, color=INK)
    axes[2].set_title("All sizes (cumulative)", fontsize=10, color=INK)
    for ax in axes:
        style(ax)
        ax.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper right" if ax is not axes[2] else "lower right")
    desc = "coding TRV sites (TR span inside one CDS block)" if args.context == "coding" else f"{args.context} TRV sites"
    fig.suptitle(f"Repeat size of non-ref PASS {desc}", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
