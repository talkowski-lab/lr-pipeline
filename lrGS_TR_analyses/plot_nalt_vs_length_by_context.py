#!/usr/bin/env python3
"""ALT alleles per TRV site vs repeat length, one figure per genomic context.

Input: merged TRV table (merge_trv_tables.sh output; columns used: FILTER, motifs, ref_len, AC, genic_context,
intergenic_category). Sites: FILTER PASS with summed AC > 0.
Contexts (one PDF each):
  coding, intronic                      genic_context
  intergenic_non_desert                 intergenic and intergenic_category in Tiny / Small / Medium
  intergenic_desert                     intergenic and intergenic_category Desert (>= 460 kb)
Intergenic TRVs outside the listed regions (between_regions, chrom_not_in_region_file) are not plotted; loci shorter than
--min-copies motif copies (REF length / shortest motif; default 2) are dropped.
Each figure: x = repeat (REF) length in bp, y = distinct ALT alleles with AC > 0, both log scale. Grey hexbin density of
all sites, plus the median (line) and IQR (band) per log-length bin for each motif class (1, 2, 3, 4-6, >= 7 bp; bins
with >= --min-n sites). Spearman rho per motif class in the legend.
Outputs: <out-prefix>.<context>.pdf and <out-prefix>.tsv (binned medians).

Usage:
  plot_nalt_vs_length_by_context.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --label hprc_hgsvc --out-prefix out/nalt_vs_length
"""
import argparse

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
MOTIF_CLASSES = [("1 bp", 1, 1), ("2 bp", 2, 2), ("3 bp", 3, 3), ("4–6 bp", 4, 6), ("≥7 bp (VNTR)", 7, 10 ** 6)]
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#4a3aa7"]
NON_DESERT = {"Tiny (<10kb)", "Small (10–100kb)", "Medium (100–460kb)"}
DESERT = {"Desert (≥460kb)"}
CONTEXTS = [("coding", "Coding"), ("intronic", "Intronic"),
            ("intergenic_non_desert", "Intergenic, non-desert (Tiny / Small / Medium)"),
            ("intergenic_desert", "Intergenic, desert (≥460 kb)")]


def load(path):
    cols = ["FILTER", "motifs", "ref_len", "AC", "genic_context", "intergenic_category"]
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=cols, dtype=str, chunksize=500000):
        ch = ch[ch["FILTER"] == "PASS"]
        ac = ch["AC"].str.split(",").map(lambda x: [int(v) for v in x if v not in (".", "")])
        n_alt = ac.map(lambda a: sum(v > 0 for v in a))
        ctx = ch["genic_context"].copy()
        inter = ctx == "intergenic"
        ctx[inter & ch["intergenic_category"].isin(NON_DESERT)] = "intergenic_non_desert"
        ctx[inter & ch["intergenic_category"].isin(DESERT)] = "intergenic_desert"
        parts.append(pd.DataFrame({"context": ctx.to_numpy(), "n_alt": n_alt.to_numpy(),
                                   "ref_len": ch["ref_len"].astype(int).to_numpy(),
                                   "motif": ch["motifs"].map(lambda m: min(len(x) for x in m.split(","))).to_numpy()}))
    d = pd.concat(parts, ignore_index=True)
    return d[(d["n_alt"] > 0) & (d["ref_len"] > 0)]


def keep_copies(d, min_copies):
    return d[d["ref_len"] / d["motif"] >= min_copies]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-n", type=int, default=20)
    p.add_argument("--min-copies", type=float, default=2, help="Drop loci shorter than this many motif copies (default 2)")
    args = p.parse_args()

    d = keep_copies(load(args.trv), args.min_copies)
    edges = np.logspace(0, np.log10(d["ref_len"].max() + 1), 31)
    xmax, ymax = d["ref_len"].max() * 1.2, d["n_alt"].max() * 1.2
    rows = []
    for ctx, title in CONTEXTS:
        s = d[d["context"] == ctx]
        fig, ax = plt.subplots(figsize=(7.5, 5.2))
        ax.hexbin(np.log10(s["ref_len"]), np.log10(s["n_alt"]), gridsize=60, bins="log", cmap="Greys", mincnt=1,
                  linewidths=0, alpha=0.55, extent=(0, np.log10(xmax), 0, np.log10(ymax)))
        for (mlab, lo, hi), color in zip(MOTIF_CLASSES, COLORS):
            m = s[(s["motif"] >= lo) & (s["motif"] <= hi)]
            if len(m) < args.min_n:
                continue
            m = m.assign(bin=pd.cut(m["ref_len"], edges, right=False))
            t = m.groupby("bin", observed=True)["n_alt"].agg(n="size", med="median",
                                                             q25=lambda v: v.quantile(.25), q75=lambda v: v.quantile(.75))
            t["x"] = m.groupby("bin", observed=True)["ref_len"].median()
            t = t[t["n"] >= args.min_n]
            rho = stats.spearmanr(m["ref_len"], m["n_alt"]).statistic if len(m) > 2 else np.nan
            ax.fill_between(np.log10(t["x"]), np.log10(t["q25"]), np.log10(t["q75"]), color=color, alpha=0.15, linewidth=0)
            ax.plot(np.log10(t["x"]), np.log10(t["med"]), color=color, linewidth=2, marker="o", markersize=3,
                    label=f"{mlab} (n = {len(m):,}; ρ = {rho:.2f})")
            for r in t.reset_index().itertuples():
                rows.append({"context": ctx, "motif_class": mlab, "bin": str(r.bin), "length_median": r.x, "n_sites": r.n,
                             "n_alt_median": r.med, "n_alt_q25": r.q25, "n_alt_q75": r.q75})
        ticks = [1, 10, 100, 1000, 10000]
        ax.set_xticks([np.log10(x) for x in ticks if x <= xmax])
        ax.set_xticklabels([f"{x:,}" for x in ticks if x <= xmax])
        yt = [1, 2, 5, 10, 20, 50, 100, 200, 500]
        ax.set_yticks([np.log10(y) for y in yt if y <= ymax])
        ax.set_yticklabels([str(y) for y in yt if y <= ymax])
        ax.set_xlim(0, np.log10(xmax))
        ax.set_ylim(0, np.log10(ymax))
        ax.set_xlabel("Repeat (REF) length, bp (log)", fontsize=9, color=INK)
        ax.set_ylabel("Distinct ALT alleles per site (log)", fontsize=9, color=INK)
        ax.set_title(f"{args.label + ': ' if args.label else ''}{title}\n{len(s):,} non-ref PASS TRVs; "
                     f"grey = site density, lines = median (IQR band) per length bin", fontsize=9.5, color=INK)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper left", title="Motif length", title_fontsize=8)
        fig.tight_layout()
        fig.savefig(f"{args.out_prefix}.{ctx}.pdf")
        plt.close(fig)
        print(f"{ctx}: {len(s):,} sites; median ALT alleles {s['n_alt'].median():g}; "
              f"Spearman rho (all motifs) {stats.spearmanr(s['ref_len'], s['n_alt']).statistic:.3f}")
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    skipped = d[d["context"] == "intergenic"]
    print(f"intergenic sites outside the listed regions (not plotted): {len(skipped):,}")


if __name__ == "__main__":
    main()
