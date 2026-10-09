#!/usr/bin/env python3
"""Plot non-ref TRV site counts vs motif size, split by genic context (coding / UTR / intronic / intergenic).

Non-ref site = FILTER PASS and summed AC across all ALT alleles > 0 (AC taken from the TRV bed).
Motif size = shortest motif across TRID components (min_motif_len from annotate_TRV_genic_context.sh).

Inputs, per cohort: LABEL=GENIC_CONTEXT_TSV,TRV_BED
Outputs: <out-prefix>.pdf (rows = cohorts, columns = contexts) and <out-prefix>.tsv (counts per bin).

Usage:
  plot_TRV_motif_counts.py --inputs hprc_hgsvc=genic.tsv.gz,trv.bed.gz AoU_I=... --out-prefix out/TRV_motif_counts
"""
import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

CONTEXTS = ["coding", "UTR", "intronic", "intergenic"]
COLORS = {"coding": "#2a78d6", "UTR": "#eb6834", "intronic": "#1baf7a", "intergenic": "#eda100"}
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
STR_BAND = "#f0efec"

# (label, lo, hi) inclusive
MOTIF_BINS = [(str(i), i, i) for i in range(1, 13)]
MOTIF_BINS += [("13-20", 13, 20), ("21-50", 21, 50), ("51-100", 51, 100), (">100", 101, np.inf)]


def nonref_ids(trv_bed):
    """IDs of sites with summed AC > 0."""
    keep = []
    for chunk in pd.read_csv(trv_bed, sep="\t", usecols=["ID", "AC"], dtype=str, chunksize=500000):
        s = chunk["AC"].str.split(",").map(lambda x: sum(int(v) for v in x))
        keep.append(chunk.loc[s > 0, "ID"])
    return set(pd.concat(keep))


def count_table(label, genic_tsv, trv_bed):
    ids = nonref_ids(trv_bed)
    g = pd.read_csv(genic_tsv, sep="\t", usecols=["ID", "FILTER", "min_motif_len", "genic_context"],
                    dtype={"ID": str, "FILTER": str, "genic_context": str})
    n_pass = int((g["FILTER"] == "PASS").sum())
    g = g[(g["FILTER"] == "PASS") & g["ID"].isin(ids)]
    print(f"{label}: PASS {n_pass:,}; non-ref PASS {len(g):,}")
    rows = []
    for c in CONTEXTS:
        ml = g.loc[g["genic_context"] == c, "min_motif_len"].to_numpy()
        for b, lo, hi in MOTIF_BINS:
            n = int(((ml >= lo) & (ml <= hi)).sum())
            rows.append({"cohort": label, "genic_context": c, "motif_bin": b, "n_sites": n,
                         "pct_of_context": round(100 * n / max(len(ml), 1), 3)})
    return pd.DataFrame(rows)


def fmt(v):
    return f"{v / 1e6:.1f}M" if v >= 1e6 else (f"{v / 1e3:.0f}k" if v >= 1e3 else f"{v:.0f}")


def plot(df, out):
    cohorts = list(dict.fromkeys(df["cohort"]))
    fig, axes = plt.subplots(len(cohorts), len(CONTEXTS), figsize=(3.4 * len(CONTEXTS), 2.9 * len(cohorts)),
                             squeeze=False)
    x = np.arange(len(MOTIF_BINS))
    for i, coh in enumerate(cohorts):
        for j, c in enumerate(CONTEXTS):
            ax = axes[i][j]
            sub = df[(df["cohort"] == coh) & (df["genic_context"] == c)]
            v = sub["n_sites"].to_numpy()
            tot = v.sum()
            pct_str = 100 * v[:6].sum() / max(tot, 1)
            ax.axvspan(-0.5, 5.5, color=STR_BAND, zorder=0)
            ax.bar(x, v, width=0.8, color=COLORS[c], edgecolor="white", linewidth=0.5, zorder=2)
            ax.set_title(f"{coh} · {c}\nn = {tot:,}; motif 1–6 bp: {pct_str:.1f}%", fontsize=9, color=INK)
            ax.set_xticks(x)
            ax.set_xticklabels(sub["motif_bin"], rotation=60, ha="right", fontsize=7)
            ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda val, _: fmt(val)))
            for s in ["top", "right"]:
                ax.spines[s].set_visible(False)
            for s in ["left", "bottom"]:
                ax.spines[s].set_color(INK2)
            ax.tick_params(colors=INK2, labelsize=7)
            ax.grid(axis="y", color=GRID, linewidth=0.6, zorder=1)
            ax.set_axisbelow(True)
            if j == 0:
                ax.set_ylabel("Non-ref TRV sites", fontsize=9, color=INK)
            if i == len(cohorts) - 1:
                ax.set_xlabel("Motif size (bp)", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", nargs="+", required=True, help="LABEL=GENIC_CONTEXT_TSV,TRV_BED")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    tables = []
    for item in args.inputs:
        label, paths = item.split("=", 1)
        genic_tsv, trv_bed = paths.split(",")
        tables.append(count_table(label, genic_tsv, trv_bed))
    df = pd.concat(tables)
    df.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    plot(df, f"{args.out_prefix}.pdf")


if __name__ == "__main__":
    main()
