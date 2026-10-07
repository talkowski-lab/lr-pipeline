#!/usr/bin/env python3
"""Line plot of motif-length distribution of non-ref TRV sites in one genic context, cohorts overlaid.

Non-ref site = FILTER PASS and summed AC across all ALT alleles > 0 (AC taken from the TRV bed).
Motif length = shortest motif across TRID components (min_motif_len from annotate_TRV_genic_context.sh).
Bins (labelled by upper value v, covering (previous v, v]): 1..10 bp singly, then 20, 30, ..., 100,
then 200, 300, ..., 1000, and so on per decade up to the largest motif observed.

Inputs, per cohort: LABEL=GENIC_CONTEXT_TSV,TRV_BED ; colors as LABEL=HEX.
Bins are plotted at evenly spaced x positions (one slot per bin).
Outputs: <out-prefix>.pdf (left: fraction of sites, linear y; right: same on log y) and <out-prefix>.tsv (counts + fractions).

Usage:
  plot_TRV_motif_length_line.py --inputs hprc_hgsvc=genic.tsv.gz,trv.bed.gz AoU_I=... \
      --colors hprc_hgsvc=#2a78d6 AoU_I=#eb6834 --context coding --out-prefix out/TRV_coding_motif_length
"""
import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"


def bin_uppers(max_len):
    """1..10, 20..100, 200..1000, ... until the bin covering max_len."""
    uppers = list(range(1, 11))
    step = 10
    while uppers[-1] < max_len:
        uppers += list(range(uppers[-1] + step, uppers[-1] * 10 + 1, step))
        step *= 10
    i = next(k for k, v in enumerate(uppers) if v >= max_len)
    return np.array(uppers[: i + 1])


def nonref_ids(trv_bed):
    """IDs of sites with summed AC > 0."""
    keep = []
    for chunk in pd.read_csv(trv_bed, sep="\t", usecols=["ID", "AC"], dtype=str, chunksize=500000):
        s = chunk["AC"].str.split(",").map(lambda x: sum(int(v) for v in x))
        keep.append(chunk.loc[s > 0, "ID"])
    return set(pd.concat(keep))


def motif_lengths(genic_tsv, trv_bed, context):
    ids = nonref_ids(trv_bed)
    g = pd.read_csv(genic_tsv, sep="\t", usecols=["ID", "FILTER", "min_motif_len", "genic_context"],
                    dtype={"ID": str, "FILTER": str, "genic_context": str})
    g = g[(g["FILTER"] == "PASS") & (g["genic_context"] == context) & g["ID"].isin(ids)]
    return g["min_motif_len"].to_numpy()


def style(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", nargs="+", required=True, help="LABEL=GENIC_CONTEXT_TSV,TRV_BED")
    p.add_argument("--colors", nargs="+", required=True, help="LABEL=HEX, one per input label")
    p.add_argument("--context", default="coding", choices=["coding", "UTR", "intronic", "intergenic"])
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    colors = dict(c.split("=", 1) for c in args.colors)
    lengths = {}
    for item in args.inputs:
        label, paths = item.split("=", 1)
        genic_tsv, trv_bed = paths.split(",")
        lengths[label] = motif_lengths(genic_tsv, trv_bed, args.context)
        print(f"{label}: {len(lengths[label]):,} non-ref PASS {args.context} sites; max motif {lengths[label].max()} bp")

    uppers = bin_uppers(max(v.max() for v in lengths.values()))
    edges = np.concatenate([[0], uppers]) + 0.5  # integer bins (prev, v]
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8))
    for label, v in lengths.items():
        n = np.histogram(v, edges)[0]
        frac = n / n.sum()
        for u, k, f in zip(uppers, n, frac):
            rows.append({"cohort": label, "genic_context": args.context, "motif_bin_upper_bp": int(u),
                         "n_sites": int(k), "fraction": round(float(f), 6)})
        x = np.arange(len(uppers))
        for ax in axes:
            ax.plot(x, frac, color=colors[label], linewidth=2, marker="o", markersize=4,
                    label=f"{label} (n = {n.sum():,})")
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)

    for ax, yscale in zip(axes, ["linear", "log"]):
        ax.set_yscale(yscale)
        ax.set_xticks(np.arange(len(uppers)))
        ax.set_xticklabels([str(t) for t in uppers], rotation=60, ha="right", fontsize=7)
        ax.set_xlabel("Motif length (bp; bins 1-10 single bp, then per 10, per 100, ...)", fontsize=9, color=INK)
        ax.set_ylabel(f"Fraction of {args.context} TRV sites" + (" (log)" if yscale == "log" else ""), fontsize=9, color=INK)
        style(ax)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
    fig.suptitle(f"Motif length of non-ref PASS {args.context} TRV sites", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
