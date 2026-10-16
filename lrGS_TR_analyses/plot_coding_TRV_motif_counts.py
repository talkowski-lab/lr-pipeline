#!/usr/bin/env python3
"""Count of TRV sites in one genic context (default coding) vs TR motif length, cohorts overlaid.

Input: per-locus tables from build_trv_locus_table.py, given as LABEL=PATH with a colour LABEL=HEX.
Sites: FILTER PASS, >= 1 ALT allele with AC > 0, genic_context == --context (coding = TR span inside one CDS block;
also UTR, intronic, intergenic, coding_partial).
Motif length = shortest motif across TRID components. Bins (labelled by upper value v, covering (previous v, v]):
1..10 bp singly, then 20, 30, ..., 100, then 200, 300, ... up to the largest motif.

Outputs (prefix --out-prefix): .pdf (+ .png with --png; left: counts, linear y, 1- and 3-bp bins labelled; right: log y;
homopolymer bin shaded) and .tsv with counts per bin.

Usage:
  plot_coding_TRV_motif_counts.py --loci hprc_hgsvc=a.tsv.gz AoU_I=b.tsv.gz --colors hprc_hgsvc=#2a78d6 AoU_I=#eb6834 \
      [--context coding] [--png] --out-prefix out/coding_TRV_motif_counts
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


def bin_uppers(max_len):
    uppers = list(range(1, 11))
    step = 10
    while uppers[-1] < max_len:
        uppers += list(range(uppers[-1] + step, uppers[-1] * 10 + 1, step))
        step *= 10
    return np.array([u for u in uppers if u < max_len] + [next(u for u in uppers if u >= max_len)])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", nargs="+", required=True, help="LABEL=PATH per-locus tables")
    p.add_argument("--colors", nargs="+", required=True, help="LABEL=HEX")
    p.add_argument("--context", default="coding", choices=["coding", "coding_partial", "UTR", "intronic", "intergenic"])
    p.add_argument("--png", action="store_true", help="Also write a PNG")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()
    colors = dict(c.split("=", 1) for c in args.colors)

    motifs = {}
    for item in args.loci:
        label, path = item.split("=", 1)
        d = pd.read_csv(path, sep="\t", usecols=["FILTER", "n_alt", "genic_context", "motif_len"])
        d = d[(d["FILTER"] == "PASS") & (d["n_alt"] > 0) & (d["genic_context"] == args.context)]
        motifs[label] = d["motif_len"].to_numpy()
        print(f"{label}: {len(d):,} non-ref PASS {args.context} TRV sites; max motif {d['motif_len'].max()} bp")

    uppers = bin_uppers(max(v.max() for v in motifs.values()))
    edges = np.concatenate([[0], uppers]) + 0.5
    x = np.arange(len(uppers))
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for label, v in motifs.items():
        n = np.histogram(v, edges)[0]
        rows += [{"cohort": label, "genic_context": args.context, "motif_bin_upper_bp": int(u), "n_sites": int(k),
                  "fraction": round(k / n.sum(), 5)}
                 for u, k in zip(uppers, n)]
        for ax in axes:
            ax.plot(x, n, color=colors[label], linewidth=2, marker="o", markersize=4, label=f"{label} (n = {n.sum():,})")
        for xi in (0, 2):  # 1-bp and 3-bp bins
            axes[0].annotate(f"{n[xi]:,}", (x[xi], n[xi]), textcoords="offset points", xytext=(8, 0), ha="left", va="center",
                             fontsize=8, color=colors[label])
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    for ax, scale in zip(axes, ["linear", "log"]):
        ax.set_yscale(scale)
        ax.axvspan(-0.4, 0.4, color=GRID, zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels([str(u) for u in uppers], rotation=60, ha="right", fontsize=7)
        ax.set_xlabel("TR motif length (bp; bins 1–10 single bp, then per 10, per 100)", fontsize=9, color=INK)
        ax.set_ylabel(f"{args.context} TRV sites" + (" (log)" if scale == "log" else ""), fontsize=9, color=INK)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
    desc = "coding TRV sites (TR span inside one CDS block)" if args.context == "coding" else f"{args.context} TRV sites"
    fig.suptitle(f"Non-ref PASS {desc} by motif length (grey band = homopolymers)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    if args.png:
        fig.savefig(f"{args.out_prefix}.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
