#!/usr/bin/env python3
"""Count of TRV sites vs TR motif length for one cohort, one line per genic context, log10 y axis.

Input: per-locus table from build_trv_locus_table.py.
Sites: FILTER PASS, >= 1 ALT allele with AC > 0, genic_context in --contexts (coding = TR span inside one CDS block).
Motif length = shortest motif across TRID components. Bins (labelled by upper value v, covering (previous v, v]):
1..10 bp singly, then 20, 30, ..., 100, then 200, 300, ... up to the largest motif. Empty bins are not drawn on the log axis.

Outputs: <out-prefix>.pdf and <out-prefix>.tsv (counts per context x bin).

Usage:
  plot_TRV_motif_counts_by_context.py --loci hprc_hgsvc.TRV_loci.tsv.gz --label hprc_hgsvc \
      [--contexts coding intronic intergenic] --out-prefix out/hprc_hgsvc.TRV_counts_by_motif_and_context
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
COLORS = {"coding": "#2a78d6", "UTR": "#eb6834", "intronic": "#1baf7a", "intergenic": "#eda100", "coding_partial": "#4a3aa7"}


def bin_uppers(max_len):
    uppers = list(range(1, 11))
    step = 10
    while uppers[-1] < max_len:
        uppers += list(range(uppers[-1] + step, uppers[-1] * 10 + 1, step))
        step *= 10
    return np.array([u for u in uppers if u < max_len] + [next(u for u in uppers if u >= max_len)])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--contexts", nargs="+", default=["coding", "intronic", "intergenic"])
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    d = pd.read_csv(args.loci, sep="\t", usecols=["FILTER", "n_alt", "genic_context", "motif_len"])
    d = d[(d["FILTER"] == "PASS") & (d["n_alt"] > 0) & d["genic_context"].isin(args.contexts)]
    uppers = bin_uppers(d["motif_len"].max())
    edges = np.concatenate([[0], uppers]) + 0.5
    x = np.arange(len(uppers))

    rows = []
    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.axvspan(-0.4, 0.4, color=GRID, zorder=0)
    for c in args.contexts:
        n = np.histogram(d.loc[d["genic_context"] == c, "motif_len"], edges)[0]
        rows += [{"cohort": args.label, "genic_context": c, "motif_bin_upper_bp": int(u), "n_sites": int(k)}
                 for u, k in zip(uppers, n)]
        y = np.where(n > 0, n, np.nan)
        ax.plot(x, y, color=COLORS[c], linewidth=2, marker="o", markersize=4, label=f"{c} (n = {n.sum():,})")
        print(f"{args.label} {c}: {n.sum():,} sites")
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)

    ax.set_yscale("log", base=10)
    ax.set_xticks(x)
    ax.set_xticklabels([str(u) for u in uppers], rotation=60, ha="right", fontsize=7)
    ax.set_xlabel("TR motif length (bp; bins 1–10 single bp, then per 10, per 100)", fontsize=9, color=INK)
    ax.set_ylabel("Non-ref PASS TRV sites (log10)", fontsize=9, color=INK)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
    ax.set_title(f"{args.label}: non-ref PASS TRV sites by motif length and genic context "
                 f"(coding = inside one CDS block; grey band = homopolymers)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
