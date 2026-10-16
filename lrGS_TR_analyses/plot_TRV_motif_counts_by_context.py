#!/usr/bin/env python3
"""Count of TRV sites vs TR motif length for one cohort, one line per genic context, log10 y axis.

Input: per-locus table from build_trv_locus_table.py.
Sites: FILTER PASS, >= 1 ALT allele with AC > 0, genic_context in --contexts (coding = TR span inside one CDS block).
Motif length = shortest motif across TRID components. Counts are per exact motif length (1-bp resolution) in three
panels: 1-10, 11-100 and 101-1000 bp, sharing one log10 y axis; motif lengths with zero sites are not drawn;
the 101-1000 bp panel shows points only.

Outputs: <out-prefix>.pdf and <out-prefix>.tsv (counts per context x exact motif length).

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


PANELS = [(1, 10), (11, 100), (101, 1000)]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--contexts", nargs="+", default=["coding", "intronic", "intergenic"])
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    d = pd.read_csv(args.loci, sep="\t", usecols=["FILTER", "n_alt", "genic_context", "motif_len"])
    d = d[(d["FILTER"] == "PASS") & (d["n_alt"] > 0) & d["genic_context"].isin(args.contexts)]
    rows = []
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True, gridspec_kw={"width_ratios": [1, 1.6, 1.6]})
    counts = {c: d.loc[d["genic_context"] == c, "motif_len"].value_counts() for c in args.contexts}
    for c, vc in counts.items():
        rows += [{"cohort": args.label, "genic_context": c, "motif_len": int(m), "n_sites": int(n)}
                 for m, n in vc.sort_index().items()]
        print(f"{args.label} {c}: {int(vc.sum()):,} sites")
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    for ax, (lo, hi) in zip(axes, PANELS):
        if lo == 1:
            ax.axvspan(0.6, 1.4, color=GRID, zorder=0)
        for c, vc in counts.items():
            xs = np.arange(lo, hi + 1)
            y = vc.reindex(xs).to_numpy(dtype=float)
            n_panel = int(np.nansum(y))
            y[y == 0] = np.nan
            if lo >= 101:  # sparse: points only
                ax.plot(xs, y, color=COLORS[c], linestyle="none", marker="o", markersize=2.5, alpha=0.7,
                        label=f"{c} (n = {n_panel:,})")
            else:
                ax.plot(xs, y, color=COLORS[c], linewidth=1.6 if lo == 1 else 1, marker="o", markersize=4 if lo == 1 else 2.5,
                        label=f"{c} (n = {n_panel:,})")
        ax.set_yscale("log", base=10)
        ax.set_xlim(lo - 0.6 if lo == 1 else lo - 2, hi + (0.6 if lo == 1 else 2))
        if lo == 1:
            ax.set_xticks(range(1, 11))
        ax.set_xlabel(f"TR motif length (bp), {lo}–{hi}", fontsize=9, color=INK)
        ax.set_title(f"Motif {lo}–{hi} bp", fontsize=10, color=INK)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper right", title="sites in panel", title_fontsize=7.5)
    axes[0].set_ylabel("Non-ref PASS TRV sites per motif length (log10)", fontsize=9, color=INK)
    fig.suptitle(f"{args.label}: non-ref PASS TRV sites by exact motif length and genic context "
                 f"(coding = inside one CDS block; grey band = homopolymers)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
