#!/usr/bin/env python3
"""STR variability vs copy number for one motif length, split by genic context (coding / UTR / intronic / intergenic).

Inputs: TRV bed (FILTER, SOURCE, TRID, AC, AN) and the matching genic-context table from annotate_TRV_genic_context.sh
(ID -> genic_context). Loci: FILTER PASS, single-component TRID, motif length == --motif, SOURCE == --source,
copy number >= --min-copies.
  copy number  = (TRID end - start) / motif length
  n_alt        = number of distinct ALT alleles with AC > 0
  nonref_AF    = sum(AC) / AN
Loci are binned by copy number (1-copy bins up to --max-single-bin, then wider bins); points with < --min-loci loci
are not plotted.

Outputs (prefix --out-prefix):
  .pdf   panels: mean n_alt and mean nonref_AF (+-95% CI) vs copy number, one line per genic context, plus loci per bin
  .tsv   per context x copy-number bin: n_loci, mean/median n_alt, fraction multi-allelic, mean/median nonref_AF

Usage:
  plot_STR_constraint_by_context.py --trv-bed TRV.bed.gz --genic-tsv genic.tsv.gz --motif 3 --source TRExplorer \
      --label hprc_hgsvc --out-prefix out/hprc_hgsvc.motif3
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
METRICS = [("n_alt", "Mean distinct ALT alleles per locus"), ("nonref_AF", "Mean non-ref allele frequency")]


def load(trv_bed, genic_tsv, motif, source, min_copies):
    parts = []
    for ch in pd.read_csv(trv_bed, sep="\t", usecols=["ID", "FILTER", "SOURCE", "TRID", "AC", "AN"], dtype=str,
                          chunksize=500000):
        ch = ch[(ch["FILTER"] == "PASS") & (ch["SOURCE"] == source) & ~ch["TRID"].str.contains(",")]
        t = ch["TRID"].str.split("-", expand=True)
        ch = ch[t[3].str.len() == motif]
        t = t.loc[ch.index]
        ac = ch["AC"].str.split(",").map(lambda x: [int(v) for v in x])
        parts.append(pd.DataFrame({
            "ID": ch["ID"].to_numpy(),
            "copies": ((t[2].astype(int) - t[1].astype(int)) / motif).to_numpy(),
            "n_alt": ac.map(lambda x: sum(v > 0 for v in x)).to_numpy(),
            "nonref_AF": (ac.map(sum) / ch["AN"].astype(int)).to_numpy()}))
    d = pd.concat(parts, ignore_index=True)
    ctx = pd.read_csv(genic_tsv, sep="\t", usecols=["ID", "genic_context"], dtype=str)
    d = d.merge(ctx, on="ID", how="left")
    assert d["genic_context"].notna().all()
    return d[d["copies"] >= min_copies]


def copy_edges(max_single, max_copies):
    edges = list(range(2, max_single + 1))
    for e in [max_single + 5, max_single + 10, max_single + 20, max_single + 40]:
        if e < max_copies:
            edges.append(e)
    return edges + [np.inf]


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
    p.add_argument("--trv-bed", required=True)
    p.add_argument("--genic-tsv", required=True)
    p.add_argument("--motif", type=int, required=True, help="Motif length (bp)")
    p.add_argument("--source", default="TRExplorer", help="TR catalog (SOURCE column) to use (default TRExplorer)")
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-copies", type=float, default=2)
    p.add_argument("--max-single-bin", type=int, default=12, help="1-copy bins up to this copy number (default 12)")
    p.add_argument("--min-loci", type=int, default=30)
    args = p.parse_args()

    d = load(args.trv_bed, args.genic_tsv, args.motif, args.source, args.min_copies)
    edges = copy_edges(args.max_single_bin, d["copies"].max())
    d["bin"] = pd.cut(d["copies"], edges, right=False)
    tab = d.groupby(["genic_context", "bin"], observed=True).agg(
        n_loci=("n_alt", "size"), copies_median=("copies", "median"),
        n_alt_mean=("n_alt", "mean"), n_alt_median=("n_alt", "median"), n_alt_sd=("n_alt", "std"),
        frac_multiallelic=("n_alt", lambda x: float((x >= 2).mean())),
        nonref_AF_mean=("nonref_AF", "mean"), nonref_AF_median=("nonref_AF", "median"),
        nonref_AF_sd=("nonref_AF", "std")).reset_index()
    tab["bin"] = tab["bin"].astype(str)
    tab.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"{args.label} {args.source} motif {args.motif} bp: {len(d):,} PASS loci")
    print(d.groupby("genic_context").agg(n_loci=("n_alt", "size"), n_alt_mean=("n_alt", "mean"),
                                         nonref_AF_mean=("nonref_AF", "mean"), copies_median=("copies", "median"))
          .reindex(CONTEXTS).round(3).to_string())

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for c in CONTEXTS:
        t = tab[(tab["genic_context"] == c)]
        shown = t[t["n_loci"] >= args.min_loci]
        n_c = int(t["n_loci"].sum())
        for ax, (metric, ylabel) in zip(axes[:2], METRICS):
            y = shown[f"{metric}_mean"].to_numpy()
            ci = 1.96 * shown[f"{metric}_sd"].fillna(0).to_numpy() / np.sqrt(shown["n_loci"].to_numpy())
            ax.fill_between(shown["copies_median"], y - ci, y + ci, color=COLORS[c], alpha=0.18, linewidth=0)
            ax.plot(shown["copies_median"], y, color=COLORS[c], linewidth=2, marker="o", markersize=3.5,
                    label=f"{c} (n = {n_c:,})")
            ax.set_ylabel(ylabel, fontsize=9, color=INK)
        axes[2].plot(t["copies_median"], t["n_loci"], color=COLORS[c], linewidth=2, marker="o", markersize=3.5, label=c)
    for ax in axes:
        ax.set_xscale("log")
        ticks = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]
        ax.set_xticks([x for x in ticks if x <= d["copies"].max()])
        ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel(f"Copy number ({args.motif}-bp motif)", fontsize=9, color=INK)
        style(ax)
    axes[2].set_yscale("log")
    axes[2].set_ylabel("Loci per bin", fontsize=9, color=INK)
    axes[0].legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left")
    fig.suptitle(f"{args.label}: {args.motif}-bp motif STRs ({args.source}, PASS) by genic context; "
                 f"mean ± 95% CI, bins with ≥{args.min_loci} loci", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
