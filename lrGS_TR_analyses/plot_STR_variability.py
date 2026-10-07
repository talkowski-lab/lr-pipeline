#!/usr/bin/env python3
"""STR locus variability vs STR size and motif size.

Input: TRV bed (gnomAD_LR.<cohort>.vep_parsed.annotated.TRV.bed.gz; columns include FILTER, TRID, AC, AN).
Loci: FILTER PASS, single-component TRID, motif length (from TRID) <= --max-motif, copy number >= --min-copies.
Everything is split by SOURCE (TR catalog: TRExplorer vs Vamos). The catalogs cover different size ranges and differ
several-fold in variability at the same size, so pooled curves show artifactual steps where the catalog mix changes.
  STR size     = TRID end - start (bp, repeat span without VCF padding);  copy number = size / motif length.
  n_alt        = number of distinct ALT alleles with AC > 0.
  nonref_AF    = sum(AC) / AN  (fraction of genotyped haplotypes carrying any ALT allele).
Loci are binned by size (bp) and by copy number within each motif length; bins with < --min-loci loci are not plotted.

Outputs (prefix --out-prefix):
  .lines.pdf    rows = n_alt / nonref_AF (mean +-95% CI); columns = source x (size bp, copy number); one line per motif
  .by_motif.pdf x = motif length; rows = n_alt / nonref_AF; columns = source x (lines = copy-number bins, lines = size bins);
                tests motif effect at fixed copy number vs at fixed bp size
  .heatmap.pdf  rows = source; motif x size-bin heatmaps of mean n_alt and mean nonref_AF (loci counts in the TSV)
  .tsv          per source x motif x bin: n_loci, mean/median n_alt, fraction with >= 2 ALT alleles, mean/median nonref_AF

Usage:
  plot_STR_variability.py --trv-bed gnomAD_LR.hprc_hgsvc.vep_parsed.annotated.TRV.bed.gz --label hprc_hgsvc --out-prefix out/x
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
MOTIF_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SIZE_EDGES = [1, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 250, 500, np.inf]
COPY_EDGES = [0, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, np.inf]
BYMOTIF_COPY_EDGES = [2, 3, 4, 5, 7, 10, 15, 25, np.inf]
BYMOTIF_SIZE_EDGES = [10, 15, 20, 30, 50, 100, np.inf]
METRICS = [("n_alt", "Mean distinct ALT alleles per locus"), ("nonref_AF", "Mean non-ref allele frequency")]


def load(trv_bed, max_motif, min_copies):
    parts = []
    for ch in pd.read_csv(trv_bed, sep="\t", usecols=["FILTER", "SOURCE", "TRID", "AC", "AN"], dtype=str, chunksize=500000):
        ch = ch[(ch["FILTER"] == "PASS") & ~ch["TRID"].str.contains(",")]
        t = ch["TRID"].str.split("-", expand=True)
        motif = t[3].str.len()
        size = t[2].astype(int) - t[1].astype(int)
        ac = ch["AC"].str.split(",").map(lambda x: [int(v) for v in x])
        parts.append(pd.DataFrame({
            "source": ch["SOURCE"].to_numpy(), "motif": motif.to_numpy(), "size": size.to_numpy(),
            "n_alt": ac.map(lambda x: sum(v > 0 for v in x)).to_numpy(),
            "nonref_AF": (ac.map(sum) / ch["AN"].astype(int)).to_numpy()}))
    d = pd.concat(parts, ignore_index=True)
    d = d[d["motif"] <= max_motif]
    d["copies"] = d["size"] / d["motif"]
    return d[d["copies"] >= min_copies]


def binned(d, col, edges):
    d = d.assign(bin=pd.cut(d[col], edges, right=False))
    g = d.groupby(["source", "motif", "bin"], observed=True)
    out = g.agg(n_loci=("n_alt", "size"), x_median=(col, "median"),
                n_alt_mean=("n_alt", "mean"), n_alt_median=("n_alt", "median"), n_alt_sd=("n_alt", "std"),
                frac_multiallelic=("n_alt", lambda x: float((x >= 2).mean())),
                nonref_AF_mean=("nonref_AF", "mean"), nonref_AF_median=("nonref_AF", "median"),
                nonref_AF_sd=("nonref_AF", "std")).reset_index()
    out["bin"] = out["bin"].astype(str)
    out["x_variable"] = col
    return out


def style(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def plot_lines(tabs, label, min_loci, out):
    sources = sorted(tabs["size"]["source"].unique())
    cols = [(col, xlabel, src) for col, xlabel in [("size", "STR size (bp)"), ("copies", "Copy number (size / motif length)")]
            for src in sources]
    fig, axes = plt.subplots(2, len(cols), figsize=(4.3 * len(cols), 7.6), squeeze=False)
    for j, (col, xlabel, src) in enumerate(cols):
        t = tabs[col]
        t = t[(t["n_loci"] >= min_loci) & (t["source"] == src)]
        for i, (metric, ylabel) in enumerate(METRICS):
            ax = axes[i][j]
            for m in sorted(t["motif"].unique()):
                s = t[t["motif"] == m]
                y = s[f"{metric}_mean"].to_numpy()
                ci = 1.96 * s[f"{metric}_sd"].fillna(0).to_numpy() / np.sqrt(s["n_loci"].to_numpy())
                color = MOTIF_COLORS[m - 1]
                ax.fill_between(s["x_median"], y - ci, y + ci, color=color, alpha=0.18, linewidth=0)
                ax.plot(s["x_median"], y, color=color, linewidth=2, marker="o", markersize=3.5, label=f"{m} bp")
            ax.set_xscale("log")
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=9, color=INK)
            if i == 0:
                ax.set_title(src, fontsize=10, color=INK)
            else:
                ax.set_xlabel(xlabel, fontsize=9, color=INK)
            style(ax)
    for i in range(2):
        lo = min(a.get_ylim()[0] for a in axes[i])
        hi = max(a.get_ylim()[1] for a in axes[i])
        for a in axes[i]:
            a.set_ylim(lo, hi)
    axes[0][0].legend(title="Motif", frameon=False, fontsize=8, title_fontsize=8, labelcolor=INK, loc="upper left")
    fig.suptitle(f"{label}: STR variability vs size and motif by TR catalog "
                 f"(PASS loci; mean ± 95% CI; bins with ≥{min_loci} loci)",
                 fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def bin_label(b):
    lo, hi = b.left, b.right
    if np.isinf(hi):
        return f"≥{lo:g}"
    return f"{lo:g}" if hi - lo == 1 else f"{lo:g}–{hi - 1:g}"


def plot_by_motif(d, label, min_loci, out):
    sources = sorted(d["source"].unique())
    strata = [("copies", BYMOTIF_COPY_EDGES, "copies"), ("size", BYMOTIF_SIZE_EDGES, "bp")]
    cols = [(col, edges, unit, src) for col, edges, unit in strata for src in sources]
    fig, axes = plt.subplots(2, len(cols), figsize=(4.3 * len(cols), 7.6), squeeze=False)
    motifs = sorted(d["motif"].unique())
    for j, (col, edges, unit, src) in enumerate(cols):
        sub = d[d["source"] == src].assign(bin=pd.cut(d.loc[d["source"] == src, col], edges, right=False))
        g = sub.groupby(["bin", "motif"], observed=True).agg(
            n=("n_alt", "size"), n_alt_mean=("n_alt", "mean"), n_alt_sd=("n_alt", "std"),
            nonref_AF_mean=("nonref_AF", "mean"), nonref_AF_sd=("nonref_AF", "std")).reset_index()
        g = g[g["n"] >= min_loci]
        bins = sorted(g["bin"].unique(), key=lambda b: b.left)
        cmap = plt.get_cmap("Blues")
        for i, (metric, ylabel) in enumerate(METRICS):
            ax = axes[i][j]
            for k, b in enumerate(bins):
                t = g[g["bin"] == b].sort_values("motif")
                if len(t) < 2:
                    continue
                y = t[f"{metric}_mean"].to_numpy()
                ci = 1.96 * t[f"{metric}_sd"].fillna(0).to_numpy() / np.sqrt(t["n"].to_numpy())
                color = cmap(0.35 + 0.65 * k / max(len(bins) - 1, 1))
                ax.fill_between(t["motif"], y - ci, y + ci, color=color, alpha=0.2, linewidth=0)
                ax.plot(t["motif"], y, color=color, linewidth=2, marker="o", markersize=3.5, label=f"{bin_label(b)} {unit}")
            ax.set_xticks(motifs)
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=9, color=INK)
            if i == 0:
                ax.set_title(f"{src}: fixed {'copy number' if col == 'copies' else 'STR size'}", fontsize=10, color=INK)
                ax.legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper left",
                          title="Copy number" if col == "copies" else "STR size", title_fontsize=7)
            else:
                ax.set_xlabel("Motif length (bp)", fontsize=9, color=INK)
            style(ax)
    for i in range(2):
        hi = max(a.get_ylim()[1] for a in axes[i])
        for a in axes[i]:
            a.set_ylim(0, hi)
    fig.suptitle(f"{label}: STR variability vs motif length at fixed copy number or fixed size "
                 f"(PASS loci; mean ± 95% CI; cells with ≥{min_loci} loci)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_heatmap(tab, label, min_loci, out):
    bins = sorted(tab["bin"].unique(), key=lambda b: float(b.strip("[)").split(",")[0]))
    motifs = sorted(tab["motif"].unique())
    sources = sorted(tab["source"].unique())
    fig, axes = plt.subplots(len(sources), 2, figsize=(13, 3.6 * len(sources)), squeeze=False)
    for r_i, src in enumerate(sources):
        t = tab[tab["source"] == src]
        for ax, (metric, title) in zip(axes[r_i], METRICS):
            mat = np.full((len(motifs), len(bins)), np.nan)
            for r in t.itertuples():
                if r.n_loci >= min_loci:
                    mat[motifs.index(r.motif), bins.index(r.bin)] = getattr(r, f"{metric}_mean")
            vmax = 0.9 if metric == "nonref_AF" else 30
            im = ax.imshow(mat, aspect="auto", cmap="Blues", origin="lower", vmin=0, vmax=vmax)
            for a in range(mat.shape[0]):
                for b in range(mat.shape[1]):
                    if not np.isnan(mat[a, b]):
                        v = mat[a, b]
                        txt = f"{v:.1f}" if metric == "n_alt" else f"{v:.2f}"
                        ax.text(b, a, txt, ha="center", va="center", fontsize=6.5, color="white" if v > vmax * 0.6 else INK)
            ax.set_xticks(range(len(bins)))
            ax.set_xticklabels([b.replace("inf", "∞") for b in bins], rotation=45, ha="right", fontsize=7)
            ax.set_yticks(range(len(motifs)))
            ax.set_yticklabels([f"{m} bp" for m in motifs], fontsize=8)
            ax.set_xlabel("STR size bin (bp)", fontsize=9, color=INK)
            ax.set_ylabel("Motif length", fontsize=9, color=INK)
            ax.set_title(f"{src}: {title}", fontsize=10, color=INK)
            fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle(f"{label}: PASS STR loci (blank = fewer than {min_loci} loci)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv-bed", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--max-motif", type=int, default=6)
    p.add_argument("--min-loci", type=int, default=100)
    p.add_argument("--min-copies", type=float, default=2, help="Drop loci shorter than this many motif copies (default 2)")
    args = p.parse_args()

    d = load(args.trv_bed, args.max_motif, args.min_copies)
    print(f"{args.label}: {len(d):,} PASS single-component STR loci (motif <= {args.max_motif} bp, >= {args.min_copies} copies)")
    tabs = {"size": binned(d, "size", SIZE_EDGES), "copies": binned(d, "copies", COPY_EDGES)}
    pd.concat(tabs.values()).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    overall = d.groupby(["source", "motif"]).agg(
        n_loci=("n_alt", "size"), n_alt_mean=("n_alt", "mean"),
        frac_multiallelic=("n_alt", lambda x: (x >= 2).mean()),
        nonref_AF_mean=("nonref_AF", "mean"),
        spearman_size_n_alt=("n_alt", lambda x: x.corr(d.loc[x.index, "size"], method="spearman")))
    print(overall.round(3).to_string())
    plot_lines(tabs, args.label, args.min_loci, f"{args.out_prefix}.lines.pdf")
    plot_heatmap(tabs["size"], args.label, args.min_loci, f"{args.out_prefix}.heatmap.pdf")
    plot_by_motif(d, args.label, args.min_loci, f"{args.out_prefix}.by_motif.pdf")


if __name__ == "__main__":
    main()
