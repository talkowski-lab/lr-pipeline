#!/usr/bin/env python3
"""Plot TRV site counts and size distributions by genic context (coding / UTR / intronic / intergenic).

Input: one or more per-cohort tables from annotate_TRV_genic_context.sh, given as LABEL=PATH.
Outputs (all prefixed with --out-prefix):
  .counts.pdf               site counts per genic context, one panel per cohort
  .locus_length.pdf         ECDF of reference locus span (END-START)
  .allele_length_change.pdf per-ALT-allele length change (len(ALT)-len(REF)), each distinct allele once,
                            in fixed integer bins, plus fraction of non-zero changes that are a multiple of 3 bp
  .motif_length.pdf         shortest motif length across TRID components
  .summary.tsv              per cohort x context: counts and size summary stats

Usage:
  plot_TRV_genic_context.py --inputs hprc_hgsvc=a.tsv.gz AoU_I=b.tsv.gz --out-prefix out/TRV_genic_context [--all-filters]
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

MAX_DIFF = 10000
# Integer allele-length-change bins: (label, lo, hi) inclusive
DIFF_BINS = [("<=-100", -np.inf, -100), ("-99..-10", -99, -10), ("-9..-2", -9, -2), ("-1", -1, -1), ("0", 0, 0),
             ("+1", 1, 1), ("+2..+9", 2, 9), ("+10..+99", 10, 99), (">=+100", 100, np.inf)]
MOTIF_BINS = list(range(1, 11)) + [">10"]


def style(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def init_acc():
    return {c: {"n": 0, "locus_vals": [],
                "diff": np.zeros(len(DIFF_BINS)), "n_allele": 0, "n_zero": 0, "n_nonzero": 0, "n_mult3": 0,
                "n_clip": 0, "diff_abs_vals": [], "motif": dict.fromkeys(MOTIF_BINS, 0)} for c in CONTEXTS}


def accumulate(path, pass_only):
    acc = init_acc()
    cols = ["FILTER", "min_motif_len", "locus_len", "alt_len_diffs", "genic_context"]
    for chunk in pd.read_csv(path, sep="\t", usecols=cols, chunksize=500000,
                             dtype={"FILTER": str, "alt_len_diffs": str, "genic_context": str}):
        if pass_only:
            chunk = chunk[chunk["FILTER"] == "PASS"]
        for c, sub in chunk.groupby("genic_context"):
            a = acc[c]
            a["n"] += len(sub)
            locus = sub["locus_len"].to_numpy()
            a["locus_vals"].append(locus)
            d = np.array(",".join(sub["alt_len_diffs"]).split(","), dtype=np.int64)
            a["n_allele"] += len(d)
            a["n_zero"] += int((d == 0).sum())
            nz = d[d != 0]
            a["n_nonzero"] += len(nz)
            a["n_mult3"] += int((nz % 3 == 0).sum())
            a["n_clip"] += int((np.abs(d) > MAX_DIFF).sum())
            a["diff"] += [int(((d >= lo) & (d <= hi)).sum()) for _, lo, hi in DIFF_BINS]
            a["diff_abs_vals"].append(np.abs(nz))
            ml = sub["min_motif_len"].to_numpy()
            for b in MOTIF_BINS[:-1]:
                a["motif"][b] += int((ml == b).sum())
            a["motif"][">10"] += int((ml > 10).sum())
    for c in CONTEXTS:
        a = acc[c]
        a["locus_vals"] = np.concatenate(a["locus_vals"]) if a["locus_vals"] else np.array([])
        a["diff_abs_vals"] = np.concatenate(a["diff_abs_vals"]) if a["diff_abs_vals"] else np.array([])
    return acc


def legend(ax, **kw):
    leg = ax.legend(frameon=False, fontsize=8, labelcolor=INK, **kw)
    return leg


def plot_counts(data, out):
    fig, axes = plt.subplots(1, len(data), figsize=(3.6 * len(data), 3.4), squeeze=False)
    for ax, (label, acc) in zip(axes[0], data.items()):
        n = [acc[c]["n"] for c in CONTEXTS]
        tot = sum(n)
        bars = ax.bar(CONTEXTS, n, color=[COLORS[c] for c in CONTEXTS], width=0.7, edgecolor="white", linewidth=1)
        for b, v in zip(bars, n):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v:,}\n({100 * v / tot:.1f}%)",
                    ha="center", va="bottom", fontsize=7, color=INK)
        ax.set_ylim(0, max(n) * 1.22)
        ax.set_title(f"{label}  (n = {tot:,})", fontsize=10, color=INK)
        ax.set_ylabel("TRV sites", fontsize=9, color=INK)
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x / 1e6:.1f}M"))
        style(ax)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_locus_ecdf(data, out):
    fig, axes = plt.subplots(1, len(data), figsize=(4.2 * len(data), 3.4), squeeze=False, sharey=True)
    for ax, (label, acc) in zip(axes[0], data.items()):
        for c in CONTEXTS:
            v = np.sort(acc[c]["locus_vals"])
            if len(v) == 0:
                continue
            ax.step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=COLORS[c], linewidth=2,
                    label=f"{c} (median {np.median(v):.0f} bp)")
        ax.set_xscale("log")
        ax.set_ylim(0, 1.01)
        ax.set_xlabel("Reference locus length (bp)", fontsize=9, color=INK)
        ax.set_title(label, fontsize=10, color=INK)
        style(ax)
        ax.grid(axis="x", color=GRID, linewidth=0.6)
        legend(ax, loc="lower right")
    axes[0][0].set_ylabel("Cumulative fraction of TRV sites", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_allele_change(data, out):
    fig, axes = plt.subplots(2, len(data), figsize=(5.2 * len(data), 6.2), squeeze=False,
                             gridspec_kw={"height_ratios": [2.2, 1]})
    x = np.arange(len(DIFF_BINS))
    w = 0.2
    for j, (label, acc) in enumerate(data.items()):
        ax = axes[0][j]
        for i, c in enumerate(CONTEXTS):
            h = acc[c]["diff"]
            ax.bar(x + (i - 1.5) * w, h / h.sum(), width=w, color=COLORS[c], label=c, edgecolor="white", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([b[0] for b in DIFF_BINS], rotation=35, ha="right")
        ax.set_xlabel("Allele length change vs REF (bp)", fontsize=9, color=INK)
        ax.set_title(label, fontsize=10, color=INK)
        style(ax)
        legend(ax, loc="upper left")
        ax2 = axes[1][j]
        frac3 = [100 * acc[c]["n_mult3"] / max(acc[c]["n_nonzero"], 1) for c in CONTEXTS]
        bars = ax2.bar(CONTEXTS, frac3, color=[COLORS[c] for c in CONTEXTS], width=0.7, edgecolor="white", linewidth=1)
        for b, v in zip(bars, frac3):
            ax2.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}%", ha="center", va="bottom", fontsize=7, color=INK)
        ax2.axhline(100 / 3, color=INK2, linewidth=0.8, linestyle="--")
        ax2.set_ylim(0, 105)
        ax2.set_ylabel("% non-zero changes\nthat are multiple of 3 bp", fontsize=8, color=INK)
        style(ax2)
    axes[0][0].set_ylabel("Fraction of ALT alleles", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_motif(data, out):
    fig, axes = plt.subplots(1, len(data), figsize=(4.8 * len(data), 3.4), squeeze=False, sharey=True)
    x = np.arange(len(MOTIF_BINS))
    w = 0.2
    for ax, (label, acc) in zip(axes[0], data.items()):
        for i, c in enumerate(CONTEXTS):
            v = np.array([acc[c]["motif"][b] for b in MOTIF_BINS], dtype=float)
            ax.bar(x + (i - 1.5) * w, v / v.sum(), width=w, color=COLORS[c], label=c, edgecolor="white", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([str(b) for b in MOTIF_BINS])
        ax.set_xlabel("Motif length (bp, shortest TRID motif)", fontsize=9, color=INK)
        ax.set_title(label, fontsize=10, color=INK)
        style(ax)
        legend(ax, loc="upper right")
    axes[0][0].set_ylabel("Fraction of TRV sites", fontsize=9, color=INK)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def summary(data, out):
    rows = []
    for label, acc in data.items():
        tot = sum(acc[c]["n"] for c in CONTEXTS)
        for c in CONTEXTS:
            a = acc[c]
            lv, dv = a["locus_vals"], a["diff_abs_vals"]
            m = a["motif"]
            rows.append({
                "cohort": label, "genic_context": c, "n_sites": a["n"], "pct_sites": round(100 * a["n"] / tot, 3),
                "locus_len_median": float(np.median(lv)) if len(lv) else np.nan,
                "locus_len_p90": float(np.percentile(lv, 90)) if len(lv) else np.nan,
                "n_alt_alleles": a["n_allele"], "mean_alt_alleles_per_site": round(a["n_allele"] / max(a["n"], 1), 2),
                "pct_alleles_len_change_0": round(100 * a["n_zero"] / max(a["n_allele"], 1), 2),
                "abs_len_change_median_nonzero": float(np.median(dv)) if len(dv) else np.nan,
                "abs_len_change_p90_nonzero": float(np.percentile(dv, 90)) if len(dv) else np.nan,
                "pct_nonzero_change_mult3": round(100 * a["n_mult3"] / max(a["n_nonzero"], 1), 2),
                "n_alleles_abs_change_gt_10kb": a["n_clip"],
                "pct_motif_1_6bp": round(100 * sum(m[b] for b in range(1, 7)) / max(a["n"], 1), 2),
                "pct_motif_3bp": round(100 * m[3] / max(a["n"], 1), 2),
            })
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    return df


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", nargs="+", required=True, help="LABEL=PATH genic-context tables")
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--all-filters", action="store_true", help="Include non-PASS sites (default: PASS only)")
    args = p.parse_args()

    data = {}
    for item in args.inputs:
        label, path = item.split("=", 1)
        data[label] = accumulate(path, pass_only=not args.all_filters)

    pre = args.out_prefix
    plot_counts(data, f"{pre}.counts.pdf")
    plot_locus_ecdf(data, f"{pre}.locus_length.pdf")
    plot_allele_change(data, f"{pre}.allele_length_change.pdf")
    plot_motif(data, f"{pre}.motif_length.pdf")
    print(summary(data, f"{pre}.summary.tsv").to_string(index=False))


if __name__ == "__main__":
    main()
