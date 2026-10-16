#!/usr/bin/env python3
"""STR size (REF allele length) distribution across coding, intronic, intergenic non-desert and intergenic desert.

Input: merged TRV table (merge_trv_tables.sh output; columns used: FILTER, motifs, AC, REF, ref_len, genic_context,
intergenic_category). Sites: FILTER PASS, summed AC > 0, shortest motif in [--min-motif, --max-motif] (default 1-6 bp).
Size: --size ref_allele (default) = len(VCF REF) - 1 (padding base); --size trid_span = catalog repeat span (ref_len column).
Contexts: coding, intronic, intergenic non-desert (Tiny / Small / Medium), intergenic desert (intergenic sites outside the
listed regions are not used).
Outputs (prefix --out-prefix):
  .pdf  left: fraction of sites at each exact size 1-(--max-exact) bp (log y; larger sizes pooled in the last point);
        right: cumulative distribution over all sizes (log x); legend: n and median / p90
  .tsv  per context: n, median, mean, p10 / p25 / p75 / p90 / p99, max

Usage:
  plot_str_ref_length_by_context.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --label hprc_hgsvc --out-prefix out/STR_ref_length
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
CONTEXTS = ["coding", "intronic", "intergenic_non_desert", "intergenic_desert"]
LABELS = {"coding": "Coding", "intronic": "Intronic", "intergenic_non_desert": "Intergenic non-desert",
          "intergenic_desert": "Intergenic desert"}
COLORS = {"coding": "#2a78d6", "intronic": "#1baf7a", "intergenic_non_desert": "#eda100", "intergenic_desert": "#8a6d3b"}
NON_DESERT = {"Tiny (<10kb)", "Small (10–100kb)", "Medium (100–460kb)"}
DESERT = {"Desert (≥460kb)"}


def load(path, min_motif, max_motif, size):
    cols = ["FILTER", "motifs", "AC", "REF", "ref_len", "genic_context", "intergenic_category"]
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=cols, dtype=str, chunksize=500000):
        ch = ch[ch["FILTER"] == "PASS"]
        m = ch["motifs"].map(lambda s: min(len(x) for x in s.split(",")))
        nonref = ch["AC"].map(lambda s: sum(int(x) for x in s.split(",") if x not in (".", "")) > 0)
        ch = ch[nonref & m.between(min_motif, max_motif)]
        ctx = ch["genic_context"].copy()
        inter = ctx == "intergenic"
        ctx[inter] = "."
        ctx[inter & ch["intergenic_category"].isin(NON_DESERT)] = "intergenic_non_desert"
        ctx[inter & ch["intergenic_category"].isin(DESERT)] = "intergenic_desert"
        length = ch["REF"].str.len() - 1 if size == "ref_allele" else ch["ref_len"].astype(int)
        parts.append(pd.DataFrame({"context": ctx.to_numpy(), "length": length.to_numpy()}))
    d = pd.concat(parts, ignore_index=True)
    return d[d["context"].isin(CONTEXTS)]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--min-motif", type=int, default=1)
    p.add_argument("--max-motif", type=int, default=6)
    p.add_argument("--size", choices=["ref_allele", "trid_span"], default="ref_allele")
    p.add_argument("--max-exact", type=int, default=150)
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    d = load(args.trv, args.min_motif, args.max_motif, args.size)
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8), gridspec_kw={"width_ratios": [1.5, 1]})
    xs = np.arange(1, args.max_exact + 1)
    for c in CONTEXTS:
        v = d.loc[d["context"] == c, "length"].to_numpy()
        if len(v) == 0:
            continue
        q = np.percentile(v, [10, 25, 50, 75, 90, 99])
        rows.append({"context": c, "n": len(v), "median": q[2], "mean": v.mean(), "p10": q[0], "p25": q[1], "p75": q[3],
                     "p90": q[4], "p99": q[5], "max": v.max()})
        frac = pd.Series(np.clip(v, 1, args.max_exact)).value_counts().reindex(xs).fillna(0).to_numpy() / len(v)
        y = np.where(frac > 0, frac, np.nan)
        axes[0].plot(xs, y, color=COLORS[c], linewidth=1.2, marker="o", markersize=2,
                     label=f"{LABELS[c]} (n = {len(v):,}; median {q[2]:g} bp; p90 {q[4]:g} bp)")
        sv = np.sort(v)
        axes[1].step(sv, np.arange(1, len(sv) + 1) / len(sv), where="post", color=COLORS[c], linewidth=1.8, label=LABELS[c])
    axes[0].set_yscale("log")
    axes[0].set_xlabel(f"STR size (bp), exact; last point = ≥{args.max_exact}", fontsize=9, color=INK)
    axes[0].set_ylabel("Fraction of STR sites (log)", fontsize=9, color=INK)
    axes[0].legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper right")
    axes[1].set_xscale("log")
    axes[1].set_ylim(0, 1.01)
    axes[1].set_xlabel("STR size (bp, log)", fontsize=9, color=INK)
    axes[1].set_ylabel("Cumulative fraction of STR sites", fontsize=9, color=INK)
    axes[1].legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower right")
    for a in axes:
        for sp in ["top", "right"]:
            a.spines[sp].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8)
        a.grid(color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
    what = "REF allele length" if args.size == "ref_allele" else "catalog repeat span"
    fig.suptitle(f"{args.label + ': ' if args.label else ''}STR size ({what}) by genomic context "
                 f"(motif {args.min_motif}–{args.max_motif} bp; non-ref PASS)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    print(t.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
