#!/usr/bin/env python3
"""Per-site ALT allele size range (max ALT length - min ALT length) by genomic context for one motif size.

Input: merged TRV table (merge_trv_tables.sh output; columns used: FILTER, motifs, ref_len, alt_len, AC, genic_context,
intergenic_category). Sites: FILTER PASS, >= 1 ALT allele with AC > 0, shortest motif == --motif, REF length / motif >=
--min-copies. Per site: range = max - min of the lengths of ALT alleles with AC > 0 (REF not included; 0 when a single ALT
allele is carried). Contexts: coding, intronic, intergenic non-desert (Tiny / Small / Medium), intergenic desert; intergenic
sites outside the listed regions are not used.
Outputs (prefix --out-prefix):
  .pdf  left: histogram of the range (bp, per-context fraction, ranges above --max-range pooled into the last bin);
        right: cumulative distribution (log x of range + 1); legend: n, median, % range = 0, % range multiple of motif
  .tsv  per context: n_sites, median / mean / p90 range, pct range 0, pct range multiple of motif (among range > 0)

Usage:
  plot_alt_size_range_by_context.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --motif 3 --label hprc_hgsvc --out-prefix out/x
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


def load(path, motif, min_copies):
    cols = ["FILTER", "motifs", "ref_len", "alt_len", "AC", "genic_context", "intergenic_category"]
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=cols, dtype=str, chunksize=500000):
        ch = ch[ch["FILTER"] == "PASS"]
        m = ch["motifs"].map(lambda s: min(len(x) for x in s.split(",")))
        ch = ch[(m == motif) & (ch["ref_len"].astype(int) / motif >= min_copies)]
        if ch.empty:
            continue
        rng = []
        for lens, acs in zip(ch["alt_len"], ch["AC"]):
            carried = [int(L) for L, a in zip(lens.split(","), acs.split(",")) if a not in (".", "") and int(a) > 0]
            rng.append(max(carried) - min(carried) if carried else np.nan)
        ctx = ch["genic_context"].copy()
        inter = ctx == "intergenic"
        ctx[inter] = "."
        ctx[inter & ch["intergenic_category"].isin(NON_DESERT)] = "intergenic_non_desert"
        ctx[inter & ch["intergenic_category"].isin(DESERT)] = "intergenic_desert"
        parts.append(pd.DataFrame({"context": ctx.to_numpy(), "range": rng}))
    d = pd.concat(parts, ignore_index=True).dropna()
    return d[d["context"].isin(CONTEXTS)]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--motif", type=int, default=3)
    p.add_argument("--min-copies", type=float, default=2)
    p.add_argument("--max-range", type=int, default=60, help="Histogram x limit; larger ranges pooled into the last bin")
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    d = load(args.trv, args.motif, args.min_copies)
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    edges = np.arange(-0.5, args.max_range + 1.5, 1)
    for c in CONTEXTS:
        v = d.loc[d["context"] == c, "range"].to_numpy()
        if len(v) == 0:
            continue
        pos = v[v > 0]
        st = {"context": c, "n_sites": len(v), "median_range": float(np.median(v)), "mean_range": float(v.mean()),
              "p90_range": float(np.percentile(v, 90)), "pct_range_0": 100 * float((v == 0).mean()),
              "pct_multiple_of_motif_among_range_gt0": 100 * float((pos % args.motif == 0).mean()) if len(pos) else np.nan}
        rows.append(st)
        lab = (f"{LABELS[c]} (n = {len(v):,}; median {st['median_range']:g} bp; {st['pct_range_0']:.0f}% = 0; "
               f"{st['pct_multiple_of_motif_among_range_gt0']:.0f}% of >0 are ×{args.motif})")
        h = np.histogram(np.clip(v, 0, args.max_range), edges)[0] / len(v)
        axes[0].step(edges[:-1] + 0.5, h, where="mid", color=COLORS[c], linewidth=1.6, label=lab)
        sv = np.sort(v)
        axes[1].step(sv + 1, np.arange(1, len(sv) + 1) / len(sv), where="post", color=COLORS[c], linewidth=1.8,
                     label=LABELS[c])
    axes[0].set_yscale("log")
    axes[0].set_xlabel(f"Max ALT length − min ALT length (bp); last bin = ≥{args.max_range}", fontsize=9, color=INK)
    axes[0].set_ylabel("Fraction of sites (log)", fontsize=9, color=INK)
    axes[0].set_xticks(range(0, args.max_range + 1, args.motif * 2 if args.motif < 5 else args.motif))
    axes[0].legend(frameon=False, fontsize=7, labelcolor=INK, loc="upper right")
    axes[1].set_xscale("log")
    ticks = [1, 2, 4, 7, 10, 31, 101, 1001, 10001]
    axes[1].set_xticks(ticks)
    axes[1].set_xticklabels([str(t - 1) for t in ticks])
    axes[1].set_ylim(0, 1.01)
    axes[1].set_xlabel("Max ALT length − min ALT length (bp)", fontsize=9, color=INK)
    axes[1].set_ylabel("Cumulative fraction of sites", fontsize=9, color=INK)
    axes[1].legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower right")
    for a in axes:
        for sp in ["top", "right"]:
            a.spines[sp].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8)
        a.grid(color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
    fig.suptitle(f"{args.label + ': ' if args.label else ''}{args.motif}-bp motif STRs: ALT allele size range per site "
                 f"(carried ALT alleles only; PASS, non-ref, ≥{args.min_copies:g} copies)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
