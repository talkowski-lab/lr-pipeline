#!/usr/bin/env python3
"""Association between TRV allele length and CpG methylation at one locus (per sample and per haplotype).

Inputs:
  --vcf       VCF containing the TRV record (per-sample GT; ALT/REF sequences)
  --trv-id    TRV variant ID (e.g. chr6-30949309-TRV-2636)
  --meth-bed  per-haplotype methylation BED for the region (#chrom start end <sample>_hap1 <sample>_hap2 ...; % methylation per
              CpG, "." missing), e.g. streamed from hprc_methylated.<chrom>.haplotype.bed.gz
Methylation per haplotype = mean over CpGs inside the TRV span (POS-1 .. POS-1+len(REF)) with a value; haplotypes with fewer
than --min-cpg CpGs are dropped. Allele length = len(allele) - 1 (VCF padding base); length change = allele length - REF length.
Analyses:
  per sample      mean methylation of the two haplotypes vs mean allele-length change of the two alleles (phase-free)
  per haplotype   phased GTs only: methylation hap1/hap2 paired with GT allele 1/2 (assumes the same haplotype order), with a
                  swapped-pairing control; and, for samples heterozygous for allele length, the paired difference
                  (longer-allele haplotype minus shorter) under both pairings
Outputs (prefix --out-prefix): .pdf, .per_haplotype.tsv, .per_sample.tsv, and printed statistics.

Usage:
  trv_length_vs_methylation.py --vcf locus.vcf.gz --trv-id chr6-30949309-TRV-2636 --meth-bed meth.bed.gz --out-prefix out/x
"""
import argparse
import re
import subprocess

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
BLUE = "#2a78d6"
ORANGE = "#eb6834"


def read_trv(vcf, trv_id):
    q = subprocess.run(["bcftools", "query", "-i", f'ID="{trv_id}"', "-f", "%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n", vcf],
                       check=True, capture_output=True, text=True).stdout.strip().split("\t")
    pos, ref, alts = int(q[0]), q[1], q[2].split(",")
    lengths = [len(a) - 1 for a in [ref] + alts]
    rows = []
    for item in q[3:]:
        sample, gt = item.split("=")
        idx = re.split(r"[|/]", gt)
        phased = "|" in gt
        for h, i in enumerate(idx[:2], start=1):
            rows.append({"sample": sample, "gt_allele_order": h, "allele": None if i == "." else int(i), "phased": phased})
    g = pd.DataFrame(rows)
    g["allele_len"] = g["allele"].map(lambda i: np.nan if pd.isna(i) else lengths[int(i)])
    g["len_change"] = g["allele_len"] - lengths[0]
    return pos, len(ref) - 1, g


def read_meth(path, start, end, min_cpg):
    m = pd.read_csv(path, sep="\t", na_values=".")
    m = m[(m["start"] >= start) & (m["end"] <= end)]
    hap_cols = m.columns[3:]
    n = m[hap_cols].notna().sum()
    mean = m[hap_cols].mean()
    out = pd.DataFrame({"col": hap_cols, "meth": mean.values, "n_cpg": n.values})
    out["sample"] = out["col"].str.replace(r"_hap[12]$", "", regex=True)
    out["meth_hap"] = out["col"].str.extract(r"_hap([12])$")[0].astype(int)
    out.loc[out["n_cpg"] < min_cpg, "meth"] = np.nan
    return out[["sample", "meth_hap", "meth", "n_cpg"]], len(m), m


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True)
    p.add_argument("--trv-id", required=True)
    p.add_argument("--meth-bed", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-cpg", type=int, default=5)
    args = p.parse_args()

    pos, ref_len, g = read_trv(args.vcf, args.trv_id)
    start, end = pos - 1, pos - 1 + ref_len + 1
    meth, n_cpg_region, mraw = read_meth(args.meth_bed, start, end, args.min_cpg)
    print(f"{args.trv_id}: REF repeat {ref_len} bp; CpGs inside span: {n_cpg_region}; "
          f"VCF samples {g['sample'].nunique()}, methylation samples {meth['sample'].nunique()}, "
          f"shared {len(set(g['sample']) & set(meth['sample']))}")
    print("haplotype allele-length change (bp) distribution, VCF samples:")
    print(g["len_change"].value_counts(dropna=False).sort_index().to_string())

    # per haplotype: GT allele h <-> methylation hap h (and swapped control)
    hap = g.merge(meth.rename(columns={"meth_hap": "gt_allele_order"}), on=["sample", "gt_allele_order"], how="inner")
    sw = g.assign(gt_allele_order=3 - g["gt_allele_order"]).merge(
        meth.rename(columns={"meth_hap": "gt_allele_order"}), on=["sample", "gt_allele_order"], how="inner")
    hap.to_csv(f"{args.out_prefix}.per_haplotype.tsv", sep="\t", index=False, float_format="%.4g")
    ph = hap[hap["phased"]].dropna(subset=["meth", "len_change"])
    phs = sw[sw["phased"]].dropna(subset=["meth", "len_change"])
    r1 = stats.spearmanr(ph["len_change"], ph["meth"])
    r2 = stats.spearmanr(phs["len_change"], phs["meth"])
    print(f"\nper haplotype (phased GTs, n = {len(ph)} haplotypes): Spearman rho = {r1.statistic:.3f} (p = {r1.pvalue:.2g}); "
          f"swapped pairing rho = {r2.statistic:.3f} (p = {r2.pvalue:.2g})")
    groups = [s["meth"].to_numpy() for _, s in ph.groupby("len_change") if len(s) >= 5]
    if len(groups) >= 2:
        kw = stats.kruskal(*groups)
        print(f"Kruskal-Wallis across length classes with >= 5 haplotypes ({len(groups)} classes): p = {kw.pvalue:.2g}")
    print(ph.groupby("len_change")["meth"].agg(["size", "median", "mean"]).round(2).to_string())

    # within heterozygous-length phased samples: longer minus shorter haplotype methylation
    diffs = {}
    for name, d in [("as paired", ph), ("swapped", phs)]:
        w = d.pivot_table(index="sample", columns="gt_allele_order", values=["len_change", "meth"])
        w = w.dropna()
        het = w[w[("len_change", 1)] != w[("len_change", 2)]]
        longer_is_1 = het[("len_change", 1)] > het[("len_change", 2)]
        dd = np.where(longer_is_1, het[("meth", 1)] - het[("meth", 2)], het[("meth", 2)] - het[("meth", 1)])
        diffs[name] = dd
        if len(dd) >= 5:
            wt = stats.wilcoxon(dd)
            print(f"length-heterozygous phased samples ({name}): n = {len(dd)}, median(longer - shorter) = "
                  f"{np.median(dd):.2f} pct points, Wilcoxon p = {wt.pvalue:.2g}")

    # per sample, phase-free
    s_len = g.groupby("sample")["len_change"].mean()
    s_meth = meth.groupby("sample")["meth"].mean()
    ps = pd.DataFrame({"mean_len_change": s_len, "mean_meth": s_meth}).dropna()
    ps.to_csv(f"{args.out_prefix}.per_sample.tsv", sep="\t", float_format="%.4g")
    r3 = stats.spearmanr(ps["mean_len_change"], ps["mean_meth"])
    print(f"\nper sample (phase-free, n = {len(ps)}): Spearman rho = {r3.statistic:.3f} (p = {r3.pvalue:.2g}); "
          f"methylation median {ps['mean_meth'].median():.1f}%, IQR {ps['mean_meth'].quantile(.25):.1f}–"
          f"{ps['mean_meth'].quantile(.75):.1f}%")

    # figure
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    ax = axes[0]
    classes = [c for c, s in ph.groupby("len_change") if len(s) >= 3]
    data = [ph.loc[ph["len_change"] == c, "meth"].to_numpy() for c in classes]
    bp = ax.boxplot(data, positions=range(len(classes)), widths=0.6, showfliers=False, patch_artist=True)
    for b in bp["boxes"]:
        b.set(facecolor=BLUE, alpha=0.3, edgecolor=BLUE)
    rng = np.random.default_rng(1)
    for i, v in enumerate(data):
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=8, color=BLUE, alpha=0.6, linewidths=0)
        ax.text(i, ax.get_ylim()[1] if False else np.nanmax(v) + 1.5, f"n={len(v)}", ha="center", fontsize=7, color=INK2)
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels([f"{int(c):+d}" for c in classes], fontsize=8)
    ax.set_xlabel("Allele length change vs REF (bp)", fontsize=9, color=INK)
    ax.set_ylabel("Haplotype methylation in TRV (% , mean over CpGs)", fontsize=9, color=INK)
    ax.set_title(f"Per haplotype (phased GTs, classes with ≥3)\nSpearman ρ = {r1.statistic:.2f} (p = {r1.pvalue:.2g}); "
                 f"swapped ρ = {r2.statistic:.2f}", fontsize=9.5, color=INK)
    ax = axes[1]
    ax.scatter(ps["mean_len_change"], ps["mean_meth"], s=12, color=BLUE, alpha=0.6, linewidths=0)
    ax.set_xlabel("Mean allele length change of the two alleles (bp)", fontsize=9, color=INK)
    ax.set_ylabel("Mean methylation of the two haplotypes (%)", fontsize=9, color=INK)
    ax.set_title(f"Per sample, phase-free (n = {len(ps)})\nSpearman ρ = {r3.statistic:.2f} (p = {r3.pvalue:.2g})",
                 fontsize=9.5, color=INK)
    ax = axes[2]
    hcols = mraw.columns[3:]
    centers = (mraw["start"].to_numpy() + mraw["end"].to_numpy()) / 2
    groups = [([f"{r.sample}_hap{r.gt_allele_order}" for r in ph.itertuples() if r.len_change < 0], ORANGE, "shorter than REF"),
              ([f"{r.sample}_hap{r.gt_allele_order}" for r in ph.itertuples() if r.len_change == 0], BLUE, "REF length"),
              ([f"{r.sample}_hap{r.gt_allele_order}" for r in ph.itertuples() if r.len_change > 0], "#1baf7a", "longer than REF")]
    for cols, color, lab in groups:
        cols = [c for c in cols if c in hcols]
        if cols:
            ax.plot(centers, mraw[cols].mean(axis=1), color=color, marker="o", markersize=3.5, linewidth=1.6,
                    label=f"{lab} ({len(cols)} haplotypes)")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.set_xlabel(f"{mraw['#chrom'].iloc[0]} position (TRV span)", fontsize=9, color=INK)
    ax.set_ylabel("Mean methylation per CpG (%)", fontsize=9, color=INK)
    ax.set_title(f"CpG profile inside the TRV ({n_cpg_region} CpGs)", fontsize=9.5, color=INK)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)
    for a in axes:
        for sp in ["top", "right"]:
            a.spines[sp].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8)
        a.grid(color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
    fig.suptitle(f"{args.trv_id}: allele length vs CpG methylation (hprc_hgsvc)", fontsize=10.5, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
