#!/usr/bin/env python3
"""Flag outlier STR alleles per site: haplotype alleles outside median +/- K x IQR of that site's allele distribution.

Input: extract_trv_genotypes.py .genotypes.tsv (one row per sample and TRV; hapN_allele, hapN_len_bp, hapN_SR_units, AL, SD, AP)
and, optionally, the per-dataset loci table from split_sr_matched_strs.py (--loci; adds gene, match_flag) and the SR
allele distribution (--sr-allele-dist, X.sr.allele_dist.tsv; adds SR max and above_SR_max).
Model per TRV: all called haplotypes (missing haplotypes, e.g. chrX in males, dropped); value = --measure (SR_units, default,
or len_bp). IQR is floored at --min-iqr (in the measure's units) so monomorphic sites do not flag every 1-unit change.
Outlier: value > median + K x IQR (expansion) or < median - K x IQR (contraction).
Per outlier the TRGT spanning-read support (SD) and purity (AP) of that allele are reported (AL / SD / AP are listed per allele
in GT order).
Outputs (prefix --out-prefix):
  .sites.tsv     per TRV: gene, match_flag, n haplotypes, median, IQR, bounds, n outlier haplotypes / samples, max, SR max
  .outliers.tsv  per outlier haplotype: sample, GT, haplotype, value, length, fold over median, SD, AP, above_SR_max
  .pdf           per site strip plot of haplotype values (log y), bounds dashed, outliers red; sites with outliers first

Usage:
  flag_str_outliers.py --genotypes STR_SRmatch.hgsvc_hprc.genotypes.tsv --loci STR_SRmatch.hgsvc_hprc.loci.tsv \\
      --sr-allele-dist STR_SRmatch.sr.allele_dist.tsv --k 10 \\
      --out-prefix STR_SRmatch.hgsvc_hprc.outliers_10xIQR
"""
import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
BLUE = "#2a78d6"
RED = "#d03b3b"


def haplotypes(g, measure):
    rows = []
    for r in g.itertuples(index=False):
        al = str(r.AL).split(",")
        sd = str(r.SD).split(",")
        ap = str(r.AP).split(",")
        for h in (1, 2):
            v = getattr(r, f"hap{h}_{measure}")
            if pd.isna(v) or v == "":
                continue
            # AL / SD / AP are per allele in GT order; haploid calls list one value
            i = h - 1 if len(al) > 1 else 0
            rows.append({"sample": r.sample, "GT": r.GT, "haplotype": h, "value": float(v),
                         "len_bp": getattr(r, f"hap{h}_len_bp"), "SD": sd[i] if i < len(sd) else ".",
                         "AP": ap[i] if i < len(ap) else "."})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genotypes", required=True)
    p.add_argument("--loci")
    p.add_argument("--sr-allele-dist")
    p.add_argument("--measure", choices=["SR_units", "len_bp"], default="SR_units")
    p.add_argument("--k", type=float, default=10)
    p.add_argument("--min-iqr", type=float, default=1)
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    gt = pd.read_csv(args.genotypes, sep="\t", dtype=str)
    loci = pd.read_csv(args.loci, sep="\t").set_index("id") if args.loci else None
    sr_max_of = {}
    if args.sr_allele_dist:
        sr = pd.read_csv(args.sr_allele_dist, sep="\t")
        sr_max_of = sr[sr["allele_count"] > 0].groupby("id")["repunit_count"].max().to_dict()
    sites, outs, plot = [], [], []
    for (gid, trv), g in gt.groupby(["gnomAD_id", "ID"], sort=False):
        h = haplotypes(g, args.measure)
        q1, med, q3 = np.percentile(h["value"], [25, 50, 75])
        iqr = max(q3 - q1, args.min_iqr)
        lo, hi = med - args.k * iqr, med + args.k * iqr
        o = h[(h["value"] > hi) | (h["value"] < lo)].copy()
        sr_max = sr_max_of.get(gid, np.nan)
        o.insert(0, "ID", trv)
        o.insert(0, "gnomAD_id", gid)
        o["direction"] = np.where(o["value"] > hi, "expansion", "contraction")
        o["fold_over_median"] = (o["value"] / med).round(2) if med else np.nan
        o["above_SR_max"] = o["value"] > sr_max if not np.isnan(sr_max) else np.nan
        outs.append(o)
        s = {"gnomAD_id": gid, "ID": trv, "n_haplotypes": len(h), "median": med, "Q1": q1, "Q3": q3, "IQR_used": iqr,
             "lower": lo, "upper": hi, "n_outlier_haplotypes": len(o), "n_outlier_samples": o["sample"].nunique(),
             "n_expansion": int((o["direction"] == "expansion").sum()),
             "n_contraction": int((o["direction"] == "contraction").sum()), "max": h["value"].max(), "min": h["value"].min()}
        if loci is not None:
            s.update({"gene": loci.at[gid, "gene"], "match_flag": loci.at[gid, "match_flag"]})
        s.update({"SR_max": sr_max, "n_outlier_above_SR_max": int((o["value"] > sr_max).sum()) if sr_max_of else np.nan})
        sites.append(s)
        plot.append((gid, trv, h, med, lo, hi, len(o)))
    sites = pd.DataFrame(sites).sort_values(["n_outlier_haplotypes", "gnomAD_id"], ascending=[False, True])
    sites.to_csv(f"{args.out_prefix}.sites.tsv", sep="\t", index=False, float_format="%.4g")
    outs = pd.concat(outs).sort_values(["gnomAD_id", "value"], ascending=[True, False])
    outs.to_csv(f"{args.out_prefix}.outliers.tsv", sep="\t", index=False)
    print(sites[sites["n_outlier_haplotypes"] > 0].drop(columns=["ID"]).to_string(index=False))
    print(f"{(sites['n_outlier_haplotypes'] > 0).sum()} of {len(sites)} sites with outliers; {len(outs)} outlier haplotypes "
          f"in {outs['sample'].nunique()} samples")

    order = {g: i for i, g in enumerate(sites["gnomAD_id"])}
    plot.sort(key=lambda x: order[x[0]])
    unit = "repeat units (gnomAD SR scale)" if args.measure == "SR_units" else "allele length (bp)"
    rng = np.random.default_rng(0)
    with PdfPages(f"{args.out_prefix}.pdf") as pdf:
        for start in range(0, len(plot), 20):
            fig, axes = plt.subplots(4, 5, figsize=(17, 13))
            for ax, (gid, trv, h, med, lo, hi, n) in zip(axes.flat, plot[start:start + 20]):
                v = h["value"].to_numpy()
                out = (v > hi) | (v < lo)
                x = rng.uniform(-0.35, 0.35, len(v))
                ax.scatter(x[~out], v[~out], s=6, color=BLUE, alpha=0.5, linewidths=0)
                ax.scatter(x[out], v[out], s=14, color=RED, linewidths=0)
                ax.axhline(med, color=INK2, linewidth=0.8)
                for b in (lo, hi):
                    if b > 0:
                        ax.axhline(b, color=INK2, linewidth=0.8, linestyle="--")
                ax.set_yscale("log")
                ax.set_xlim(-0.6, 0.6)
                ax.set_xticks([])
                ax.set_title(f"{gid}  ({n} outlier hap{'s' if n != 1 else ''})\n{trv}", fontsize=8, color=INK)
                ax.tick_params(colors=INK2, labelsize=7)
                ax.grid(axis="y", color=GRID, linewidth=0.6)
                for sp in ["top", "right", "bottom"]:
                    ax.spines[sp].set_visible(False)
            for ax in axes.flat[len(plot[start:start + 20]):]:
                ax.axis("off")
            for ax in axes[:, 0]:
                ax.set_ylabel(unit, fontsize=8, color=INK)
            fig.suptitle(f"{args.label + ': ' if args.label else ''}per-site STR outliers (haplotype alleles outside median ± "
                         f"{args.k:g} × IQR, IQR floor {args.min_iqr:g}; red = outlier, dashed = bounds)", fontsize=11, color=INK)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    main()
