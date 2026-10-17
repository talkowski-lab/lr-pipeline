#!/usr/bin/env python3
"""Flag outlier STR alleles per site: alleles outside median +/- K x IQR of that site's allele distribution, per cohort.

Inputs (--input LABEL=KIND=path, one or more cohorts):
  KIND genotypes    extract_trv_genotypes.py .genotypes.tsv (one row per sample and TRV; hapN_SR_units / hapN_len_bp, AL, SD,
                    AP); outliers are reported per sample haplotype with TRGT spanning-read support (SD) and purity (AP)
  KIND allele_dist  split_sr_matched_strs.py .allele_dist.tsv (id gene repunit_count allele_count; aggregate counts only, e.g.
                    for cohorts whose genotypes cannot leave a secure workspace); outliers are reported per allele size with
                    their allele count. The TRV ID is taken from --loci LABEL=path when given.
Optional: --loci LABEL=path (split_sr_matched_strs.py .loci.tsv; adds TRV ID, gene, match_flag), --sr-allele-dist
(X.sr.allele_dist.tsv; adds the gnomAD SR max and the number of outlier alleles above it).
Model per site and cohort: all called alleles (missing haplotypes, e.g. chrX in males, dropped); value = repeat units on the
gnomAD SR scale (genotypes: --measure SR_units, default, or len_bp). IQR is floored at --min-iqr so monomorphic sites do not
flag every 1-unit change. Outlier: value > median + K x IQR (expansion) or < median - K x IQR (contraction).
Outputs (prefix --out-prefix):
  .sites.tsv               per site: for each cohort n alleles, median, IQR used, bounds, n outlier alleles (expansion /
                           contraction), n above SR max, max
  .LABEL.outliers.tsv      per cohort outlier list (genotypes: per sample haplotype; allele_dist: per allele size with count)
  .pdf                     per site strip plots, one column per cohort (colors --colors); alleles within bounds
                           half-transparent, outliers solid; bounds dashed;
                           sites ordered by total outlier alleles

The gnomAD SR distribution can be added as a cohort (e.g. SR=allele_dist=X.sr.allele_dist.tsv, grey) to plot it alongside;
--rank-by then keeps the site order based on the long-read cohorts only.

Usage:
  flag_str_outliers.py --input SR=allele_dist=STR_SRmatch.sr.allele_dist.tsv \\
      hprc_hgsvc=genotypes=STR_SRmatch.hgsvc_hprc.genotypes.tsv \\
      AoU_I=allele_dist=STR_SRmatch.aou_i.allele_dist.tsv \\
      --loci hprc_hgsvc=STR_SRmatch.hgsvc_hprc.loci.tsv AoU_I=STR_SRmatch.aou_i.loci.tsv \\
      --sr-allele-dist STR_SRmatch.sr.allele_dist.tsv --colors '#8a8984' '#2a78d6' '#eb6834' \\
      --rank-by hprc_hgsvc AoU_I --k 10 \\
      --out-prefix STR_SRmatch.outliers_10xIQR
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


def genotype_alleles(path, measure):
    """Per site: DataFrame of haplotype alleles (sample, GT, haplotype, value, len_bp, SD, AP)."""
    gt = pd.read_csv(path, sep="\t", dtype=str)
    out = {}
    for (gid, trv), g in gt.groupby(["gnomAD_id", "ID"], sort=False):
        rows = []
        for r in g.itertuples(index=False):
            al, sd, ap = str(r.AL).split(","), str(r.SD).split(","), str(r.AP).split(",")
            for h in (1, 2):
                v = getattr(r, f"hap{h}_{measure}")
                if pd.isna(v) or v == "":
                    continue
                # AL / SD / AP are per allele in GT order; haploid calls list one value
                i = h - 1 if len(al) > 1 else 0
                rows.append({"sample": r.sample, "GT": r.GT, "haplotype": h, "value": float(v),
                             "len_bp": getattr(r, f"hap{h}_len_bp"), "SD": sd[i] if i < len(sd) else ".",
                             "AP": ap[i] if i < len(ap) else ".", "count": 1})
        out[gid] = (trv, pd.DataFrame(rows))
    return out


def dist_alleles(path, loci):
    """Per site: DataFrame of distinct allele sizes with counts (value, count)."""
    d = pd.read_csv(path, sep="\t")
    out = {}
    for gid, g in d.groupby("id", sort=False):
        trv = loci.at[gid, "TRV"] if loci is not None and gid in loci.index else "."
        out[gid] = (trv, pd.DataFrame({"value": g["repunit_count"].astype(float).to_numpy(),
                                       "count": g["allele_count"].astype(int).to_numpy()}))
    return out


def wpercentile(v, w, qs):
    o = np.argsort(v)
    v, c = np.asarray(v)[o], np.cumsum(np.asarray(w)[o])
    return [v[np.searchsorted(c, q / 100 * c[-1])] for q in qs]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", nargs="+", required=True, help="LABEL=KIND=path; KIND genotypes or allele_dist")
    p.add_argument("--loci", nargs="*", default=[], help="LABEL=path")
    p.add_argument("--sr-allele-dist")
    p.add_argument("--measure", choices=["SR_units", "len_bp"], default="SR_units")
    p.add_argument("--k", type=float, default=10)
    p.add_argument("--min-iqr", type=float, default=1)
    p.add_argument("--rank-by", nargs="*", help="Cohort labels whose outlier counts order the sites (default: all)")
    p.add_argument("--colors", nargs="+", default=["#2a78d6", "#eb6834", "#1baf7a"])
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    loci = {lab: pd.read_csv(path, sep="\t").set_index("id") for lab, path in (x.split("=", 1) for x in args.loci)}
    sr_max_of = {}
    if args.sr_allele_dist:
        sr = pd.read_csv(args.sr_allele_dist, sep="\t")
        sr_max_of = sr[sr["allele_count"] > 0].groupby("id")["repunit_count"].max().to_dict()
    cohorts = []
    for spec in args.input:
        lab, kind, path = spec.split("=", 2)
        data = genotype_alleles(path, args.measure) if kind == "genotypes" else dist_alleles(path, loci.get(lab))
        cohorts.append((lab, kind, data))
    colors = dict(zip([c[0] for c in cohorts], args.colors))

    site_ids = list(dict.fromkeys(g for _, _, data in cohorts for g in data))
    sites, plot = [], {}
    outs = {lab: [] for lab, _, _ in cohorts}
    for gid in site_ids:
        s = {"gnomAD_id": gid}
        for lab, L in loci.items():
            if gid in L.index:
                s.setdefault("gene", L.at[gid, "gene"])
                s[f"{lab}_match_flag"] = L.at[gid, "match_flag"]
        sr_max = sr_max_of.get(gid, np.nan)
        s["SR_max"] = sr_max
        for lab, kind, data in cohorts:
            if gid not in data or data[gid][1].empty:
                continue
            trv, h = data[gid]
            q1, med, q3 = wpercentile(h["value"], h["count"], [25, 50, 75])
            iqr = max(q3 - q1, args.min_iqr)
            lo, hi = med - args.k * iqr, med + args.k * iqr
            o = h[(h["value"] > hi) | (h["value"] < lo)].copy()
            o.insert(0, "ID", trv)
            o.insert(0, "gnomAD_id", gid)
            o["direction"] = np.where(o["value"] > hi, "expansion", "contraction")
            o["fold_over_median"] = (o["value"] / med).round(2) if med else np.nan
            o["above_SR_max"] = (o["value"] > sr_max) if not np.isnan(sr_max) else np.nan
            outs[lab].append(o)
            n_out = int(o["count"].sum())
            s.update({f"{lab}_TRV": trv, f"{lab}_n_alleles": int(h["count"].sum()), f"{lab}_median": med,
                      f"{lab}_IQR_used": iqr, f"{lab}_lower": lo, f"{lab}_upper": hi, f"{lab}_n_outlier": n_out,
                      f"{lab}_n_expansion": int(o.loc[o["direction"] == "expansion", "count"].sum()),
                      f"{lab}_n_contraction": int(o.loc[o["direction"] == "contraction", "count"].sum()),
                      f"{lab}_n_above_SR_max": int(o.loc[o["value"] > sr_max, "count"].sum()) if sr_max_of else np.nan,
                      f"{lab}_max": h["value"].max()})
            if kind == "genotypes":
                s[f"{lab}_n_outlier_samples"] = o["sample"].nunique()
            plot.setdefault(gid, {})[lab] = (trv, h, med, lo, hi, n_out)
        sites.append(s)
    sites = pd.DataFrame(sites)
    ncols = [f"{lab}_n_outlier" for lab, _, _ in cohorts]
    rank = [f"{lab}_n_outlier" for lab in (args.rank_by or [c[0] for c in cohorts])]
    sites["_tot"] = sites[rank].fillna(0).sum(axis=1)
    sites = sites.sort_values(["_tot", "gnomAD_id"], ascending=[False, True]).drop(columns="_tot")
    sites.to_csv(f"{args.out_prefix}.sites.tsv", sep="\t", index=False, float_format="%.4g")
    for lab, kind, _ in cohorts:
        o = pd.concat(outs[lab]).sort_values(["gnomAD_id", "value"], ascending=[True, False])
        if kind == "allele_dist":
            o = o.rename(columns={"count": "allele_count"})
        else:
            o = o.drop(columns="count")
        o.to_csv(f"{args.out_prefix}.{lab}.outliers.tsv", sep="\t", index=False)
        n_sites = int((sites[f"{lab}_n_outlier"] > 0).sum())
        print(f"{lab}: {n_sites} of {int(sites[f'{lab}_n_alleles'].notna().sum())} sites with outliers; "
              f"{int(sites[f'{lab}_n_outlier'].sum())} outlier alleles; {int(sites[f'{lab}_n_above_SR_max'].sum())} above SR max")
    fields = ["match_flag", "median", "n_alleles", "n_outlier", "n_above_SR_max", "max"]
    show = ["gnomAD_id"] + [f"{lab}_{f}" for lab, _, _ in cohorts for f in fields if f"{lab}_{f}" in sites] + ["SR_max"]
    print(sites[sites[ncols].fillna(0).sum(axis=1) > 0][show].to_string(index=False))

    unit = "repeat units (gnomAD SR scale)" if args.measure == "SR_units" else "allele length (bp)"
    rng = np.random.default_rng(0)
    labs = [c[0] for c in cohorts]
    order = list(sites["gnomAD_id"])
    with PdfPages(f"{args.out_prefix}.pdf") as pdf:
        for start in range(0, len(order), 20):
            fig, axes = plt.subplots(4, 5, figsize=(20, 14))
            for ax, gid in zip(axes.flat, order[start:start + 20]):
                title = []
                for j, lab in enumerate(labs):
                    if lab not in plot.get(gid, {}):
                        title.append(f"{lab}: no TRV")
                        continue
                    trv, h, med, lo, hi, n = plot[gid][lab]
                    v, c = h["value"].to_numpy(), h["count"].to_numpy()
                    out = (v > hi) | (v < lo)
                    # allele_dist rows carry counts: point size scales with log count
                    size = 6 + 6 * np.log10(c)
                    x = j + rng.uniform(-0.3, 0.3, len(v))
                    ax.scatter(x[~out], v[~out], s=size[~out], color=colors[lab], alpha=0.35, linewidths=0)
                    ax.scatter(x[out], v[out], s=size[out] + 8, color=colors[lab], alpha=1, edgecolors=INK, linewidths=0.3)
                    ax.hlines(med, j - 0.4, j + 0.4, color=INK2, linewidth=0.9)
                    for b in (lo, hi):
                        if b > 0:
                            ax.hlines(b, j - 0.4, j + 0.4, color=INK2, linewidth=0.8, linestyle="--")
                    title.append(f"{lab} {n:,}")
                if gid in sr_max_of:
                    ax.axhline(sr_max_of[gid], color="#8a8984", linewidth=0.8, linestyle=":")
                ax.set_yscale("log")
                ax.set_xlim(-0.6, len(labs) - 0.4)
                ax.set_xticks(range(len(labs)))
                ax.set_xticklabels(labs, fontsize=7.5, color=INK)
                ax.set_title(f"{gid}\noutlier alleles: " + "; ".join(title), fontsize=8, color=INK)
                ax.tick_params(colors=INK2, labelsize=7)
                ax.grid(axis="y", color=GRID, linewidth=0.6)
                for sp in ["top", "right"]:
                    ax.spines[sp].set_visible(False)
            for ax in axes.flat[len(order[start:start + 20]):]:
                ax.axis("off")
            for ax in axes[:, 0]:
                ax.set_ylabel(unit, fontsize=8, color=INK)
            fig.suptitle(f"Per-site STR outliers: alleles outside median ± {args.k:g} × IQR (IQR floor {args.min_iqr:g}); "
                         "faded = within bounds, solid = outlier; line = median, dashed = bounds, dotted = gnomAD SR max",
                         fontsize=11, color=INK)
            fig.tight_layout(rect=(0, 0, 1, 0.97))
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    main()
