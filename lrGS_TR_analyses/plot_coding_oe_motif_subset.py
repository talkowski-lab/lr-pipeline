#!/usr/bin/env python3
"""Obs/exp ALT alleles of coding TRVs by gene LOEUF decile for a chosen set of motif lengths, with a comparison set.

Uses the constraint model of make_trv_deck_figures.py (Poisson GLM of distinct ALT alleles on log size, log size^2,
motif class x catalog, motif class x log size, GC and purity, fitted on all coding TRVs of the cohort), then aggregates
obs/exp per gene and per LOEUF decile (10 equal-size bins over MANE genes) for loci whose shortest motif is in --motifs.
Loci from all catalogs are used (catalog is a model covariate). Multi-gene loci count toward each gene.

Outputs (prefix --out-prefix): .pdf (left: per-gene log2 obs/exp, genes with >= 3 expected; right: pooled obs/exp with
500 gene bootstraps, main set vs --compare-motifs; y capped at 2.2) and .tsv (per-decile values for both sets).

Usage:
  plot_coding_oe_motif_subset.py --loci hprc_hgsvc.TRV_loci.tsv.gz --constraint gnomad.v4.1.1.constraint_metrics.tsv.bgz \
      --motifs 1 2 4 5 7 8 10 --compare-motifs 3 6 9 --label hprc_hgsvc --out-prefix out/x
"""
import argparse

import matplotlib
import numpy as np
import pandas as pd

from make_trv_deck_figures import draw_oe, fit_model, gene_obs_exp, load_loeuf, style

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
MAIN = "#2a78d6"
COMPARE = "#8a8984"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", required=True)
    p.add_argument("--constraint", required=True)
    p.add_argument("--motifs", nargs="+", type=int, required=True)
    p.add_argument("--compare-motifs", nargs="*", type=int, default=[])
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    rng = np.random.default_rng(args.seed)

    d = pd.read_csv(args.loci, sep="\t", dtype={"genes": str})
    coding, model = fit_model(d)
    loeuf = load_loeuf(args.constraint)
    main_set = coding[coding["motif_len"].isin(args.motifs)]
    per_gene, dec = gene_obs_exp(main_set, loeuf, rng)
    name = "motifs " + ",".join(map(str, args.motifs)) + " bp"
    tabs = [dec.assign(motif_set=name)]
    print(f"{args.label}: model on {len(coding):,} coding TRVs; {name}: {len(main_set):,} loci, "
          f"{per_gene['gene'].nunique():,} genes with LOEUF; "
          f"obs {main_set['n_alt'].sum():,} vs exp {main_set['expected'].sum():,.0f}")
    print(main_set["motif_len"].value_counts().sort_index().to_string())

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    draw_oe(axes[0], axes[1], per_gene, dec, MAIN, f"{name} ({len(main_set):,} loci)")
    if args.compare_motifs:
        cmp_set = coding[coding["motif_len"].isin(args.compare_motifs)]
        _, dec_c = gene_obs_exp(cmp_set, loeuf, rng)
        cname = "motifs " + ",".join(map(str, args.compare_motifs)) + " bp"
        tabs.append(dec_c.assign(motif_set=cname))
        axes[1].fill_between(dec_c["decile"], dec_c["lo"], dec_c["hi"], color=COMPARE, alpha=0.15, linewidth=0)
        axes[1].plot(dec_c["decile"], dec_c["obs_exp"], color=COMPARE, linewidth=1.6, linestyle="--", marker="o",
                     markersize=3.5, label=f"{cname} ({len(cmp_set):,} loci)")
    for _, r in dec.iterrows():
        axes[1].text(r["decile"], min(r["hi"], 2.0), f"{int(r['n_loci'])}", fontsize=6.5, color=MAIN, ha="center", va="bottom")
    axes[0].set_title(f"Per-gene obs/exp, {name} ({int((per_gene['exp'] >= 3).sum()):,} genes with ≥3 expected)",
                      fontsize=9.5, color=INK)
    axes[1].set_title("Pooled obs/exp by LOEUF decile (numbers = loci per decile)", fontsize=9.5, color=INK)
    axes[1].set_ylim(0.5, 2.2)
    axes[1].text(10.4, 2.15, "y capped at 2.2; full CIs in TSV", fontsize=7, color=INK, ha="right", va="top")
    axes[1].legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left")
    style(axes[1], grid="both")
    fig.suptitle(f"{args.label}: coding TRVs (inside one CDS block, PASS, single-component, ≥2 copies), "
                 f"obs/exp distinct ALT alleles by gene LOEUF decile (1 = most constrained)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)
    out = pd.concat(tabs)
    out.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")
    cols = ["motif_set", "decile", "LOEUF_median", "n_genes", "n_loci", "obs_exp", "lo", "hi"]
    print(out[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
