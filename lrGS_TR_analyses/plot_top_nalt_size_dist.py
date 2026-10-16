#!/usr/bin/env python3
"""Allele size distribution (ALT size - REF size) at the TRV loci with the most distinct ALT alleles.

Input: merged TRV table (merge_trv_tables.sh output; columns used: #chrom start end ID FILTER motifs ref_len AC AN
alt_size_diff genic_context genes SOURCE intergenic_category). Loci: FILTER PASS, ranked by distinct ALT alleles with
AC > 0 (ties broken by REF length); top --top loci, optionally restricted with --max-motif / --contexts.
Per locus: histogram of the size difference across haplotypes (each ALT allele weighted by its AC; REF haplotypes,
AN - sum(AC), at 0), with the median and range. Panels are titled with ID, gene / intergenic category and motifs.
Outputs: <out-prefix>.pdf (grid) and <out-prefix>.tsv (one row per ALT allele: locus, size difference, AC).

Usage:
  plot_top_nalt_size_dist.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --top 10 --label hprc_hgsvc --out-prefix out/top10
"""
import argparse
import math

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
BAR = "#2a78d6"
REFC = "#eb6834"
COLS = ["#chrom", "start", "end", "ID", "FILTER", "motifs", "ref_len", "AC", "AN", "alt_size_diff", "genic_context", "genes",
        "SOURCE", "intergenic_category"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--top", type=int, default=10)
    p.add_argument("--max-motif", type=int, help="Only loci whose shortest motif is <= this")
    p.add_argument("--contexts", nargs="*", help="Only these genic_context values")
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    parts = []
    for ch in pd.read_csv(args.trv, sep="\t", usecols=COLS, dtype=str, chunksize=500000):
        ch = ch[ch["FILTER"] == "PASS"].copy()
        if args.contexts:
            ch = ch[ch["genic_context"].isin(args.contexts)]
        ch["motif_len"] = ch["motifs"].map(lambda m: min(len(x) for x in m.split(",")))
        if args.max_motif:
            ch = ch[ch["motif_len"] <= args.max_motif]
        ch["n_alt"] = ch["AC"].map(lambda s: sum(int(x) > 0 for x in s.split(",") if x not in (".", "")))
        ch["ref_len_i"] = ch["ref_len"].astype(int)
        parts.append(ch.nlargest(args.top * 3, ["n_alt", "ref_len_i"]))
    top = pd.concat(parts).nlargest(args.top, ["n_alt", "ref_len_i"]).reset_index(drop=True)

    ncol = min(5, args.top)
    nrow = math.ceil(len(top) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.1 * nrow), squeeze=False)
    rows = []
    for k, r in top.iterrows():
        ax = axes[k // ncol][k % ncol]
        ac = np.array([int(x) for x in r["AC"].split(",")])
        diff = np.array([int(x) for x in r["alt_size_diff"].split(",")])
        an = int(r["AN"])
        n_ref = an - ac.sum()
        for dlt, c in zip(diff, ac):
            rows.append({"ID": r["ID"], "size_diff": int(dlt), "AC": int(c)})
        vals = np.concatenate([diff, [0]])
        wts = np.concatenate([ac, [n_ref]]).astype(float)
        lo, hi = int(vals[wts > 0].min()), int(vals[wts > 0].max())
        span = max(hi - lo, 1)
        width = max(1, int(np.ceil(span / 40)))
        bins = np.arange(lo - 0.5, hi + width + 0.5, width)
        ax.hist(diff, bins=bins, weights=ac, color=BAR, alpha=0.85, label="ALT haplotypes")
        if n_ref > 0:
            ax.bar(0, n_ref, width=width, color=REFC, label=f"REF haplotypes ({n_ref})")
        ax.axvline(0, color=INK2, linewidth=0.8, linestyle=":")
        hap = np.repeat(np.concatenate([diff, [0]]), np.concatenate([ac, [max(n_ref, 0)]]).astype(int))
        med = float(np.median(hap)) if len(hap) else np.nan
        where = r["genes"] if r["genic_context"] in ("coding", "UTR", "intronic") else r["intergenic_category"]
        mot = r["motifs"].split(",")
        mot_show = ",".join(m if len(m) <= 8 else m[:8] + "…" for m in mot[:2]) + (f" +{len(mot) - 2}" if len(mot) > 2 else "")
        ax.set_title(f"{r['ID']}\n{r['genic_context']}: {where}; {r['SOURCE']}\nmotif {mot_show}; REF {r['ref_len']} bp",
                     fontsize=7.5, color=INK)
        ax.text(0.98, 0.97, f"{int(r['n_alt'])} ALT alleles / AN {an}\nmedian {med:+g} bp\nrange {lo:+d} to {hi:+d} bp",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color=INK)
        ax.set_xlabel("ALT size − REF size (bp)", fontsize=8, color=INK)
        if k % ncol == 0:
            ax.set_ylabel("Haplotypes", fontsize=8, color=INK)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.tick_params(colors=INK2, labelsize=7)
        ax.grid(axis="y", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
    for k in range(len(top), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    lab = f"{args.label}: " if args.label else ""
    fig.suptitle(f"{lab}allele size distribution at the {len(top)} PASS TRV loci with the most distinct ALT alleles "
                 f"(ALT alleles weighted by AC; REF haplotypes at 0)", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    print(top[["ID", "n_alt", "AN", "ref_len", "motifs", "SOURCE", "genic_context", "genes", "intergenic_category"]]
          .to_string(index=False, max_colwidth=30))


if __name__ == "__main__":
    main()
