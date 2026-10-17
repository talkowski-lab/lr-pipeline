#!/usr/bin/env python3
"""Motif-size distribution of STR sites with extreme outlier alleles, stacked by genomic context.

Inputs: sr_envelope_outliers.py .sites.tsv.gz (--sites; ID chrom start end motif_len genic_context n_outlier_alleles
n_motif_expansion n_other_repeat ...), the merged TRV table (--trv; ID -> intergenic_category) and the GTF used for the genic
context (--gtf; CDS records).
Contexts: coding (inside one CDS block), coding-intron junction (overlaps CDS but crosses a CDS boundary), UTR, intronic,
intergenic non-desert (Tiny / Small / Medium), intergenic desert (>= 460 kb); intergenic sites outside the listed regions are
"intergenic other".
--set outlier (default): sites with >= 1 outlier allele; repeat: sites whose outlier alleles include a motif_expansion or
other_repeat allele.
Outputs (prefix --out-prefix):
  .pdf  left: outlier sites per motif length (stacked by context); middle: all STR sites in the same set-up (stacked);
        right: % of STR sites with an outlier per motif length and context
  .tsv  counts per motif length x context (all sites, outlier sites, %)

Usage:
  plot_outlier_motif_by_context.py --sites sr_envelope/hprc_hgsvc.sites.tsv.gz --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz \\
      --gtf r3.gencode.v39.ensembl.105.gtf.gz --label hprc_hgsvc --out-prefix sr_envelope/hprc_hgsvc.outlier_motif_by_context
"""
import argparse
import bisect
import collections
import gzip

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
CONTEXTS = ["coding", "coding_intron_junction", "UTR", "intronic", "intergenic_non_desert", "intergenic_desert",
            "intergenic_other"]
LABELS = {"coding": "Coding", "coding_intron_junction": "Coding-intron junction", "UTR": "UTR", "intronic": "Intronic",
          "intergenic_non_desert": "Intergenic non-desert", "intergenic_desert": "Intergenic desert",
          "intergenic_other": "Intergenic, other"}
COLORS = {"coding": "#2a78d6", "coding_intron_junction": "#8e5bd0", "UTR": "#17a2b8", "intronic": "#1baf7a",
          "intergenic_non_desert": "#eda100", "intergenic_desert": "#8a6d3b", "intergenic_other": "#c9c6bf"}
NON_DESERT = {"Tiny (<10kb)", "Small (10–100kb)", "Medium (100–460kb)"}
DESERT = {"Desert (≥460kb)"}


def cds_blocks(gtf):
    blocks = collections.defaultdict(set)
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            x = line.split("\t", 5)
            if x[2] == "CDS":
                blocks[x[0]].add((int(x[3]) - 1, int(x[4])))
    return {c: (sorted(b), [s for s, _ in sorted(b)]) for c, b in blocks.items()}


def inside_one_block(blocks, chrom, start, end):
    if chrom not in blocks:
        return False
    b, starts = blocks[chrom]
    i = bisect.bisect_right(starts, start)
    # CDS blocks can overlap (several transcripts): check every block starting at or before the site
    return any(s <= start and end <= e for s, e in b[max(0, i - 200):i])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", required=True)
    p.add_argument("--trv", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--set", choices=["outlier", "repeat"], default="outlier")
    p.add_argument("--label", default="")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    s = pd.read_csv(args.sites, sep="\t")
    cat = {}
    for ch in pd.read_csv(args.trv, sep="\t", usecols=["ID", "intergenic_category"], dtype=str, chunksize=500000):
        cat.update(zip(ch["ID"], ch["intergenic_category"]))
    blocks = cds_blocks(args.gtf)
    ctx = s["genic_context"].copy()
    cod = ctx == "coding"
    junction = [not inside_one_block(blocks, c, a, b) for c, a, b in zip(s.loc[cod, "chrom"], s.loc[cod, "start"],
                                                                         s.loc[cod, "end"])]
    ctx[cod] = np.where(junction, "coding_intron_junction", "coding")
    inter = ctx == "intergenic"
    ic = s["ID"].map(cat)
    ctx[inter] = "intergenic_other"
    ctx[inter & ic.isin(NON_DESERT)] = "intergenic_non_desert"
    ctx[inter & ic.isin(DESERT)] = "intergenic_desert"
    s["context"] = ctx
    hit = s["n_outlier_alleles"] > 0
    if args.set == "repeat":
        hit &= (s["n_motif_expansion"] > 0) | (s["n_other_repeat"] > 0)
    s["hit"] = hit

    ks = sorted(s["motif_len"].unique())
    allc = s.groupby(["motif_len", "context"]).size().unstack(fill_value=0).reindex(index=ks, columns=CONTEXTS, fill_value=0)
    hitc = s[s["hit"]].groupby(["motif_len", "context"]).size().unstack(fill_value=0).reindex(
        index=ks, columns=CONTEXTS, fill_value=0)
    pct = (100 * hitc / allc.replace(0, np.nan)).round(3)
    rows = []
    for k in ks:
        for c in CONTEXTS:
            rows.append({"motif_len": k, "context": c, "n_sites": int(allc.at[k, c]), "n_outlier_sites": int(hitc.at[k, c]),
                         "pct": pct.at[k, c]})
    t = pd.DataFrame(rows)
    t.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False)
    tot = t.groupby("context")[["n_sites", "n_outlier_sites"]].sum().reindex(CONTEXTS)
    tot["pct"] = (100 * tot["n_outlier_sites"] / tot["n_sites"]).round(3)
    print(tot.to_string())

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
    for a, data, ttl in [(axes[0], hitc, f"STR sites with an outlier ({int(hitc.values.sum()):,})"),
                         (axes[1], allc, f"All STR sites ({int(allc.values.sum()):,})")]:
        bottom = np.zeros(len(ks))
        for c in CONTEXTS:
            v = data[c].to_numpy()
            if v.sum() == 0:
                continue
            a.bar(ks, v, bottom=bottom, color=COLORS[c], width=0.75, label=f"{LABELS[c]} ({int(v.sum()):,})")
            bottom += v
        a.set_title(ttl, fontsize=9.5, color=INK)
        a.set_ylabel("STR sites", fontsize=9, color=INK)
        a.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper right")
    a = axes[2]
    w = 0.8 / len(CONTEXTS)
    for i, c in enumerate(CONTEXTS):
        a.bar(np.array(ks) + (i - len(CONTEXTS) / 2 + 0.5) * w, pct[c].fillna(0).to_numpy(), width=w, color=COLORS[c],
              label=LABELS[c])
    a.set_title("% of STR sites with an outlier", fontsize=9.5, color=INK)
    a.set_ylabel("% of sites", fontsize=9, color=INK)
    a.legend(frameon=False, fontsize=7.5, labelcolor=INK, loc="upper left")
    for a in axes:
        a.set_xticks(ks)
        a.set_xlabel("Motif length (bp)", fontsize=9, color=INK)
        for sp in ["top", "right"]:
            a.spines[sp].set_visible(False)
        a.tick_params(colors=INK2, labelsize=8)
        a.grid(axis="y", color=GRID, linewidth=0.6)
        a.set_axisbelow(True)
    what = "any outlier allele" if args.set == "outlier" else "a repeat-expansion outlier allele"
    fig.suptitle(f"{args.label + ': ' if args.label else ''}STR sites with {what} beyond the SR-calibrated reach, by motif "
                 "length and genomic context (single-motif STRs, motif 2-6 bp)", fontsize=10.5, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
