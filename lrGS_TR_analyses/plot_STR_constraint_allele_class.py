#!/usr/bin/env python3
"""STR variability vs copy number for one motif length, split by ALT-allele class and genic context.

Inputs: TRV bed (ID, FILTER, SOURCE, TRID, AC, AN), the genic-context table from annotate_TRV_genic_context.sh
(ID, genic_context, alt_len_diffs in ALT order) and a GENCODE GTF (CDS blocks).
Loci: FILTER PASS, single-component TRID, motif length == --motif, SOURCE == --source, copy number >= --min-copies.
Coding loci are kept only if the TR span lies entirely inside one merged CDS block (so frame is defined); coding loci
that straddle a CDS boundary are dropped (count printed).

Per ALT allele, delta = len(ALT) - len(REF):
  inframe        delta != 0 and delta % 3 == 0
  frameshift     delta % 3 != 0
  no_len_change  delta == 0
Per locus and class: n_alt = distinct ALT alleles of that class with AC > 0; AF = sum(AC of that class) / AN.

Outputs (prefix --out-prefix):
  .pdf  rows = mean n_alt / mean AF (+-95% CI); columns = allele class (independent y axes); one line per genic context
  .tsv  per context x copy-number bin x class: n_loci, mean n_alt, mean AF, fraction of loci with >= 1 allele of class

Usage:
  plot_STR_constraint_allele_class.py --trv-bed TRV.bed.gz --genic-tsv genic.tsv.gz --gtf gencode.gtf.gz --motif 3 \
      --label hprc_hgsvc --out-prefix out/hprc_hgsvc.motif3.allele_class
"""
import argparse
import gzip

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

CONTEXTS = ["coding", "UTR", "intronic", "intergenic"]
CLASSES = ["inframe", "frameshift", "no_len_change"]
CLASS_TITLES = {"inframe": "In-frame alleles (Δ ≠ 0, Δ % 3 = 0)", "frameshift": "Frameshift alleles (Δ % 3 ≠ 0)",
                "no_len_change": "No length change (Δ = 0)"}
COLORS = {"coding": "#2a78d6", "UTR": "#eb6834", "intronic": "#1baf7a", "intergenic": "#eda100"}
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"


def cds_blocks(gtf):
    rows = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.split("\t", 5)
            if t[2] == "CDS":
                rows.append((t[0], int(t[3]) - 1, int(t[4])))
    blocks = {}
    for chrom, sub in pd.DataFrame(rows, columns=["chrom", "s", "e"]).sort_values(["chrom", "s"]).groupby("chrom"):
        merged = []
        for s, e in zip(sub["s"], sub["e"]):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        m = np.array(merged)
        blocks[chrom] = (m[:, 0], m[:, 1])
    return blocks


def cds_contained(chroms, starts, ends, blocks):
    out = np.zeros(len(starts), dtype=bool)
    for chrom in np.unique(chroms):
        if chrom not in blocks:
            continue
        idx = np.where(chroms == chrom)[0]
        bs, be = blocks[chrom]
        k = np.searchsorted(bs, starts[idx], side="right") - 1
        ok = k >= 0
        out[idx[ok]] = ends[idx[ok]] <= be[k[ok]]
    return out


def load(args):
    parts = []
    for ch in pd.read_csv(args.trv_bed, sep="\t", usecols=["#CHROM", "ID", "FILTER", "SOURCE", "TRID", "AC", "AN"], dtype=str,
                          chunksize=500000):
        ch = ch[(ch["FILTER"] == "PASS") & (ch["SOURCE"] == args.source) & ~ch["TRID"].str.contains(",")]
        t = ch["TRID"].str.split("-", expand=True)
        keep = t[3].str.len() == args.motif
        ch, t = ch[keep], t[keep]
        parts.append(pd.DataFrame({"chrom": ch["#CHROM"].to_numpy(), "ID": ch["ID"].to_numpy(), "AC": ch["AC"].to_numpy(),
                                   "AN": ch["AN"].astype(int).to_numpy(), "start": t[1].astype(int).to_numpy(),
                                   "end": t[2].astype(int).to_numpy()}))
    d = pd.concat(parts, ignore_index=True)
    d["copies"] = (d["end"] - d["start"]) / args.motif
    d = d[d["copies"] >= args.min_copies]
    ctx = pd.read_csv(args.genic_tsv, sep="\t", usecols=["ID", "genic_context", "alt_len_diffs"], dtype=str)
    d = d.merge(ctx, on="ID", how="left")
    assert d["genic_context"].notna().all()

    contained = cds_contained(d["chrom"].to_numpy(), d["start"].to_numpy(), d["end"].to_numpy(), cds_blocks(args.gtf))
    drop = (d["genic_context"] == "coding").to_numpy() & ~contained
    print(f"coding loci straddling a CDS boundary dropped: {int(drop.sum()):,}")
    d = d[~drop].reset_index(drop=True)

    ac = d["AC"].str.split(",").map(lambda x: np.array(x, dtype=int))
    delta = d["alt_len_diffs"].str.split(",").map(lambda x: np.array(x, dtype=int))
    for c in CLASSES:
        n_c, s_c = [], []
        for a, dl in zip(ac, delta):
            m = (dl != 0) & (dl % 3 == 0) if c == "inframe" else ((dl % 3 != 0) if c == "frameshift" else (dl == 0))
            m &= a > 0
            n_c.append(int(m.sum()))
            s_c.append(int(a[m].sum()))
        d[f"n_alt_{c}"] = n_c
        d[f"AF_{c}"] = np.array(s_c) / d["AN"].to_numpy()
    return d


def style(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv-bed", required=True)
    p.add_argument("--genic-tsv", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--motif", type=int, required=True)
    p.add_argument("--source", default="TRExplorer")
    p.add_argument("--label", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--min-copies", type=float, default=2)
    p.add_argument("--min-loci", type=int, default=30)
    args = p.parse_args()

    d = load(args)
    edges = list(range(2, 13)) + [17, 22, np.inf]
    d["bin"] = pd.cut(d["copies"], edges, right=False)
    agg = {"n_loci": ("ID", "size"), "copies_median": ("copies", "median")}
    for c in CLASSES:
        agg.update({f"n_alt_{c}_mean": (f"n_alt_{c}", "mean"), f"n_alt_{c}_sd": (f"n_alt_{c}", "std"),
                    f"AF_{c}_mean": (f"AF_{c}", "mean"), f"AF_{c}_sd": (f"AF_{c}", "std"),
                    f"frac_loci_with_{c}": (f"n_alt_{c}", lambda x: float((x > 0).mean()))})
    tab = d.groupby(["genic_context", "bin"], observed=True).agg(**agg).reset_index()
    tab["bin"] = tab["bin"].astype(str)
    tab.to_csv(f"{args.out_prefix}.tsv", sep="\t", index=False, float_format="%.4g")

    overall = d.groupby("genic_context")[[f"AF_{c}" for c in CLASSES] + [f"n_alt_{c}" for c in CLASSES]].mean()
    print(f"{args.label} {args.source} motif {args.motif} bp: {len(d):,} PASS loci")
    print(pd.concat([d.groupby("genic_context").size().rename("n_loci"), overall], axis=1).reindex(CONTEXTS).round(4).to_string())

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for j, c in enumerate(CLASSES):
        for i, (metric, ylabel) in enumerate([(f"n_alt_{c}", "Mean distinct ALT alleles per locus"),
                                              (f"AF_{c}", "Mean summed AF of class")]):
            ax = axes[i][j]
            for ctx in CONTEXTS:
                t = tab[(tab["genic_context"] == ctx) & (tab["n_loci"] >= args.min_loci)]
                y = t[f"{metric}_mean"].to_numpy()
                ci = 1.96 * t[f"{metric}_sd"].fillna(0).to_numpy() / np.sqrt(t["n_loci"].to_numpy())
                n_ctx = int((d["genic_context"] == ctx).sum())
                ax.fill_between(t["copies_median"], y - ci, y + ci, color=COLORS[ctx], alpha=0.18, linewidth=0)
                ax.plot(t["copies_median"], y, color=COLORS[ctx], linewidth=2, marker="o", markersize=3.5,
                        label=f"{ctx} (n = {n_ctx:,})")
            ax.set_xscale("log")
            ax.set_xticks([2, 3, 4, 5, 6, 8, 10, 15, 20])
            ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            if j == 0:
                ax.set_ylabel(ylabel, fontsize=9, color=INK)
            if i == 0:
                ax.set_title(CLASS_TITLES[c], fontsize=10, color=INK)
            else:
                ax.set_xlabel(f"Copy number ({args.motif}-bp motif)", fontsize=9, color=INK)
            style(ax)
    for a in axes.flat:
        a.set_ylim(0, a.get_ylim()[1])
    axes[0][0].legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left")
    fig.suptitle(f"{args.label}: {args.motif}-bp motif STRs ({args.source}, PASS) by allele class and genic context "
                 f"(coding = locus inside one CDS block); mean ± 95% CI, bins with ≥{args.min_loci} loci", fontsize=10, color=INK)
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
