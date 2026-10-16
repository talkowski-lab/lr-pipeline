#!/usr/bin/env python3
"""Locus plot of STR / TR loci across a gene or region: catalog tracks on top, observed TRVs, gene structure at the bottom.

Panels (shared x axis, genomic coordinates):
  1. One track per TR catalog (--catalog LABEL=path, TRGT-format BED: chrom start end ID=...;MOTIFS=...). Catalog loci in
     the window are drawn as rectangles: filled with the catalog colour when the locus has a non-ref PASS TRV in --trv whose
     SOURCE equals LABEL (matched on the TRID span), light grey otherwise. Right of each track: total ALT alleles (distinct
     ALT alleles with AC > 0, summed over that catalog's observed TRVs in the window) and observed / catalog loci.
  2. Observed TRVs: lollipops at the TR span midpoint, height = distinct ALT alleles (AC > 0), coloured by catalog.
  3. Gene structure from --gtf: exons as rectangles (CDS tall, UTR / non-coding exon short), introns as lines with strand
     arrows, gene names; overlapping genes are stacked.
--trv is the merged TRV table (merge_trv_tables.sh output; tabix-indexed; columns #chrom start end ID FILTER TRID motifs AC
SOURCE). Window: --gene NAME (transcript span from the GTF +- --flank) or --region chr:start-end; or batch mode --gene-table T
(trv_gene_table.py output) with --count-column (default n_coding) >= --min-count (default 10): one PDF per gene in
--out-dir plus all genes in --out (multi-page PDF).

Usage:
  plot_trv_locus.py --trv gnomAD_LR.hprc_hgsvc.TRV.bed.gz --gtf gencode.gtf.gz \
      --catalog TRExplorer=variation_clusters_and_isolated_TRs_v1.0.1.hg38.TRGT.bed.gz Vamos=vamos_catalog.ori.v2.1.TRGT.bed.gz \
      --gene MUCL3 --flank 5000 --label hprc_hgsvc --out MUCL3.trv_locus.pdf
  plot_trv_locus.py ... --gene-table STR_gene_table.tsv --min-count 10 --max-motif 6 --out-dir plots/ --out all_genes.pdf
"""
import argparse
import gzip
import os
import re
import subprocess

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42

INK = "#1d1f24"
INK2 = "#5b5e66"
GRID = "#e6e4df"
UNOBSERVED = "#d9d7d2"
CAT_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
GENE_COLOR = "#3a3f4b"


def gene_span(gtf, name):
    """Span of the gene's transcript record(s) (the MANE / canonical transcript in the gnomAD-SV GTF), 0-based half-open."""
    spans = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            t = line.split("\t")
            if len(t) > 8 and t[2] == "transcript" and f'gene_name "{name}";' in t[8]:
                spans.append((t[0], int(t[3]) - 1, int(t[4])))
    if not spans:
        raise SystemExit(f"gene {name} not found in {gtf}")
    return spans[0][0], min(s[1] for s in spans), max(s[2] for s in spans)


def read_gtf_window(gtf, chrom, start, end):
    feats = []
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            t = line.rstrip("\n").split("\t")
            if t[0] != chrom or t[2] not in ("transcript", "exon", "CDS"):
                continue
            s, e = int(t[3]) - 1, int(t[4])
            if e <= start or s >= end:
                continue
            gene = re.search(r'gene_name "([^"]+)"', t[8]).group(1)
            tx = re.search(r'transcript_id "([^"]+)"', t[8]).group(1)
            feats.append((gene, tx, t[2], s, e, t[6]))
    return pd.DataFrame(feats, columns=["gene", "tx", "feature", "start", "end", "strand"])


def read_catalog(path, chrom, start, end):
    alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    awk = f'($1=="{chrom}" || $1=="{alt}") && $2 < {end} && $3 > {start}'
    out = subprocess.run(f"gzip -cd '{path}' | awk -F'\\t' '{awk}'", shell=True, check=True, capture_output=True,
                         text=True).stdout
    rows = []
    for line in out.strip().split("\n"):
        if line:
            t = line.split("\t")
            rows.append((int(t[1]), int(t[2]), t[3]))
    return pd.DataFrame(rows, columns=["start", "end", "info"])


def read_trv(path, chrom, start, end):
    header = subprocess.run(f"gzip -cd '{path}' | sed -n '1p'", shell=True, check=True, capture_output=True,
                            text=True).stdout.rstrip("\n").split("\t")
    out = subprocess.run(["tabix", path, f"{chrom}:{start + 1}-{end}"], check=True, capture_output=True, text=True).stdout
    rows = [line.split("\t") for line in out.strip().split("\n") if line]
    d = pd.DataFrame(rows, columns=header) if rows else pd.DataFrame(columns=header)
    if d.empty:
        return d
    d["start"] = d["start"].astype(int)
    d["end"] = d["end"].astype(int)
    d["n_alt"] = d["AC"].map(lambda s: sum(int(x) > 0 for x in s.split(",") if x not in (".", "")))
    d = d[(d["FILTER"] == "PASS") & (d["n_alt"] > 0)]
    return d


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--catalog", nargs="+", required=True, help="LABEL=path; LABEL must equal the TRV SOURCE value")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--gene")
    g.add_argument("--region", help="chr:start-end (1-based, inclusive)")
    g.add_argument("--gene-table", help="trv_gene_table.py output; plot every gene passing --count-column >= --min-count")
    p.add_argument("--count-column", default="n_coding")
    p.add_argument("--min-count", type=int, default=10)
    p.add_argument("--flank", type=int, default=2000, help="bp added on each side of a gene (default 2000)")
    p.add_argument("--max-motif", type=int, help="Only draw loci whose shortest motif is <= this (e.g. 6 for STRs)")
    p.add_argument("--label", default="")
    p.add_argument("--out", required=True, help="PDF (multi-page in --gene-table mode)")
    p.add_argument("--out-dir", help="--gene-table mode: also write one PDF per gene here")
    args = p.parse_args()

    if args.gene_table:
        t = pd.read_csv(args.gene_table, sep="\t")
        sel = t[t[args.count_column] >= args.min_count].sort_values(args.count_column, ascending=False)
        print(f"{len(sel)} genes with {args.count_column} >= {args.min_count}: " +
              ", ".join(f"{r.gene} ({getattr(r, args.count_column)})" for r in sel.itertuples()))
        if args.out_dir:
            os.makedirs(args.out_dir, exist_ok=True)
        with PdfPages(args.out) as pdf:
            for r in sel.itertuples():
                chrom, gs, ge = gene_span(args.gtf, r.gene)
                note = f"{args.count_column} = {getattr(r, args.count_column)}"
                fig = plot_locus(args, chrom, max(gs - args.flank, 0), ge + args.flank,
                                 f"{r.gene} ({chrom}:{gs + 1:,}-{ge:,}, ±{args.flank:,} bp; {note})")
                pdf.savefig(fig, bbox_inches="tight")
                if args.out_dir:
                    fig.savefig(os.path.join(args.out_dir, f"{r.gene}.trv_locus.pdf"), bbox_inches="tight")
                plt.close(fig)
        return
    if args.gene:
        chrom, gs, ge = gene_span(args.gtf, args.gene)
        start, end = max(gs - args.flank, 0), ge + args.flank
        title = f"{args.gene} ({chrom}:{gs + 1:,}-{ge:,}, ±{args.flank:,} bp)"
    else:
        chrom, rng = args.region.split(":")
        a, b = (int(x.replace(",", "")) for x in rng.split("-"))
        start, end = a - 1, b
        title = f"{chrom}:{a:,}-{b:,}"
    fig = plot_locus(args, chrom, start, end, title)
    fig.savefig(args.out, bbox_inches="tight")
    plt.close(fig)


def plot_locus(args, chrom, start, end, title):
    trv = read_trv(args.trv, chrom, start, end)
    if args.max_motif and not trv.empty:
        trv = trv[trv["motifs"].map(lambda m: min(len(x) for x in m.split(","))) <= args.max_motif]
    cats = [c.split("=", 1) for c in args.catalog]
    genes = read_gtf_window(args.gtf, chrom, start, end)

    gene_names = list(dict.fromkeys(genes.loc[genes.feature == "transcript", "gene"])) if not genes.empty else []
    gene_rows = []  # stack overlapping genes
    for name in gene_names:
        tx = genes[(genes.gene == name) & (genes.feature == "transcript")].iloc[0]
        for r, last_end in enumerate(gene_rows):
            if tx.start > last_end:
                gene_rows[r] = tx.end
                break
        else:
            gene_rows.append(tx.end)
    n_gene_rows = max(len(gene_rows), 1)

    heights = [0.55] * len(cats) + [1.6, 0.6 * n_gene_rows + 0.3]
    fig, axes = plt.subplots(len(heights), 1, figsize=(12, sum(heights) + 1.2), sharex=True,
                             gridspec_kw={"height_ratios": heights, "hspace": 0.15})
    summary = []
    for i, (label, path) in enumerate(cats):
        ax = axes[i]
        color = CAT_COLORS[i % len(CAT_COLORS)]
        cat = read_catalog(path, chrom, start, end)
        if args.max_motif and not cat.empty:
            mot = cat["info"].str.extract(r"MOTIFS=([^;]+)")[0].fillna("")
            cat = cat[mot.map(lambda m: min((len(x) for x in m.split(",") if x), default=999)) <= args.max_motif]
        obs = trv[trv["SOURCE"] == label] if not trv.empty else trv
        obs_spans = set(zip(obs["start"], obs["end"])) if not obs.empty else set()
        n_obs_loci = 0
        for r in cat.itertuples():
            seen = (r.start, r.end) in obs_spans
            n_obs_loci += seen
            ax.add_patch(Rectangle((r.start, 0.15), max(r.end - r.start, (end - start) / 1500), 0.7,
                                   facecolor=color if seen else UNOBSERVED, edgecolor="none"))
        n_alt = int(obs["n_alt"].sum()) if not obs.empty else 0
        summary.append((label, len(cat), n_obs_loci, n_alt))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(label, rotation=0, ha="right", va="center", fontsize=9, color=INK)
        ax.text(1.01, 0.5, f"{n_alt:,} ALT alleles\n{n_obs_loci}/{len(cat)} loci observed", transform=ax.transAxes,
                ha="left", va="center", fontsize=8.5, color=INK)
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color(GRID)
    # observed TRVs (lollipops)
    ax = axes[len(cats)]
    for i, (label, _) in enumerate(cats):
        obs = trv[trv["SOURCE"] == label] if not trv.empty else trv
        if obs.empty:
            continue
        mid = (obs["start"] + obs["end"]) / 2
        ax.vlines(mid, 0, obs["n_alt"], color=CAT_COLORS[i % len(CAT_COLORS)], linewidth=1)
        ax.scatter(mid, obs["n_alt"], s=14, color=CAT_COLORS[i % len(CAT_COLORS)], zorder=3, label=label)
    ax.set_ylabel("ALT alleles\nper TRV", fontsize=9, color=INK)
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    if not trv.empty:
        ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
    # gene structure
    ax = axes[-1]
    row_of, rows_end = {}, []
    for name in gene_names:
        tx = genes[(genes.gene == name) & (genes.feature == "transcript")].iloc[0]
        for r, last_end in enumerate(rows_end):
            if tx.start > last_end:
                rows_end[r] = tx.end
                row_of[name] = r
                break
        else:
            rows_end.append(tx.end)
            row_of[name] = len(rows_end) - 1
    for name in gene_names:
        gsub = genes[genes.gene == name]
        tx = gsub[gsub.feature == "transcript"].iloc[0]
        y = n_gene_rows - row_of[name] - 0.5
        ax.hlines(y, max(tx.start, start), min(tx.end, end), color=GENE_COLOR, linewidth=1)
        n_arrows = 8
        step = (min(tx.end, end) - max(tx.start, start)) / (n_arrows + 1)
        for k in range(1, n_arrows + 1):
            x = max(tx.start, start) + k * step
            ax.annotate("", xy=(x + (1 if tx.strand == "+" else -1) * step * 0.15, y), xytext=(x, y),
                        arrowprops=dict(arrowstyle="->", color=GENE_COLOR, lw=0.8))
        for ex in gsub[gsub.feature == "exon"].itertuples():
            ax.add_patch(Rectangle((ex.start, y - 0.12), ex.end - ex.start, 0.24, facecolor=GENE_COLOR, edgecolor="none"))
        for cds in gsub[gsub.feature == "CDS"].itertuples():
            ax.add_patch(Rectangle((cds.start, y - 0.25), cds.end - cds.start, 0.5, facecolor=GENE_COLOR, edgecolor="none"))
        lx = min(max((tx.start + tx.end) / 2, start), end)
        ax.text(lx, y + 0.32, f"{name} ({tx.strand})", ha="center", va="bottom", fontsize=8.5, color=INK, style="italic")
    ax.set_ylim(0, n_gene_rows + 0.1)
    ax.set_yticks([])
    ax.set_ylabel("Genes", rotation=0, ha="right", va="center", fontsize=9, color=INK)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.set_xlim(start, end)
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_xlabel(f"{chrom} position", fontsize=9, color=INK)
    for a in axes:
        a.tick_params(colors=INK2, labelsize=8)
    motif_note = f"; loci with motif ≤ {args.max_motif} bp" if args.max_motif else ""
    fig.suptitle(f"{args.label + ': ' if args.label else ''}TR loci across {title}{motif_note}\n"
                 f"coloured = non-ref PASS TRV observed; grey = catalog locus not observed", fontsize=10, color=INK)
    fig.subplots_adjust(left=0.1, right=0.84, top=0.9)
    for label, n_cat, n_obs, n_alt in summary:
        print(f"{label}: {n_cat} catalog loci in window, {n_obs} observed non-ref, {n_alt} ALT alleles")
    print(f"genes in window: {', '.join(gene_names) if gene_names else 'none'}; TRVs plotted: {len(trv)}")
    return fig


if __name__ == "__main__":
    main()
