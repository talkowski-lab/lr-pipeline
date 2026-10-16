#!/usr/bin/env python3
"""Tandem repeats in one or more genes and their lrGS vs srGS cohort coverage.

For each gene, finds TR catalog loci overlapping the gene (+/- flank), and
reports for each locus and each 100bp bin:
  - srGS (gnomAD v3) mean depth, normalized depth, low5_frac (fraction of
    samples not over 5x)
  - lrGS (HPRC+HGSVC) mean depth, normalized depth, low5_frac (1 - n_gt5 / N)
Depth is normalized to the genome-wide mean of the same chromosome class
(autosomes, or chrX with --sr-norm-x/--lr-norm-x).

TR class uses the dominant motif: for compound loci (TRExplorer variation
clusters) the component motif spanning the most bp; STR if its length is
<= --max-str-motif, else VNTR.

Bins with srGS norm < --sr-max and lrGS norm >= --lr-min are shaded.

Outputs: <prefix>.bins.tsv, <prefix>.TR_loci.tsv, <prefix>.pdf (one page per gene)
"""
import argparse
import gzip

import matplotlib.ticker
import pysam
from matplotlib import font_manager
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

LR_COLOR = "#2a78d6"
SR_COLOR = "#eb6834"
FLAG_SHADE = "#fbe1d6"
COV_KEYS = ["sr_mean", "sr_norm", "sr_low5_frac", "lr_mean", "lr_norm", "lr_low5_frac"]


def gene_models(gtf, genes, transcripts):
    """Gene interval and exons of one transcript per gene: the transcript given
    in `transcripts` if any, else Ensembl_canonical, else APPRIS principal."""
    wanted = set(genes)
    ivs, tx_exons, tx_rank = {}, {}, {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] not in ("gene", "exon") or 'gene_name "' not in c[8]:
                continue
            gene = c[8].split('gene_name "')[1].split('"')[0]
            if gene not in wanted:
                continue
            if c[2] == "gene":
                ivs[gene] = (c[0], int(c[3]) - 1, int(c[4]), c[6])
            else:
                tx = c[8].split('transcript_name "')[1].split('"')[0]
                tx_exons.setdefault(gene, {}).setdefault(tx, []).append((int(c[3]) - 1, int(c[4])))
                tx_rank.setdefault(gene, {})[tx] = 0 if "Ensembl_canonical" in c[8] else (
                    1 if "appris_principal" in c[8] else 2)
    models = {}
    for gene in genes:
        if gene not in ivs:
            print(f"WARNING: {gene} not found in {gtf}; skipped")
            continue
        tx = transcripts.get(gene) or min(tx_rank[gene], key=lambda t: (tx_rank[gene][t], t))
        models[gene] = (ivs[gene], tx, sorted(tx_exons[gene][tx]))
    return models


def dominant_motif(info):
    fields = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
    bp_by_motif = {}
    for comp in fields["ID"].split(","):
        parts = comp.split("-")
        if len(parts) >= 4 and parts[1].isdigit() and parts[2].isdigit():
            bp_by_motif[parts[3]] = bp_by_motif.get(parts[3], 0) + int(parts[2]) - int(parts[1])
    if not bp_by_motif:
        return fields["MOTIFS"].split(",")[0], fields["MOTIFS"]
    return max(bp_by_motif, key=bp_by_motif.get), fields["MOTIFS"]


def catalog_loci(specs, windows, max_str):
    """windows: {gene: (chrom, lo, hi)}; returns {gene: [locus, ...]}."""
    by_chrom = {}
    for gene, (chrom, lo, hi) in windows.items():
        by_chrom.setdefault(chrom, []).append((lo, hi, gene))
    out = {gene: [] for gene in windows}
    for spec in specs:
        label, path = spec.split("=", 1)
        with gzip.open(path, "rt") as f:
            for line in f:
                c = line.split("\t", 4)
                if c[0] not in by_chrom:
                    continue
                s, e = int(c[1]), int(c[2])
                hits = [g for lo, hi, g in by_chrom[c[0]] if e > lo and s < hi]
                if not hits:
                    continue
                motif, motifs = dominant_motif(c[3].rstrip("\n"))
                for g in hits:
                    out[g].append({
                        "gene": g, "chrom": c[0], "start": s, "end": e, "catalog": label,
                        "class": "STR" if len(motif) <= max_str else "VNTR",
                        "dominant_motif": motif, "motif_len": len(motif), "all_motifs": motifs,
                    })
    for g in out:
        out[g].sort(key=lambda x: (x["start"], x["end"]))
    return out


def read_bins(sr_tbx, lr_tbx, gene, chrom, start, end, sr_norm, lr_norm, lr_n):
    sr = {}
    for r in sr_tbx.fetch(chrom, start, end):
        c = r.split("\t")
        sr[int(c[1])] = (int(c[2]), float(c[3]), int(c[4]), float(c[5]))
    lr = {}
    for r in lr_tbx.fetch(chrom, start, end):
        c = r.split("\t")
        lr[int(c[1])] = (float(c[3]), int(c[6]))
    bins = []
    for s in sorted(sr):
        e, sr_mean, n_loci, sr_low5 = sr[s]
        lr_mean, lr_gt5 = lr.get(s, (float("nan"), 0))
        bins.append({
            "gene": gene, "chrom": chrom, "start": s, "end": e,
            "sr_mean": sr_mean, "sr_norm": sr_mean / sr_norm, "sr_low5_frac": sr_low5, "sr_n_loci": n_loci,
            "lr_mean": lr_mean, "lr_norm": lr_mean / lr_norm, "lr_low5_frac": 1 - lr_gt5 / lr_n,
        })
    return bins


def locus_coverage(locus, bins):
    tot, acc = 0, dict.fromkeys(COV_KEYS, 0.0)
    for b in bins:
        ov = min(locus["end"], b["end"]) - max(locus["start"], b["start"])
        if ov > 0:
            tot += ov
            for k in COV_KEYS:
                acc[k] += b[k] * ov
    return {k: acc[k] / tot if tot else float("nan") for k in COV_KEYS}


def write_tsv(path, rows, cols):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")


def plot_gene(pdf, gene, gene_iv, tx, exons, bins, loci, lo, hi, sr_max, lr_min):
    fig, axes = plt.subplots(3, 1, figsize=(10, 6.5), sharex=True,
                             gridspec_kw={"height_ratios": [3, 2, 1.6]})
    xs = [b["start"] / 1e6 for b in bins] + [bins[-1]["end"] / 1e6]
    flagged = [b for b in bins if b["sr_norm"] < sr_max and b["lr_norm"] >= lr_min]
    for ax, key, ylab in ((axes[0], "norm", "Depth / genome-wide mean"),
                          (axes[1], "low5_frac", "Fraction of samples < 5x")):
        for b in flagged:
            ax.axvspan(b["start"] / 1e6, b["end"] / 1e6, color=FLAG_SHADE, linewidth=0, zorder=0)
        for pfx, col, lab in (("lr", LR_COLOR, "lrGS (HPRC+HGSVC)"), ("sr", SR_COLOR, "srGS (gnomAD v3)")):
            ys = [b[f"{pfx}_{key}"] for b in bins]
            ax.stairs(ys, xs, color=col, linewidth=2 if len(bins) < 400 else 1, label=lab, baseline=None)
        ax.set_ylabel(ylab, fontsize=10)
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].axhline(1, color="#9a9a9a", linewidth=1, linestyle=(0, (3, 3)))
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylim(-0.02, 1.02)
    handles, labels = axes[0].get_legend_handles_labels()
    if flagged:
        handles.append(plt.Rectangle((0, 0), 1, 1, color=FLAG_SHADE))
        labels.append(f"srGS < {sr_max}x and lrGS >= {lr_min}x")
    axes[0].legend(handles, labels, frameon=False, fontsize=9, loc="lower left", ncol=3)

    ax = axes[2]
    ax.plot([gene_iv[1] / 1e6, gene_iv[2] / 1e6], [2, 2], color="#555555", linewidth=1)
    for s, e in exons:
        ax.add_patch(plt.Rectangle((s / 1e6, 1.7), (e - s) / 1e6, 0.6, color="#555555", linewidth=0))
    ax.text(gene_iv[1] / 1e6, 2.55, f"{gene} ({gene_iv[3]}, {tx})", fontsize=9, va="bottom")
    rows = {"STR": 1.0, "VNTR": 0.3}
    min_w = (hi - lo) / 600
    for lc in loci:
        y = rows[lc["class"]]
        ax.add_patch(plt.Rectangle((lc["start"] / 1e6, y - 0.18), max(lc["end"] - lc["start"], min_w) / 1e6, 0.36,
                                   color="#222222" if lc["class"] == "VNTR" else "#8a8a8a", linewidth=0))
        if lc["class"] == "VNTR" and len(loci) <= 40:
            ax.text((lc["start"] + lc["end"]) / 2e6, y - 0.32, f"{lc['motif_len']} bp motif",
                    fontsize=7.5, ha="center", va="top")
    ax.set_yticks([rows["STR"], rows["VNTR"], 2])
    ax.set_yticklabels(["STR", "VNTR", "exons"], fontsize=9)
    ax.set_ylim(-0.4, 3.0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(lo / 1e6, hi / 1e6)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.3f"))
    ax.set_xlabel(f"{gene_iv[0]} position (Mb, GRCh38)", fontsize=10)
    fig.suptitle(f"{gene}: lrGS vs srGS coverage over tandem repeats (100 bp bins)", fontsize=11, x=0.01, ha="left")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--gene", nargs="+", help="gene name(s)")
    g.add_argument("--gene-list", help="file with one gene name per line")
    ap.add_argument("--transcript", nargs="*", default=[], help="GENE=transcript_name overrides for the exon track")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--sr-bins", required=True, help="tabix-indexed srGS 100bp bins (extract_chrom_coverage_bins.py)")
    ap.add_argument("--lr-summary", required=True, help="tabix-indexed lrGS RD summary bed (mean in col 4, n_gt5 col 7)")
    ap.add_argument("--lr-n-samples", type=int, required=True)
    ap.add_argument("--sr-norm", type=float, required=True, help="srGS genome-wide autosomal mean depth")
    ap.add_argument("--lr-norm", type=float, required=True, help="lrGS genome-wide autosomal mean depth")
    ap.add_argument("--sr-norm-x", type=float, help="srGS chrX mean depth (default: --sr-norm)")
    ap.add_argument("--lr-norm-x", type=float, help="lrGS chrX mean depth (default: --lr-norm)")
    ap.add_argument("--catalog", nargs="+", required=True, help="LABEL=TRGT-style catalog .bed.gz")
    ap.add_argument("--max-str-motif", type=int, default=6)
    ap.add_argument("--sr-max", type=float, default=0.5)
    ap.add_argument("--lr-min", type=float, default=0.7)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    if "Arial" not in {f.name for f in font_manager.fontManager.ttflist}:
        raise SystemExit("Arial font not available to matplotlib")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42

    if args.gene_list:
        with open(args.gene_list) as f:
            genes = [x.strip() for x in f if x.strip()]
    else:
        genes = args.gene
    transcripts = dict(x.split("=", 1) for x in args.transcript)
    models = gene_models(args.gtf, genes, transcripts)
    windows = {}
    for gene, (iv, _, _) in models.items():
        windows[gene] = (iv[0], (iv[1] - args.flank) // 100 * 100, -(-(iv[2] + args.flank) // 100) * 100)
    loci_by_gene = catalog_loci(args.catalog, windows, args.max_str_motif)

    sr_tbx, lr_tbx = pysam.TabixFile(args.sr_bins), pysam.TabixFile(args.lr_summary)
    all_bins, all_loci = [], []
    with PdfPages(f"{args.out_prefix}.pdf") as pdf:
        for gene in genes:
            if gene not in models:
                continue
            iv, tx, exons = models[gene]
            chrom, lo, hi = windows[gene]
            is_x = chrom == "chrX"
            sr_norm = args.sr_norm_x if is_x and args.sr_norm_x else args.sr_norm
            lr_norm = args.lr_norm_x if is_x and args.lr_norm_x else args.lr_norm
            bins = read_bins(sr_tbx, lr_tbx, gene, chrom, lo, hi, sr_norm, lr_norm, args.lr_n_samples)
            loci = loci_by_gene[gene]
            for lc in loci:
                lc["transcript"] = tx
                lc["in_exon"] = "yes" if any(lc["start"] < e and lc["end"] > s for s, e in exons) else "no"
                lc.update(locus_coverage(lc, bins))
            all_bins += bins
            all_loci += loci
            if bins:
                plot_gene(pdf, gene, iv, tx, exons, bins, loci, lo, hi, args.sr_max, args.lr_min)

    write_tsv(f"{args.out_prefix}.bins.tsv", all_bins, ["gene", "chrom", "start", "end", "sr_n_loci"] + COV_KEYS)
    write_tsv(f"{args.out_prefix}.TR_loci.tsv", all_loci,
              ["gene", "chrom", "start", "end", "catalog", "class", "dominant_motif", "motif_len", "all_motifs",
               "transcript", "in_exon"] + COV_KEYS)
    print(f"{len(models)} genes plotted -> {args.out_prefix}.pdf")


if __name__ == "__main__":
    main()
