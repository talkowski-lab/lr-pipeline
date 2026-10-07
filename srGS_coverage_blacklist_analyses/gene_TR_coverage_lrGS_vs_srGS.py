#!/usr/bin/env python3
"""Tandem repeats in one gene and their lrGS vs srGS cohort coverage.

Finds TR catalog loci overlapping a gene (+/- flank), and reports for each
locus and each 100bp bin:
  - srGS (gnomAD v3) mean depth, normalized depth, low5_frac (fraction of
    samples not over 5x)
  - lrGS (HPRC+HGSVC) mean depth, normalized depth, low5_frac (1 - n_gt5 / N)

TR class uses the dominant motif: for compound loci (TRExplorer variation
clusters) the component motif spanning the most bp; STR if its length is
<= --max-str-motif, else VNTR.

Outputs: <prefix>.bins.tsv, <prefix>.TR_loci.tsv, <prefix>.pdf
"""
import argparse
import gzip

import matplotlib.ticker
import pysam
from matplotlib import font_manager
from matplotlib import pyplot as plt

LR_COLOR = "#2a78d6"
SR_COLOR = "#eb6834"


def gene_model(gtf, gene, transcript):
    """Gene interval and exons of one transcript: --transcript if given, else
    the Ensembl_canonical transcript, else the APPRIS principal one."""
    gene_iv, tx_exons, tx_rank = None, {}, {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if f'gene_name "{gene}";' not in c[8]:
                continue
            if c[2] == "gene":
                gene_iv = (c[0], int(c[3]) - 1, int(c[4]), c[6])
            elif c[2] == "exon":
                tx = c[8].split('transcript_name "')[1].split('"')[0]
                tx_exons.setdefault(tx, []).append((int(c[3]) - 1, int(c[4])))
                tx_rank[tx] = 0 if "Ensembl_canonical" in c[8] else (1 if "appris_principal" in c[8] else 2)
    if gene_iv is None:
        raise SystemExit(f"{gene} not found in {gtf}")
    if transcript is None:
        transcript = min(tx_rank, key=lambda t: (tx_rank[t], t))
    print(f"{gene}: exon track from {transcript}")
    return gene_iv, sorted(tx_exons[transcript])


def dominant_motif(info):
    fields = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
    bp_by_motif = {}
    for comp in fields["ID"].split(","):
        parts = comp.split("-")
        if len(parts) >= 4 and parts[1].isdigit() and parts[2].isdigit():
            bp_by_motif[parts[3]] = bp_by_motif.get(parts[3], 0) + int(parts[2]) - int(parts[1])
    if not bp_by_motif:
        m = fields["MOTIFS"].split(",")[0]
        return m, fields["MOTIFS"]
    return max(bp_by_motif, key=bp_by_motif.get), fields["MOTIFS"]


def catalog_loci(specs, chrom, start, end, max_str):
    loci = []
    for spec in specs:
        label, path = spec.split("=", 1)
        with gzip.open(path, "rt") as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if c[0] != chrom or int(c[2]) <= start or int(c[1]) >= end:
                    continue
                motif, motifs = dominant_motif(c[3])
                loci.append({
                    "chrom": c[0], "start": int(c[1]), "end": int(c[2]), "catalog": label,
                    "class": "STR" if len(motif) <= max_str else "VNTR",
                    "dominant_motif": motif, "motif_len": len(motif), "all_motifs": motifs,
                })
    return sorted(loci, key=lambda x: (x["start"], x["end"]))


def read_bins(sr_path, lr_path, chrom, start, end, sr_norm, lr_norm, lr_n):
    sr = {}
    for r in pysam.TabixFile(sr_path).fetch(chrom, start, end):
        c = r.split("\t")
        sr[int(c[1])] = (int(c[2]), float(c[3]), int(c[4]), float(c[5]))
    lr = {}
    for r in pysam.TabixFile(lr_path).fetch(chrom, start, end):
        c = r.split("\t")
        lr[int(c[1])] = (float(c[3]), int(c[6]))
    bins = []
    for s in sorted(sr):
        e, sr_mean, n_loci, sr_low5 = sr[s]
        lr_mean, lr_gt5 = lr.get(s, (float("nan"), 0))
        bins.append({
            "chrom": chrom, "start": s, "end": e,
            "sr_mean": sr_mean, "sr_norm": sr_mean / sr_norm, "sr_low5_frac": sr_low5, "sr_n_loci": n_loci,
            "lr_mean": lr_mean, "lr_norm": lr_mean / lr_norm, "lr_low5_frac": 1 - lr_gt5 / lr_n,
        })
    return bins


def locus_coverage(locus, bins):
    keys = ["sr_mean", "sr_norm", "sr_low5_frac", "lr_mean", "lr_norm", "lr_low5_frac"]
    tot, acc = 0, dict.fromkeys(keys, 0.0)
    for b in bins:
        ov = min(locus["end"], b["end"]) - max(locus["start"], b["start"])
        if ov > 0:
            tot += ov
            for k in keys:
                acc[k] += b[k] * ov
    return {k: acc[k] / tot for k in keys}


def write_tsv(path, rows, cols):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")


def plot(prefix, gene, gene_iv, exons, bins, loci, lo, hi):
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42
    fig, axes = plt.subplots(3, 1, figsize=(10, 6.5), sharex=True,
                             gridspec_kw={"height_ratios": [3, 2, 1.6]})
    xs = [b["start"] / 1e6 for b in bins] + [bins[-1]["end"] / 1e6]
    for ax, key, ylab in ((axes[0], "norm", "Depth / genome-wide mean"),
                          (axes[1], "low5_frac", "Fraction of samples < 5x")):
        for pfx, col, lab in (("lr", LR_COLOR, "lrGS (HPRC+HGSVC)"), ("sr", SR_COLOR, "srGS (gnomAD v3)")):
            ys = [b[f"{pfx}_{key}"] for b in bins]
            ax.stairs(ys, xs, color=col, linewidth=2, label=lab, baseline=None)
        ax.set_ylabel(ylab, fontsize=10)
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].axhline(1, color="#9a9a9a", linewidth=1, linestyle=(0, (3, 3)))
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylim(-0.02, 1.02)
    axes[0].legend(frameon=False, fontsize=9, loc="lower left", ncol=2)

    ax = axes[2]
    ax.plot([gene_iv[1] / 1e6, gene_iv[2] / 1e6], [2, 2], color="#555555", linewidth=1)
    for s, e in exons:
        ax.add_patch(plt.Rectangle((s / 1e6, 1.7), (e - s) / 1e6, 0.6, color="#555555", linewidth=0))
    ax.text(gene_iv[1] / 1e6, 2.55, f"{gene} ({gene_iv[3]})", fontsize=9, va="bottom")
    rows = {"STR": 1.0, "VNTR": 0.3}
    for lc in loci:
        y = rows[lc["class"]]
        ax.add_patch(plt.Rectangle((lc["start"] / 1e6, y - 0.18), max(lc["end"] - lc["start"], 15) / 1e6, 0.36,
                                   color="#222222" if lc["class"] == "VNTR" else "#8a8a8a", linewidth=0))
        if lc["class"] == "VNTR":
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
    fig.savefig(f"{prefix}.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--transcript", help="transcript_name for exon track (default: Ensembl_canonical)")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--sr-bins", required=True, help="tabix-indexed srGS 100bp bins (extract_chrom_coverage_bins.py)")
    ap.add_argument("--lr-summary", required=True, help="tabix-indexed lrGS RD summary bed (mean in col 4, n_gt5 col 7)")
    ap.add_argument("--lr-n-samples", type=int, required=True)
    ap.add_argument("--sr-norm", type=float, required=True, help="srGS genome-wide autosomal mean depth")
    ap.add_argument("--lr-norm", type=float, required=True, help="lrGS genome-wide autosomal mean depth")
    ap.add_argument("--catalog", nargs="+", required=True, help="LABEL=TRGT-style catalog .bed.gz")
    ap.add_argument("--max-str-motif", type=int, default=6)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    if "Arial" not in {f.name for f in font_manager.fontManager.ttflist}:
        raise SystemExit("Arial font not available to matplotlib")
    gene_iv, exons = gene_model(args.gtf, args.gene, args.transcript)
    chrom = gene_iv[0]
    lo = (gene_iv[1] - args.flank) // 100 * 100
    hi = -(-(gene_iv[2] + args.flank) // 100) * 100
    bins = read_bins(args.sr_bins, args.lr_summary, chrom, lo, hi, args.sr_norm, args.lr_norm, args.lr_n_samples)
    loci = catalog_loci(args.catalog, chrom, lo, hi, args.max_str_motif)
    for lc in loci:
        lc["in_exon"] = "yes" if any(lc["start"] < e and lc["end"] > s for s, e in exons) else "no"
        lc.update(locus_coverage(lc, bins))

    cov_cols = ["sr_mean", "sr_norm", "sr_low5_frac", "lr_mean", "lr_norm", "lr_low5_frac"]
    write_tsv(f"{args.out_prefix}.bins.tsv", bins,
              ["chrom", "start", "end", "sr_n_loci"] + cov_cols)
    write_tsv(f"{args.out_prefix}.TR_loci.tsv", loci,
              ["chrom", "start", "end", "catalog", "class", "dominant_motif", "motif_len", "all_motifs", "in_exon"] + cov_cols)
    plot(args.out_prefix, args.gene, gene_iv, exons, bins, loci, lo, hi)


if __name__ == "__main__":
    main()
