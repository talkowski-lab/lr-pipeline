#!/usr/bin/env python3
"""Coding tandem repeats and their srGS vs lrGS coverage.

1. Coding TR regions: TR catalog loci overlapping a CDS of a protein_coding
   gene (any transcript). Loci overlapping across catalogs are merged into one
   region. Region class uses the dominant motif of its largest constituent
   locus (for compound loci, the component motif spanning the most bp):
   STR if <= --max-str-motif bp, else VNTR. Genes = all protein_coding genes
   whose CDS the region overlaps (multi-gene regions count for each gene).
2. Coverage per 100bp bin: depth normalized to the genome-wide mean of the
   same chromosome class (autosome or chrX). chrY is excluded.
3. A region is "srGS-poor / lrGS-ok" if any overlapping bin has
   srGS norm < --sr-max and lrGS norm >= --lr-min.
4. segdup_frac: fraction of the region inside --segdup; regions with
   segdup_frac >= 0.5 are labelled SegDup (low srGS there is usually
   multi-mapping rather than the repeat itself).

Outputs (<prefix>.*):
  coding_TR.tsv.gz       all coding TR regions with coverage, flag and SegDup overlap
  flagged.tsv            flagged regions (gene, region, lrGS / srGS coverage)
  flagged_genes.txt      genes with >= 1 flagged region
  size_table.tsv         size-range breakdown (regions, genes, flagged)
  size_summary.pdf       total coding TRs by size and % flagged by size
"""
import argparse
import bisect
import gzip

import matplotlib.ticker
import pysam
from matplotlib import font_manager
from matplotlib import pyplot as plt

TOTAL_COLOR = "#2a78d6"
FLAG_COLOR = "#eb6834"


def load_cds(gtf):
    cds = {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t", 9)
            if c[2] != "CDS" or 'gene_type "protein_coding"' not in c[8]:
                continue
            gene = c[8].split('gene_name "')[1].split('"')[0]
            cds.setdefault(c[0], []).append((int(c[3]) - 1, int(c[4]), gene))
    index = {}
    for chrom, ivs in cds.items():
        ivs.sort()
        starts = [s for s, _, _ in ivs]
        max_len = max(e - s for s, e, _ in ivs)
        index[chrom] = (starts, ivs, max_len)
    return index


def cds_genes(index, chrom, start, end):
    if chrom not in index:
        return set()
    starts, ivs, max_len = index[chrom]
    i = bisect.bisect_left(starts, start - max_len)
    genes = set()
    while i < len(ivs) and ivs[i][0] < end:
        s, e, g = ivs[i]
        if e > start:
            genes.add(g)
        i += 1
    return genes


def dominant_motif(info):
    fields = dict(kv.split("=", 1) for kv in info.split(";") if "=" in kv)
    bp_by_motif = {}
    for comp in fields["ID"].split(","):
        parts = comp.split("-")
        if len(parts) >= 4 and parts[1].isdigit() and parts[2].isdigit():
            bp_by_motif[parts[3]] = bp_by_motif.get(parts[3], 0) + int(parts[2]) - int(parts[1])
    if not bp_by_motif:
        return fields["MOTIFS"].split(",")[0]
    return max(bp_by_motif, key=bp_by_motif.get)


def coding_loci(specs, cds_index, skip_chroms):
    loci = []
    for spec in specs:
        label, path = spec.split("=", 1)
        with gzip.open(path, "rt") as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if c[0] in skip_chroms:
                    continue
                s, e = int(c[1]), int(c[2])
                if not cds_genes(cds_index, c[0], s, e):
                    continue
                loci.append((c[0], s, e, label, dominant_motif(c[3])))
    loci.sort()
    return loci


def merge_regions(loci, cds_index, max_str):
    regions, cur = [], None
    for chrom, s, e, label, motif in loci:
        if cur and chrom == cur["chrom"] and s < cur["end"]:
            cur["end"] = max(cur["end"], e)
            cur["members"].append((s, e, label, motif))
            continue
        if cur:
            regions.append(cur)
        cur = {"chrom": chrom, "start": s, "end": e, "members": [(s, e, label, motif)]}
    if cur:
        regions.append(cur)
    for r in regions:
        s, e, _, motif = max(r["members"], key=lambda m: m[1] - m[0])
        r["size"] = r["end"] - r["start"]
        r["dominant_motif"] = motif
        r["motif_len"] = len(motif)
        r["class"] = "STR" if len(motif) <= max_str else "VNTR"
        r["catalogs"] = ",".join(sorted({m[2] for m in r["members"]}))
        r["genes"] = ",".join(sorted(cds_genes(cds_index, r["chrom"], r["start"], r["end"])))
    return regions


def add_coverage(regions, sr_path, lr_path, norms, lr_n, sr_max, lr_min):
    sr_tbx, lr_tbx = pysam.TabixFile(sr_path), pysam.TabixFile(lr_path)
    for r in regions:
        sr_norm, lr_norm = norms[r["chrom"] == "chrX"]
        sr = {int(x.split("\t")[1]): x.split("\t") for x in sr_tbx.fetch(r["chrom"], r["start"], r["end"])}
        lr = {int(x.split("\t")[1]): x.split("\t") for x in lr_tbx.fetch(r["chrom"], r["start"], r["end"])}
        tot, acc, worst, flagged = 0, [0.0] * 6, None, False
        for b in sorted(sr):
            if b not in lr:
                continue
            ov = min(r["end"], int(sr[b][2])) - max(r["start"], b)
            vals = [float(sr[b][3]), float(sr[b][3]) / sr_norm, float(sr[b][5]),
                    float(lr[b][3]), float(lr[b][3]) / lr_norm, 1 - int(lr[b][6]) / lr_n]
            tot += ov
            acc = [a + v * ov for a, v in zip(acc, vals)]
            if worst is None or vals[1] < worst[2]:
                worst = [b] + vals
            if vals[1] < sr_max and vals[4] >= lr_min:
                flagged = True
        keys = ["sr_mean", "sr_norm", "sr_low5_frac", "lr_mean", "lr_norm", "lr_low5_frac"]
        if tot == 0:
            r.update({k: float("nan") for k in keys})
            r.update({f"worst_{k}": float("nan") for k in keys})
            r["worst_bin"], r["flag"] = ".", "no"
            continue
        r.update({k: a / tot for k, a in zip(keys, acc)})
        r.update({f"worst_{k}": v for k, v in zip(keys, worst[1:])})
        r["worst_bin"] = f"{r['chrom']}:{worst[0]}-{worst[0] + 100}"
        r["flag"] = "yes" if flagged else "no"


def load_intervals(path):
    ivs = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            c = line.split("\t")
            ivs.setdefault(c[0], []).append((int(c[1]), int(c[2])))
    index = {}
    for chrom, v in ivs.items():
        v.sort()
        merged = []
        for st, en in v:
            if merged and st <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], en)
            else:
                merged.append([st, en])
        index[chrom] = ([m[0] for m in merged], merged)
    return index


def overlap_bp(index, chrom, start, end):
    if chrom not in index:
        return 0
    starts, ivs = index[chrom]
    i = max(bisect.bisect_right(starts, start) - 1, 0)
    bp = 0
    while i < len(ivs) and ivs[i][0] < end:
        bp += max(0, min(end, ivs[i][1]) - max(start, ivs[i][0]))
        i += 1
    return bp


def write_tsv(path, rows, cols):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wt") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")


def size_label(lo, hi):
    return f">{lo - 1}" if hi is None else f"{lo}-{hi}"


def bin_regions(regions, breaks):
    """breaks: ascending lower bounds; returns list of (label, members)."""
    out = []
    for i, lo in enumerate(breaks):
        hi = breaks[i + 1] - 1 if i + 1 < len(breaks) else None
        members = [r for r in regions if r["size"] >= lo and (hi is None or r["size"] <= hi)]
        out.append((size_label(lo, hi), members))
    return out


def gene_set(rows):
    return {g for r in rows for g in r["genes"].split(",") if g}


def style_axes(ax):
    ax.grid(axis="y", color="#e5e5e5", linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def plot_size_summary(path, binned, crit):
    """Top: total coding TRs per size bin. Bottom: % of them srGS-poor / lrGS-ok."""
    labels = [lab for lab, _ in binned]
    total = [len(m) for _, m in binned]
    flag = [sum(r["flag"] == "yes" for r in m) for _, m in binned]
    prop = [100 * f / t if t else 0 for f, t in zip(flag, total)]
    xs = list(range(len(labels)))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 11), sharex=True, gridspec_kw={"height_ratios": [1, 1]})
    ax1.bar(xs, total, width=0.7, color=TOTAL_COLOR)
    for x, t in zip(xs, total):
        ax1.text(x, t, f"{t:,}", ha="center", va="bottom", fontsize=15)
    ax1.set_ylim(0, max(total) * 1.1)
    ax1.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:,.0f}"))
    ax1.set_ylabel("Coding TR regions", fontsize=20)
    ax1.set_title("Coding TRs by size", fontsize=22, loc="left")
    ax2.bar(xs, prop, width=0.7, color=FLAG_COLOR)
    for x, p_, f in zip(xs, prop, flag):
        ax2.text(x, p_, f"{f:,}", ha="center", va="bottom", fontsize=15)
    ax2.set_ylim(0, 50)
    ax2.set_ylabel("% poorly covered in srGS,\nwell covered in lrGS", fontsize=20)
    ax2.set_title(f"Proportion poorly covered by srGS only\n({crit})", fontsize=20, loc="left")
    ax2.set_xticks(xs)
    ax2.set_xticklabels(labels, rotation=45, ha="right", fontsize=18)
    ax2.set_xlabel("TR size in reference (bp)", fontsize=20)
    for ax in (ax1, ax2):
        style_axes(ax)
        ax.tick_params(axis="y", labelsize=18)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--catalog", nargs="+", required=True, help="LABEL=TRGT-style catalog .bed.gz")
    ap.add_argument("--sr-bins", required=True, help="tabix-indexed srGS 100bp bins")
    ap.add_argument("--lr-summary", required=True, help="tabix-indexed lrGS RD summary bed")
    ap.add_argument("--lr-n-samples", type=int, required=True)
    ap.add_argument("--sr-norm", type=float, required=True, help="srGS autosomal mean depth")
    ap.add_argument("--lr-norm", type=float, required=True, help="lrGS autosomal mean depth")
    ap.add_argument("--sr-norm-x", type=float, required=True, help="srGS chrX mean depth")
    ap.add_argument("--lr-norm-x", type=float, required=True, help="lrGS chrX mean depth")
    ap.add_argument("--segdup", required=True, help="SegDup bed.gz")
    ap.add_argument("--sr-max", type=float, default=0.5)
    ap.add_argument("--lr-min", type=float, default=0.7)
    ap.add_argument("--max-str-motif", type=int, default=6)
    ap.add_argument("--plot-breaks", default="1,11,21,31,51,76,101,151,201,301,501,1001,2001")
    ap.add_argument("--table-breaks", default="1,21,51,101,201,501,1001")
    ap.add_argument("--out-prefix", required=True)
    args = ap.parse_args()

    if "Arial" not in {f.name for f in font_manager.fontManager.ttflist}:
        raise SystemExit("Arial font not available to matplotlib")
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["pdf.fonttype"] = 42

    cds_index = load_cds(args.gtf)
    loci = coding_loci(args.catalog, cds_index, {"chrY", "chrM"})
    regions = merge_regions(loci, cds_index, args.max_str_motif)
    norms = {False: (args.sr_norm, args.lr_norm), True: (args.sr_norm_x, args.lr_norm_x)}
    add_coverage(regions, args.sr_bins, args.lr_summary, norms, args.lr_n_samples, args.sr_max, args.lr_min)
    segdup = load_intervals(args.segdup)
    for r in regions:
        r["segdup_frac"] = overlap_bp(segdup, r["chrom"], r["start"], r["end"]) / r["size"]
        r["context"] = "SegDup" if r["segdup_frac"] >= 0.5 else "unique"
    print(f"{len(loci)} coding catalog loci -> {len(regions)} merged coding TR regions")

    cov = ["sr_mean", "sr_norm", "sr_low5_frac", "lr_mean", "lr_norm", "lr_low5_frac"]
    cols = ["chrom", "start", "end", "size", "class", "dominant_motif", "motif_len", "catalogs", "genes"] + cov \
        + ["worst_bin"] + [f"worst_{k}" for k in cov] + ["segdup_frac", "context", "flag"]
    write_tsv(f"{args.out_prefix}.coding_TR.tsv.gz", regions, cols)
    flagged = [r for r in regions if r["flag"] == "yes"]
    flagged.sort(key=lambda r: (r["context"] == "SegDup", r["worst_sr_norm"]))
    write_tsv(f"{args.out_prefix}.flagged.tsv", flagged,
              ["genes", "chrom", "start", "end", "size", "class", "dominant_motif", "motif_len", "context",
               "worst_bin", "worst_lr_mean", "worst_lr_norm", "worst_sr_mean", "worst_sr_norm", "worst_sr_low5_frac",
               "lr_mean", "lr_norm", "sr_mean", "sr_norm", "sr_low5_frac"])
    with open(f"{args.out_prefix}.flagged_genes.txt", "w") as f:
        f.write("\n".join(sorted(gene_set(flagged))) + "\n")

    table_breaks = [int(x) for x in args.table_breaks.split(",")]
    with open(f"{args.out_prefix}.size_table.tsv", "w") as f:
        f.write("size_bp\tclass\tTR_regions\tgenes\tflagged_TR_regions\tflagged_genes\tpct_TR_flagged"
                "\tflagged_TR_unique\tflagged_genes_unique\tflagged_TR_SegDup\n")
        for cls in ("all", "STR", "VNTR"):
            sub = regions if cls == "all" else [r for r in regions if r["class"] == cls]
            for lab, m in bin_regions(sub, table_breaks) + [("all", sub)]:
                fl = [r for r in m if r["flag"] == "yes"]
                pct = 100 * len(fl) / len(m) if m else 0
                fu = [r for r in fl if r["context"] == "unique"]
                f.write(f"{lab}\t{cls}\t{len(m)}\t{len(gene_set(m))}\t{len(fl)}\t{len(gene_set(fl))}\t{pct:.2f}"
                        f"\t{len(fu)}\t{len(gene_set(fu))}\t{len(fl) - len(fu)}\n")

    binned = bin_regions(regions, [int(x) for x in args.plot_breaks.split(",")])
    crit = f"srGS < {args.sr_max}x and lrGS >= {args.lr_min}x of genome mean in any 100 bp bin"
    plot_size_summary(f"{args.out_prefix}.size_summary.pdf", binned, crit)


if __name__ == "__main__":
    main()
