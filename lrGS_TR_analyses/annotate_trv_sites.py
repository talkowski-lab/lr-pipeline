#!/usr/bin/env python3
"""Add intergenic-region and gene-constraint annotations to the AnalyzeTRVariants site BED.

Inputs:
  --sites        <prefix>.TRV.sites.bed.gz from AnalyzeTRVariants (columns include #chrom start end ID genic_context genes)
  --intergenic   intergenic regions BED with header: chr start end size is_desert ... size_stratum
                 (e.g. gene_dessert_refs/intergenic_regions_annotated.bed.gz)
  --gene-info    gene table with gene_name, LOEUF_v4, LOEUF_v4_percentile, LOEUF_v4_decile
                 (e.g. gene_information_loeufv4_pHaplo_pTriplo.with_GATKSV_blacklist.txt.gz)
Added columns (appended; "." when not applicable or not found):
  intergenic TRVs   intergenic_region (chr:start-end of the region with the largest overlap with the TR span),
                    intergenic_region_size, intergenic_category (size_stratum), intergenic_is_desert,
                    intergenic_n_regions (regions overlapped). TRVs outside every listed region get intergenic_category
                    "between_regions" (inside a gene in the region file's gene model, which is broader than the MANE /
                    canonical transcripts used for genic_context) or "chrom_not_in_region_file" (e.g. chrX / chrY)
  coding / intronic LOEUF_v4, LOEUF_v4_percentile, LOEUF_v4_decile: one value per gene in `genes`, same order, comma-joined
                    ("NA" for a gene in the table without a score, "." for a gene not in the table)
Output: bgzipped BED (tabix-indexed when tabix is on PATH) and a summary printed.

Usage:
  annotate_trv_sites.py --sites gnomAD_LR.hprc_hgsvc.TRV.sites.bed.gz --intergenic intergenic_regions_annotated.bed.gz \
      --gene-info gene_information_loeufv4_pHaplo_pTriplo.with_GATKSV_blacklist.txt.gz \
      --out gnomAD_LR.hprc_hgsvc.TRV.sites.annotated.bed.gz
"""
import argparse
import bisect
import collections
import gzip
import shutil
import subprocess

import pandas as pd

NEW_COLS = ["intergenic_region", "intergenic_region_size", "intergenic_category", "intergenic_is_desert",
            "intergenic_n_regions", "LOEUF_v4", "LOEUF_v4_percentile", "LOEUF_v4_decile"]
GENE_CONTEXTS = {"coding", "intronic"}


def load_regions(path):
    r = pd.read_csv(path, sep="\t", dtype={"chr": str})
    idx = {}
    for chrom, s in r.sort_values(["chr", "start"]).groupby("chr"):
        rows = list(zip(s["start"], s["end"], s["size"], s["size_stratum"], s["is_desert"].astype(str)))
        idx[chrom] = ([x[0] for x in rows], rows, max(x[1] - x[0] for x in rows))
    return idx, len(r)


def overlap_regions(idx, chrom, s, e):
    if chrom not in idx:
        return []
    starts, rows, maxlen = idx[chrom]
    out = []
    j = bisect.bisect_left(starts, e) - 1
    while j >= 0 and starts[j] > s - maxlen:
        rs, re_, size, cat, desert = rows[j]
        ov = min(e, re_) - max(s, rs)
        if ov > 0 or (e == s and rs <= s < re_):
            out.append((max(ov, 0), rs, re_, size, cat, desert))
        j -= 1
    return out


def fmt(v):
    if pd.isna(v):
        return "NA"
    return f"{v:g}" if isinstance(v, float) else str(v)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", required=True)
    p.add_argument("--intergenic", required=True)
    p.add_argument("--gene-info", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    regions, n_regions = load_regions(args.intergenic)
    gi = pd.read_csv(args.gene_info, sep="\t", usecols=["gene_name", "LOEUF_v4", "LOEUF_v4_percentile", "LOEUF_v4_decile"],
                     na_values=["NA", ""])
    gi = gi.drop_duplicates("gene_name").set_index("gene_name")
    genes = {g: (fmt(r.LOEUF_v4), fmt(r.LOEUF_v4_percentile), fmt(r.LOEUF_v4_decile)) for g, r in gi.iterrows()}

    stats = collections.Counter()
    cat_counts = collections.Counter()
    gene_seen, gene_found, gene_scored = set(), set(), set()
    tmp = args.out + ".tmp"
    with gzip.open(args.sites, "rt") as fin, open(tmp, "w") as fout:
        header = fin.readline().rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        fout.write("\t".join(header + NEW_COLS) + "\n")
        for line in fin:
            t = line.rstrip("\n").split("\t")
            ctx = t[col["genic_context"]]
            add = ["."] * len(NEW_COLS)
            if ctx == "intergenic":
                stats["intergenic"] += 1
                hits = overlap_regions(regions, t[0], int(t[1]), int(t[2]))
                add[4] = str(len(hits))
                if hits:
                    best = max(hits, key=lambda h: h[0])
                    add[0:4] = [f"{t[0]}:{best[1]}-{best[2]}", str(best[3]), best[4], best[5]]
                    cat_counts[best[4]] += 1
                    stats["intergenic_in_region"] += 1
                    if len(hits) > 1:
                        stats["intergenic_multi_region"] += 1
                else:
                    add[2] = "between_regions" if t[0] in regions else "chrom_not_in_region_file"
                    cat_counts[add[2]] += 1
            elif ctx in GENE_CONTEXTS:
                stats[ctx] += 1
                gl = t[col["genes"]].split(",")
                vals = [genes.get(g, (".", ".", ".")) for g in gl]
                add[5:8] = [",".join(v[k] for v in vals) for k in range(3)]
                for g, v in zip(gl, vals):
                    gene_seen.add(g)
                    if v[0] != ".":
                        gene_found.add(g)
                    if v[0] not in (".", "NA"):
                        gene_scored.add(g)
                if all(v[0] not in (".", "NA") for v in vals):
                    stats[f"{ctx}_all_genes_scored"] += 1
            fout.write("\t".join(t + add) + "\n")

    bgzip = shutil.which("bgzip")
    if bgzip:
        subprocess.run(f"{bgzip} -c {tmp} > {args.out}", shell=True, check=True)
        if shutil.which("tabix"):
            subprocess.run(["tabix", "-f", "-p", "bed", "-S", "1", args.out], check=True)
    else:
        with open(tmp, "rb") as fi, gzip.open(args.out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    subprocess.run(["rm", "-f", tmp], check=True)

    print(f"intergenic regions loaded: {n_regions:,}; gene table rows: {len(gi):,}")
    print(f"intergenic TRVs: {stats['intergenic']:,}; in a listed region: {stats['intergenic_in_region']:,} "
          f"({100 * stats['intergenic_in_region'] / max(stats['intergenic'], 1):.2f}%); "
          f"overlapping 2+ regions: {stats['intergenic_multi_region']:,}")
    for cat, n in sorted(cat_counts.items(), key=lambda x: -x[1]):
        print(f"  {cat}: {n:,}")
    for ctx in ("coding", "intronic"):
        print(f"{ctx} TRVs: {stats[ctx]:,}; with LOEUF_v4 for every gene: {stats[ctx + '_all_genes_scored']:,} "
              f"({100 * stats[ctx + '_all_genes_scored'] / max(stats[ctx], 1):.1f}%)")
    print(f"genes in coding/intronic TRVs: {len(gene_seen):,}; in gene table: {len(gene_found):,}; "
          f"with LOEUF_v4: {len(gene_scored):,}")


if __name__ == "__main__":
    main()
