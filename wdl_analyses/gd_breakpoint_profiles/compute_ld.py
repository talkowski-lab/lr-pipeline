#!/usr/bin/env python3
"""Pairwise linkage disequilibrium (r^2) from unphased genotypes across a region, per site pair and binned for plotting.

Input: genotype table CHROM POS ID REF ALT + one GT column per sample (e.g. from extract_gt.sh), plain or .gz.
Sites kept: biallelic SNVs (single-base REF and ALT, no genotype with allele >= 2), call rate >= --min-call,
  minor allele frequency >= --min-maf, and no heterozygote excess (observed het / HWE-expected het <= --max-het-ratio;
  removes paralog (SD) differences that look heterozygous in nearly everyone); optionally only the sites in --sites
  (a table with CHROM POS columns, e.g. the clean-site BAF table of find_baf_terminal_loh.py).
Samples: all, minus those matching --exclude-pattern (regex; e.g. carriers of a deletion, whose hemizygous calls look hom).
r^2 = squared Pearson correlation of allele dosages (0/1/2) over samples called at both sites (composite LD, no phasing).

Outputs:
  <prefix>.sites.tsv      every input site with its QC values and whether it was kept
  <prefix>.pairs.tsv.gz   kept site pairs within --max-dist: pos1 pos2 dist r2 n
  <prefix>.binned.tsv     --bin-size windows: bin1 bin2 (window starts) mean_r2 max_r2 n_pairs, for every window pair

Usage:
  compute_ld.py --gt aou_II/chr16_bp4.renamed.GT.tsv.gz --out-prefix aou_II/chr16_bp4.renamed.ld \\
      [--exclude-pattern "GD_16p11.2-BP4-5"] [--sites clean_sites.baf.tsv.gz] [--min-maf 0.05] [--min-call 0.9]
      [--max-het-ratio 1.5] [--bin-size 5000] [--max-dist 0 (no limit)]
"""
import argparse
import re
import sys

import numpy as np
import pandas as pd


def dosage(gt):
    """GT string -> (dosage or nan, has allele >= 2)"""
    a = gt.replace("|", "/").split("/")
    if len(a) != 2 or "." in a:
        return np.nan, any(x not in (".", "0", "1") for x in a)
    if a[0] not in ("0", "1") or a[1] not in ("0", "1"):
        return np.nan, True
    return int(a[0]) + int(a[1]), False


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gt", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--exclude-pattern", help="regex; matching samples are left out")
    p.add_argument("--sites", help="table with CHROM and POS columns; only these sites are used")
    p.add_argument("--min-maf", type=float, default=0.05)
    p.add_argument("--min-call", type=float, default=0.9)
    p.add_argument("--max-het-ratio", type=float, default=1.5)
    p.add_argument("--bin-size", type=int, default=5000)
    p.add_argument("--max-dist", type=int, default=0, help="only pairs closer than this (0 = all pairs)")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()

    d = pd.read_csv(a.gt, sep="\t", dtype=str)
    samples = list(d.columns[5:])
    if a.exclude_pattern:
        samples = [s for s in samples if not re.search(a.exclude_pattern, s)]
    print(f"{len(d)} sites, {len(samples)} samples used", file=sys.stderr)

    g = d[samples].to_numpy()
    dos = np.full(g.shape, np.nan)
    multi = np.zeros(len(d), dtype=bool)
    for i in range(g.shape[0]):
        for j in range(g.shape[1]):
            v, m = dosage(g[i, j])
            dos[i, j] = v
            multi[i] |= m
    called = ~np.isnan(dos)
    call = called.mean(1)
    af = np.nansum(dos, 1) / np.maximum(2 * called.sum(1), 1)
    maf = np.minimum(af, 1 - af)
    het_obs = np.nansum(dos == 1, 1) / np.maximum(called.sum(1), 1)
    het_exp = 2 * af * (1 - af)
    het_ratio = np.where(het_exp > 0, het_obs / np.maximum(het_exp, 1e-9), np.nan)
    snv = (d["REF"].str.len().eq(1) & d["ALT"].str.len().eq(1)).to_numpy()
    keep = snv & ~multi & (call >= a.min_call) & (maf >= a.min_maf) & (het_ratio <= a.max_het_ratio)
    if a.sites:
        s = pd.read_csv(a.sites, sep="\t", usecols=["CHROM", "POS"], dtype=str)
        keep &= (d["CHROM"] + ":" + d["POS"]).isin(set(s["CHROM"] + ":" + s["POS"])).to_numpy()
    qc = d.iloc[:, :5].copy()
    qc["snv"], qc["multiallelic"], qc["call_rate"] = snv, multi, call.round(3)
    qc["maf"], qc["het_obs_exp_ratio"], qc["kept"] = maf.round(3), np.round(het_ratio, 2), keep
    qc.to_csv(a.out_prefix + ".sites.tsv", sep="\t", index=False)
    print(f"kept {keep.sum()} sites (SNV {snv.sum()}, not multiallelic {(snv & ~multi).sum()}, "
          f"het excess removed {(snv & ~multi & (het_ratio > a.max_het_ratio)).sum()})", file=sys.stderr)

    pos = d["POS"].astype(int).to_numpy()[keep]
    x = dos[keep]
    # r^2 with mean imputation of missing dosages (call rate is high), via standardised dosages
    mu = np.nanmean(x, 1, keepdims=True)
    x = np.where(np.isnan(x), mu, x) - mu
    sd = np.sqrt((x ** 2).sum(1, keepdims=True))
    z = x / np.where(sd > 0, sd, 1)
    r2 = (z @ z.T) ** 2
    n = called[keep].astype(float) @ called[keep].T.astype(float)

    iu = np.triu_indices(len(pos), k=1)
    dist = pos[iu[1]] - pos[iu[0]]
    sel = dist < a.max_dist if a.max_dist > 0 else np.ones(len(dist), dtype=bool)
    pairs = pd.DataFrame({"pos1": pos[iu[0]][sel], "pos2": pos[iu[1]][sel], "dist": dist[sel],
                          "r2": r2[iu][sel].round(4), "n": n[iu][sel].astype(int)})
    pairs.to_csv(a.out_prefix + ".pairs.tsv.gz", sep="\t", index=False, compression="gzip")

    b = (pos // a.bin_size) * a.bin_size
    pb = pd.DataFrame({"bin1": b[iu[0]][sel], "bin2": b[iu[1]][sel], "r2": r2[iu][sel]})
    binned = pb.groupby(["bin1", "bin2"])["r2"].agg(mean_r2="mean", max_r2="max", n_pairs="size").reset_index()
    binned.round(4).to_csv(a.out_prefix + ".binned.tsv", sep="\t", index=False)
    print(f"{len(pairs)} pairs; median r2 {np.median(pairs['r2']):.3f}; wrote {a.out_prefix}.{{sites,binned}}.tsv, "
          f".pairs.tsv.gz", file=sys.stderr)


if __name__ == "__main__":
    main()
