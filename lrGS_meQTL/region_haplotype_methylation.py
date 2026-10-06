#!/usr/bin/env python3
"""Summarise per-haplotype CpG methylation within target regions.

Subcommands:
  contig   for one contig of a wide per-haplotype methylation table (#chrom start end <sample>_hap1 ... ; percent
           methylation, '.' = missing), take every target region on that contig with >= --min-cpg CpG sites and
           write, per haplotype column, the number of called CpGs, mean, median and SD (SD needs >= 2 called CpGs)
  concat   concatenate gzipped TSVs that share a header (header written once)

contig outputs (<prefix> = --prefix):
  <prefix>.regions.tsv.gz    region_index (0-based row of the region in --regions, header excluded), the region's
                             original columns, region_n_cpg (CpG sites in the methylation table within
                             [start, end)), region_n_haplotypes_called (haplotypes with >= 1 called CpG in the region)
  <prefix>.n_called.tsv.gz, <prefix>.mean.tsv.gz, <prefix>.median.tsv.gz, <prefix>.sd.tsv.gz
                             region_index, chrom, start, end, then one column per haplotype (NA = not computable)
"""
import argparse
import gzip
import warnings

import numpy as np
import pandas as pd


def load_matrix(path, chunksize=200000):
    """Return (chrom, CpG start positions, haplotype column names, float32 matrix CpGs x haplotypes)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
    names = header[3:]
    dtypes = {c: np.float32 for c in names}
    dtypes.update({header[0]: str, header[1]: np.int64, header[2]: np.int64})
    chroms, pos, blocks = [], [], []
    reader = pd.read_csv(path, sep="\t", na_values=["."], dtype=dtypes, chunksize=chunksize, engine="c")
    for chunk in reader:
        chroms.append(chunk.iloc[:, 0].to_numpy())
        pos.append(chunk.iloc[:, 1].to_numpy(np.int64))
        blocks.append(chunk.iloc[:, 3:].to_numpy(np.float32))
    if not blocks:
        return np.array([], dtype=str), np.array([], dtype=np.int64), names, np.zeros((0, len(names)), np.float32)
    pos = np.concatenate(pos)
    mat = np.concatenate(blocks)
    order = np.argsort(pos, kind="stable")
    if not np.all(order == np.arange(len(pos))):
        pos, mat = pos[order], mat[order]
    return np.concatenate(chroms), pos, names, mat


def load_regions(path):
    """Return the region table (all columns as text) with a region_index column, and its column names."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = fh.readline().rstrip("\n")
    has_header = first.startswith("#")
    regions = pd.read_csv(path, sep="\t", header=None, skiprows=1 if has_header else 0, dtype=str,
                          keep_default_na=False)
    cols = first.lstrip("#").split("\t") if has_header else \
        ["chrom", "start", "end"] + [f"col{i}" for i in range(4, regions.shape[1] + 1)]
    regions.columns = cols
    regions.insert(0, "region_index", np.arange(len(regions)))
    return regions, cols


def cmd_contig(args):
    regions, cols = load_regions(args.regions)
    chrom_col, start_col, end_col = cols[0], cols[1], cols[2]
    regions = regions[regions[chrom_col] == args.contig].copy()
    _, pos, names, mat = load_matrix(args.bed)
    starts = regions[start_col].astype(np.int64).to_numpy()
    ends = regions[end_col].astype(np.int64).to_numpy()
    i0 = np.searchsorted(pos, starts, side="left")
    i1 = np.searchsorted(pos, ends, side="left")
    n_cpg = i1 - i0
    keep = n_cpg >= args.min_cpg
    regions, i0, i1, n_cpg = regions[keep], i0[keep], i1[keep], n_cpg[keep]

    n_hap = len(names)
    out = {k: np.full((len(regions), n_hap), np.nan, dtype=np.float32) for k in ("n_called", "mean", "median", "sd")}
    cache = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for r, (a, b) in enumerate(zip(i0, i1)):
            key = (a, b)
            if key not in cache:
                block = mat[a:b]
                n = (~np.isnan(block)).sum(axis=0).astype(np.float32)
                sd = np.nanstd(block, axis=0, ddof=1)
                sd[n < 2] = np.nan
                cache = {key: (n, np.nanmean(block, axis=0), np.nanmedian(block, axis=0), sd)}
            n, mean, median, sd = cache[key]
            out["n_called"][r], out["mean"][r], out["median"][r], out["sd"][r] = n, mean, median, sd

    regions["region_n_cpg"] = n_cpg
    regions["region_n_haplotypes_called"] = (out["n_called"] > 0).sum(axis=1)
    regions.rename(columns={chrom_col: "#" + chrom_col}).to_csv(
        f"{args.prefix}.regions.tsv.gz", sep="\t", index=False, compression="gzip")
    key_cols = regions[["region_index", chrom_col, start_col, end_col]].reset_index(drop=True)
    key_cols.columns = ["region_index", "#chrom", "start", "end"]
    for stat, arr in out.items():
        df = pd.concat([key_cols, pd.DataFrame(arr, columns=names)], axis=1)
        df.to_csv(f"{args.prefix}.{stat}.tsv.gz", sep="\t", index=False, na_rep="NA",
                  float_format="%.0f" if stat == "n_called" else "%.2f", compression="gzip")
    print(f"{args.contig}: {int(keep.sum()):,} of {len(keep):,} regions with >= {args.min_cpg} CpG sites; "
          f"{len(pos):,} CpG sites x {n_hap} haplotypes")


def cmd_concat(args):
    with gzip.open(args.out, "wt") as out:
        wrote_header = False
        for path in args.inputs:
            with gzip.open(path, "rt") as fh:
                header = fh.readline()
                if not wrote_header:
                    out.write(header)
                    wrote_header = True
                for line in fh:
                    out.write(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("contig")
    c.add_argument("--bed", required=True, help="per-contig wide per-haplotype methylation table (bgzip)")
    c.add_argument("--regions", required=True, help="target regions (BED-like, optional '#' header, optional gzip)")
    c.add_argument("--contig", required=True)
    c.add_argument("--min-cpg", type=int, default=1)
    c.add_argument("--prefix", required=True)
    c.set_defaults(func=cmd_contig)
    m = sub.add_parser("concat")
    m.add_argument("--inputs", nargs="+", required=True)
    m.add_argument("--out", required=True)
    m.set_defaults(func=cmd_concat)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
