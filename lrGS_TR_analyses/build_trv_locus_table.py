#!/usr/bin/env python3
"""Per-locus TRV table: size, motif, sequence composition, variability and genic context.

Inputs: TRV bed (gnomAD_LR.<cohort>.vep_parsed.annotated.TRV.bed.gz: #CHROM, ID, REF, ALT, FILTER, SOURCE, TRID, AC, AN),
the genic-context table from annotate_TRV_genic_context.sh (ID, genic_context, genes) and a GENCODE GTF (CDS blocks).

One row per TRV site with columns:
  ID chrom start end FILTER source n_components motif_len motifs size copies gc purity
  n_alt nonref_AF max_abs_delta n_alt_inframe n_alt_frameshift AF_inframe AF_frameshift genic_context genes
Definitions:
  start/end      TRID span (0-based half-open; min start / max end across components); size = end - start
  motif_len      shortest motif across TRID components; copies = size / motif_len
  gc             G+C fraction of the REF repeat sequence (REF trimmed to the TRID span)
  purity         fraction of REF repeat bases covered by a greedy exact tiling with the TRID motifs
  n_alt          distinct ALT alleles with AC > 0;  nonref_AF = sum(AC) / AN
  max_abs_delta  max |len(ALT) - len(REF)| over ALT alleles with AC > 0 (0 if none)
  inframe / frameshift  ALT alleles with AC > 0 and delta != 0 multiple of 3 / delta not multiple of 3
  genic_context  coding (TR span inside one merged CDS block) / coding_partial (overlaps CDS but crosses a boundary) /
                 UTR / intronic / intergenic (from the genic-context table)

Usage:
  build_trv_locus_table.py --trv-bed TRV.bed.gz --genic-tsv genic.tsv.gz --gtf gencode.gtf.gz --out loci.tsv.gz
"""
import argparse
import gzip

import numpy as np
import pandas as pd


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


def purity(seq, motifs):
    if not seq:
        return np.nan
    i = covered = 0
    while i < len(seq):
        for m in motifs:
            if seq.startswith(m, i):
                covered += len(m)
                i += len(m)
                break
        else:
            i += 1
    return covered / len(seq)


def locus_rows(ch):
    out = []
    cols = [ch[c] for c in ["#CHROM", "ID", "REF", "ALT", "FILTER", "SOURCE", "TRID", "AC", "AN"]]
    for chrom, vid, ref, alt, filt, src, trid, ac, an in zip(*cols):
        comps = [c.split("-") for c in trid.split(",")]
        s0 = min(int(c[1]) for c in comps)
        e0 = max(int(c[2]) for c in comps)
        motifs = sorted({c[-1] for c in comps}, key=lambda m: (-len(m), m))
        pos0 = int(vid.split("-")[1]) - 1
        left = max(s0 - pos0, 0)
        right = max(pos0 + len(ref) - e0, 0)
        rep = ref[left:len(ref) - right]
        acs = np.array(ac.split(","), dtype=int)
        delta = np.array([len(a) - len(ref) for a in alt.split(",")])
        carried = acs > 0
        inframe = carried & (delta != 0) & (delta % 3 == 0)
        frameshift = carried & (delta % 3 != 0)
        an = int(an)
        mlen = min(len(m) for m in motifs)
        out.append((vid, chrom, s0, e0, filt, src, len(comps), mlen, ",".join(motifs), e0 - s0, (e0 - s0) / mlen,
                    (rep.count("G") + rep.count("C")) / len(rep) if rep else np.nan, purity(rep, motifs),
                    int(carried.sum()), acs.sum() / an, int(np.abs(delta[carried]).max()) if carried.any() else 0,
                    int(inframe.sum()), int(frameshift.sum()), acs[inframe].sum() / an, acs[frameshift].sum() / an))
    return out


COLUMNS = ["ID", "chrom", "start", "end", "FILTER", "source", "n_components", "motif_len", "motifs", "size", "copies", "gc",
           "purity", "n_alt", "nonref_AF", "max_abs_delta", "n_alt_inframe", "n_alt_frameshift", "AF_inframe", "AF_frameshift"]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trv-bed", required=True)
    p.add_argument("--genic-tsv", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    rows = []
    for ch in pd.read_csv(args.trv_bed, sep="\t", dtype=str, chunksize=200000,
                          usecols=["#CHROM", "ID", "REF", "ALT", "FILTER", "SOURCE", "TRID", "AC", "AN"]):
        rows.extend(locus_rows(ch))
    d = pd.DataFrame(rows, columns=COLUMNS)
    ctx = pd.read_csv(args.genic_tsv, sep="\t", usecols=["ID", "genic_context", "genes"], dtype=str)
    d = d.merge(ctx, on="ID", how="left")
    assert d["genic_context"].notna().all()
    inside = cds_contained(d["chrom"].to_numpy(), d["start"].to_numpy(), d["end"].to_numpy(), cds_blocks(args.gtf))
    d.loc[(d["genic_context"] == "coding") & ~inside, "genic_context"] = "coding_partial"
    d.to_csv(args.out, sep="\t", index=False, float_format="%.5g", compression="gzip")
    print(f"{args.out}: {len(d):,} loci")
    print(d.groupby(["genic_context"]).size().to_string())


if __name__ == "__main__":
    main()
