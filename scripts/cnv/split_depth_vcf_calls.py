#!/usr/bin/env python3
"""Split a merged gCNV depth VCF (DepthPreprocessing merged_vcf) into one BED row per (call, carrier), so that each
sample's call is refined and filtered on that sample's own depth.

Only DEL and DUP records with END - POS + 1 >= --min-size are written; a record is dropped when its interval overlaps
any --gaps interval (assembly gaps, centromeres, telomeres). Carriers are samples whose GT holds a non-reference allele.

Output BED (gzipped, with header): #chrom start end SVID carriers SVTYPE SVLEN orig_id, where start = POS - 1,
end = END, SVID = <record ID>__<sample> and carriers = that sample.

Usage:
  split_depth_vcf_calls.py --vcf merged_cnvs.vcf.gz --out-bed calls.per_sample.bed.gz [--gaps gaps.bed] [--min-size 10000]
"""

import argparse
import bisect
import gzip
import sys
from collections import defaultdict

import pysam


def load_intervals(path):
    iv = defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or line.startswith(("#", "track", "browser")) or not f[1].isdigit():
                continue
            iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in iv.items():
        v.sort()
        out = []
        for s, e in v:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[c] = out
    return merged


def overlaps(iv, c, s, e):
    v = iv.get(c, [])
    i = bisect.bisect_right(v, (e, e)) - 1
    while i >= 0 and v[i][1] > s:
        if v[i][0] < e:
            return True
        i -= 1
    return False


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="merged gCNV depth VCF")
    p.add_argument("--out-bed", required=True, help="output BED (.gz)")
    p.add_argument("--gaps", help="BED of assembly gaps, centromeres and telomeres; overlapping calls are dropped")
    p.add_argument("--min-size", type=int, default=10000, help="minimum END - POS + 1 (default 10000)")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    gaps = load_intervals(a.gaps) if a.gaps else {}
    n_rec = n_small = n_gap = n_rows = 0
    with pysam.VariantFile(a.vcf) as vcf, gzip.open(a.out_bed, "wt") as out:
        out.write("#chrom\tstart\tend\tSVID\tcarriers\tSVTYPE\tSVLEN\torig_id\n")
        for rec in vcf:
            svtype = rec.info.get("SVTYPE")
            if svtype not in ("DEL", "DUP"):
                continue
            n_rec += 1
            start, end = rec.pos - 1, rec.stop
            if end - start < a.min_size:
                n_small += 1
                continue
            if overlaps(gaps, rec.chrom, start, end):
                n_gap += 1
                continue
            for s, smp in rec.samples.items():
                gt = smp.get("GT")
                if gt and any(x not in (None, 0) for x in gt):
                    out.write(f"{rec.chrom}\t{start}\t{end}\t{rec.id}__{s}\t{s}\t{svtype}\t{end - start}\t{rec.id}\n")
                    n_rows += 1
    print(
        f"{n_rec} DEL/DUP records: {n_small} < {a.min_size} bp, {n_gap} overlapping gaps; {n_rows} sample-calls written",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
