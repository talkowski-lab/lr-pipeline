#!/usr/bin/env python3
"""Look up precomputed CADD scores (GRCh38 v1.7) for a set of variants without downloading the score files: the
bgzipped, tabix-indexed files are read remotely over https with pysam, fetching only the blocks around the variants.

SNVs are looked up in whole_genome_SNVs.tsv.gz (every possible SNV), other variants (indels, MNVs) in
gnomad.genomes.r4.0.indel.tsv.gz (indels seen in gnomAD v4 genomes); variants absent from the files get NA. Variants are
grouped into regions (consecutive variants within --gap bp), so each part of a file is fetched once. Alleles are matched
as given (left-aligned, VCF-style); a 'chr' prefix is dropped for the lookup.

Input: a VCF/BCF (multi-allelic records split per ALT) or a TSV with a header containing chrom, pos, ref, alt (any
other columns are carried over).
Output: TSV with the input columns (VCF: chrom pos id ref alt) + cadd_raw cadd_phred cadd_source.

Region mode (--region, repeatable, instead of --input): every score in the region(s), from both files:
  chrom pos ref alt cadd_raw cadd_phred cadd_source (snv: all 3 possible SNVs per base; indel: gnomAD v4 indels).

Usage:
  cadd_remote_lookup.py --input variants.vcf.gz --out variants.cadd.tsv [--gap 10000] [--cadd-url URL]
  cadd_remote_lookup.py --input variants.tsv --out variants.cadd.tsv
  cadd_remote_lookup.py --region chr16:29324176-29638676 [--region ...] --out bp4.cadd.tsv [--snv-only | --indel-only]
"""
import argparse
import csv
import gzip
import sys

import pysam

BASE = "https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38"


def read_variants(path):
    if path.endswith((".vcf", ".vcf.gz", ".vcf.bgz", ".bcf")):
        vcf = pysam.VariantFile(path)
        cols = ["chrom", "pos", "id", "ref", "alt"]
        rows = []
        for rec in vcf:
            for alt in rec.alts or ():
                rows.append(dict(chrom=rec.chrom, pos=rec.pos, id=rec.id or ".", ref=rec.ref, alt=alt))
        return cols, rows
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        cols = r.fieldnames
        need = {"chrom", "pos", "ref", "alt"}
        if not need <= set(cols):
            sys.exit(f"TSV needs columns {sorted(need)}; has {cols}")
        rows = [dict(row, pos=int(row["pos"])) for row in r]
    return cols, rows


def regions(rows, gap):
    """Groups of row indices: same chromosome, consecutive positions within gap."""
    order = sorted(range(len(rows)), key=lambda i: (rows[i]["chrom"], rows[i]["pos"]))
    out, cur = [], []
    for i in order:
        if cur and (rows[i]["chrom"] != rows[cur[-1]]["chrom"] or rows[i]["pos"] - rows[cur[-1]]["pos"] > gap):
            out.append(cur)
            cur = []
        cur.append(i)
    if cur:
        out.append(cur)
    return out


def dump_regions(a):
    """Every CADD score in the regions, SNV file then indel file, sorted by position."""
    kinds = [k for k in ("snv", "indel") if not (k == "snv" and a.indel_only) and not (k == "indel" and a.snv_only)]
    urls = {"snv": a.snv_url, "indel": a.indel_url}
    n = 0
    with open(a.out, "w") as fh:
        fh.write("chrom\tpos\tref\talt\tcadd_raw\tcadd_phred\tcadd_source\n")
        for reg in a.region:
            c, rng = reg.rsplit(":", 1)
            s, e = (int(x.replace(",", "")) for x in rng.split("-"))
            cc = c[3:] if c.startswith("chr") else c
            lines = []
            for k in kinds:
                for line in pysam.TabixFile(urls[k]).fetch(cc, s - 1, e):
                    f = line.split("\t")
                    lines.append((int(f[1]), f"{c}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}\t{f[5]}\t{k}\n"))
            for _, l in sorted(lines, key=lambda x: x[0]):
                fh.write(l)
            n += len(lines)
            print(f"{reg}: {len(lines):,} scores", file=sys.stderr)
    print(f"{n:,} scores -> {a.out}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", help="VCF/BCF or TSV of variants to score")
    p.add_argument("--region", action="append", help="chrom:start-end: write every CADD score in it (repeatable)")
    p.add_argument("--snv-only", action="store_true")
    p.add_argument("--indel-only", action="store_true")
    p.add_argument("--out", required=True)
    p.add_argument("--gap", type=int, default=10000, help="variants closer than this are fetched as one region")
    p.add_argument("--snv-url", default=f"{BASE}/whole_genome_SNVs.tsv.gz")
    p.add_argument("--indel-url", default=f"{BASE}/gnomad.genomes.r4.0.indel.tsv.gz")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    if bool(a.input) == bool(a.region):
        sys.exit("give either --input or --region")
    if a.region:
        dump_regions(a)
        return
    cols, rows = read_variants(a.input)
    print(f"{len(rows)} variants", file=sys.stderr)
    tabs = {"snv": pysam.TabixFile(a.snv_url), "indel": pysam.TabixFile(a.indel_url)}
    for r in rows:
        r["kind"] = "snv" if len(r["ref"]) == 1 and len(r["alt"]) == 1 else "indel"
        r["cadd_raw"] = r["cadd_phred"] = "NA"
        r["cadd_source"] = "not_in_cadd_file"
    for kind in ("snv", "indel"):
        idx = [i for i, r in enumerate(rows) if r["kind"] == kind]
        sub = [rows[i] for i in idx]
        groups = regions(sub, a.gap)
        for k, g in enumerate(groups, 1):
            c = sub[g[0]]["chrom"]
            cc = c[3:] if c.startswith("chr") else c
            lo, hi = sub[g[0]]["pos"], sub[g[-1]]["pos"]
            want = {(sub[j]["pos"], sub[j]["ref"].upper(), sub[j]["alt"].upper()): j for j in g}
            try:
                for line in tabs[kind].fetch(cc, lo - 1, hi):
                    f = line.split("\t")
                    key = (int(f[1]), f[2], f[3])
                    if key in want:
                        r = sub[want[key]]
                        r["cadd_raw"], r["cadd_phred"], r["cadd_source"] = f[4], f[5], kind
            except ValueError:
                pass  # chromosome not in the CADD file (e.g. alt contigs)
            if k % 200 == 0:
                print(f"  {kind}: {k}/{len(groups)} regions", file=sys.stderr)
    out_cols = [c for c in cols] + ["cadd_raw", "cadd_phred", "cadd_source"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in out_cols) + "\n")
    n = sum(r["cadd_source"] != "not_in_cadd_file" for r in rows)
    print(f"{n}/{len(rows)} variants scored -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
