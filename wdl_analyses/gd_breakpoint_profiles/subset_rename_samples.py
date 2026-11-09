#!/usr/bin/env python3
"""Keep only the samples in an ID map and rename them, in a VCF or a wide per-sample table (BAF TSV, binned depth matrix).

ID map: 2 columns, tab- or space-separated: real sample ID, new sample ID (a header line is skipped if its first field is
not a sample in the input). Output samples follow the order of the map. Samples in the map but missing from the input
stop the run (listed) unless --allow-missing.

VCF (.vcf, .vcf.gz, .bcf; needs bcftools on PATH):
  bcftools view -s <real IDs> | bcftools reheader -s <map>; FORMAT fields other than --keep-format (default GT,DP,AD)
  are removed (--keep-format all keeps every field); --drop-info also removes all INFO fields. Output is bgzipped and
  indexed (.tbi) when it ends in .vcf.gz.
Table (anything else, plain or .gz): the first --fixed-cols columns are kept as is (5 for a BAF table: CHROM POS ID REF
  ALT; 3 for a depth matrix: #Chr Start End), the rest are sample columns.

Usage:
  subset_rename_samples.py --input calls.vcf.gz --id-map id_map.tsv --out calls.renamed.vcf.gz [--keep-format GT,DP,AD]
  subset_rename_samples.py --input chr16_bp4.baf.tsv.gz --id-map id_map.tsv --out chr16_bp4.renamed.baf.tsv.gz --fixed-cols 5
  subset_rename_samples.py --input region.RD.txt.gz --id-map id_map.tsv --out region.renamed.RD.txt.gz --fixed-cols 3
"""
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile


def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") or path.endswith(".bgz") else open(path, mode.replace("t", ""))


def read_map(path):
    pairs = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 2 and not line.startswith("#"):
                pairs.append((f[0], f[1]))
    return pairs


def check(pairs, available, allow_missing, what):
    if pairs and pairs[0][0] not in available and pairs[0][0].lower() in ("sample", "sample_id", "id", "real_id"):
        pairs = pairs[1:]  # header line
    missing = [r for r, _ in pairs if r not in available]
    if missing and not allow_missing:
        sys.exit(f"{len(missing)} of {len(pairs)} mapped samples are not in the {what}: " + ", ".join(missing[:20])
                 + (" ..." if len(missing) > 20 else "") + "  (use --allow-missing to skip them)")
    kept = [(r, n) for r, n in pairs if r in available]
    dup = {n for _, n in kept if [m for _, m in kept].count(n) > 1}
    if dup:
        sys.exit("new IDs used more than once: " + ", ".join(sorted(dup)))
    print(f"keeping {len(kept)} samples ({len(missing)} mapped but absent)", file=sys.stderr)
    return kept


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def do_vcf(a, pairs):
    if not shutil.which("bcftools"):
        sys.exit("bcftools not found on PATH")
    samples = set(subprocess.run(["bcftools", "query", "-l", a.input], capture_output=True, text=True, check=True).stdout.split())
    kept = check(pairs, samples, a.allow_missing, "VCF")
    with tempfile.TemporaryDirectory() as tmp:
        keep = os.path.join(tmp, "keep.txt")
        ren = os.path.join(tmp, "rename.txt")
        with open(keep, "w") as k, open(ren, "w") as r:
            for real, new in kept:
                k.write(real + "\n")
                r.write(f"{real} {new}\n")
        strip = []
        if a.keep_format != "all":
            strip.append("^" + ",".join("FORMAT/" + x for x in a.keep_format.split(",")))
        if a.drop_info:
            strip.append("INFO")
        annot = f"| bcftools annotate -x '{','.join(strip)}' -Ou " if strip else ""
        otype = "z" if a.out.endswith(".gz") else ("b" if a.out.endswith(".bcf") else "v")
        run(f"bcftools view -S {keep} --force-samples -Ou '{a.input}' {annot}"
            f"| bcftools reheader -s {ren} | bcftools view -O{otype} -o '{a.out}'")
    if otype in "zb":
        run(f"bcftools index -f {'-t ' if otype == 'z' else ''}'{a.out}'")
    print(f"wrote {a.out}", file=sys.stderr)


def do_table(a, pairs):
    with opener(a.input) as fh, opener(a.out, "wt") as out:
        header = fh.readline().rstrip("\n").split("\t")
        fixed, samples = header[:a.fixed_cols], header[a.fixed_cols:]
        kept = check(pairs, set(samples), a.allow_missing, "table header")
        idx = [a.fixed_cols + samples.index(r) for r, _ in kept]
        out.write("\t".join(fixed + [n for _, n in kept]) + "\n")
        n = 0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.write("\t".join(f[:a.fixed_cols] + [f[i] for i in idx]) + "\n")
            n += 1
    print(f"wrote {a.out}: {n} rows x {len(kept)} samples", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="VCF/BCF or wide per-sample table")
    p.add_argument("--id-map", required=True, help="2 columns: real sample ID, new sample ID")
    p.add_argument("--out", required=True)
    p.add_argument("--fixed-cols", type=int, default=5, help="leading non-sample columns of a table (5 BAF, 3 depth)")
    p.add_argument("--keep-format", default="GT,DP,AD", help="VCF FORMAT fields to keep, or 'all'")
    p.add_argument("--drop-info", action="store_true", help="VCF: remove all INFO fields")
    p.add_argument("--allow-missing", action="store_true", help="skip mapped samples absent from the input")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    pairs = read_map(a.id_map)
    if any(a.input.endswith(x) for x in (".vcf", ".vcf.gz", ".vcf.bgz", ".bcf")):
        do_vcf(a, pairs)
    else:
        do_table(a, pairs)


if __name__ == "__main__":
    main()
