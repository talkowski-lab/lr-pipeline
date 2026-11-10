#!/usr/bin/env python3
"""Convert an SV site table (svtk-style BED: '#chrom start end name svtype ...', 0-based start, with header) to the
interval BED used for CADD summaries: chrom start end name svtype.

Insertions (svtype INS or INS:*, e.g. INS:ME:ALU) become a window around the insertion position POS (= start + 1):
POS - --ins-flank to POS + --ins-flank (2 x flank + 1 bp, whatever the length of the insertion's interval in the
table); all other SV types keep start and end.
Rows are written in input order; starts are clipped at 0.

Usage:
  sv_bed_for_cadd.py --input gnomAD_SV_v3.PASS.key_info.gz --out gnomAD_SV_v3.PASS.cadd_intervals.bed.gz [--ins-flank 10]
"""
import argparse
import gzip
import sys
from collections import Counter


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="SV BED table with a header line (columns chrom start end name svtype)")
    p.add_argument("--out", required=True, help="output BED (.gz for gzip)")
    p.add_argument("--ins-flank", type=int, default=10)
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    opener_in = gzip.open if a.input.endswith((".gz", ".bgz")) else open
    opener_out = gzip.open if a.out.endswith(".gz") else open
    n, types = 0, Counter()
    with opener_in(a.input, "rt") as fh, opener_out(a.out, "wt") as out:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        col = {c: i for i, c in enumerate(header)}
        for c in ("chrom", "start", "end", "name", "svtype"):
            if c not in col:
                sys.exit(f"column {c} missing from the header of {a.input}: {header[:8]}")
        out.write("#chrom\tstart\tend\tname\tsvtype\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, name, svtype = f[col["chrom"]], f[col["name"]], f[col["svtype"]]
            start, end = int(f[col["start"]]), int(f[col["end"]])
            if svtype == "INS" or svtype.startswith("INS:"):
                pos = start + 1  # VCF POS (BED start = POS - 1)
                start, end = max(0, pos - 1 - a.ins_flank), pos + a.ins_flank
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{svtype}\n")
            n += 1
            types[svtype] += 1
    print(f"{n:,} SVs -> {a.out}; " + ", ".join(f"{k} {v:,}" for k, v in types.most_common()), file=sys.stderr)


if __name__ == "__main__":
    main()
