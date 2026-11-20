#!/usr/bin/env python3
"""Write the inserted sequences of a remap TSV as FASTA (ALT minus anchor base, named by ID) and a table of the
reference flanks around each insertion point (for TSD / poly-A checks), using a local reference FASTA (.fai index).

Flanks: left = the --flank bp ending at POS (the anchor base is the last base), right = the --flank bp after POS.
"""
import argparse

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--remap", required=True, help="Truvari remap TSV, no header (CHROM POS REF ALT ID ...)")
    ap.add_argument("--ref-fa", required=True, help="Reference FASTA with .fai (e.g. hg38 chr22)")
    ap.add_argument("--flank", type=int, default=50, help="Flank length, bp [50]")
    ap.add_argument("--min-len", type=int, default=1, help="Skip inserted sequences shorter than this [1]")
    ap.add_argument("--out-fa", required=True, help="Output FASTA of inserted sequences")
    ap.add_argument("--out-flanks", required=True, help="Output TSV: ID CHROM POS len left_flank right_flank")
    a = ap.parse_args()

    ref = pysam.FastaFile(a.ref_fa)
    n = 0
    with open(a.remap) as f, open(a.out_fa, "w") as fa, open(a.out_flanks, "w") as fl:
        fl.write("ID\tCHROM\tPOS\tlen\tleft_flank\tright_flank\n")
        for line in f:
            c = line.rstrip("\n").split("\t")
            chrom, pos, ref_allele, alt, iid = c[0], int(c[1]), c[2], c[3].upper(), c[4]
            ins = alt[len(ref_allele):] if alt.startswith(ref_allele.upper()) else alt[1:]
            if len(ins) < a.min_len:
                continue
            left = ref.fetch(chrom, max(0, pos - a.flank), pos).upper()
            right = ref.fetch(chrom, pos, pos + a.flank).upper()
            fa.write(f">{iid}\n")
            fa.writelines(ins[i:i + 80] + "\n" for i in range(0, len(ins), 80))
            fl.write(f"{iid}\t{chrom}\t{pos}\t{len(ins)}\t{left}\t{right}\n")
            n += 1
    print(f"wrote {n} inserted sequences")


if __name__ == "__main__":
    main()
