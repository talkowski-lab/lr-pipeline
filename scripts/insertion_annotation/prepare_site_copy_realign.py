#!/usr/bin/env python3
"""Prepare the realignment of insertion pieces to the SD copy at the insertion site (last pipeline step).
For each insertion in --ids with AlignQueriesIterative pieces: when the site lies in an SD copy (build_sd_reference.py
loci) the target is the union of the site's SD copies, padded on each side by max(--min-flank, insertion length),
and every aligned piece (>= --min-piece bp) is written as a query.
Outputs in --out-dir: <ID>.insertion.fa, <ID>.right.fa / <ID>.left.fa (the reference right of POS and left of it,
insertion length + max(--min-pad, --pad-frac x length) each, for the direct tandem test), <ID>.pieces.fa
(>ID|piece|q_start|q_end), <ID>.target.fa (>chrom:start-end, 1-based) and
manifest.tsv (ID, site, len, site_copies, target, n_pieces). For insertions whose site is not in an SD
(site_copies "none") the target is the site +/- max(--min-flank, insertion length).
"""
import argparse
import collections
import csv
import os

import pysam


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True, help="File of insertion IDs")
    ap.add_argument("--segments", required=True, help="AlignQueriesIterative segments TSV (site reference)")
    ap.add_argument("--sequences", required=True, help="Inserted sequences FASTA")
    ap.add_argument("--loci", required=True, help="build_sd_reference.py <prefix>.loci.bed")
    ap.add_argument("--ref-fa", required=True, help="Reference FASTA with .fai (site contig)")
    ap.add_argument("--min-flank", type=int, default=10000, help="Minimum flank around the site copies [10000]")
    ap.add_argument("--min-piece", type=int, default=30, help="Minimum piece length [30]")
    ap.add_argument("--min-pad", type=int, default=1000, help="Extra flank beyond the insertion length, minimum [1000]")
    ap.add_argument("--pad-frac", type=float, default=0.1, help="Extra flank beyond the insertion length, x length [0.1]")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    ids = [x.strip() for x in open(a.ids) if x.strip()]
    seg = collections.defaultdict(list)
    for r in csv.DictReader(open(a.segments), delimiter="\t"):
        if r["ID"] in ids:
            seg[r["ID"]].append(r)
    seqs, name = {}, None
    for line in open(a.sequences):
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in seg:
                seqs[name] = []
        elif name in seqs:
            seqs[name].append(line.strip())
    loci = collections.defaultdict(list)
    with open(a.loci) as f:
        f.readline()
        for line in f:
            c = line.rstrip("\n").split("\t")
            loci[c[0]].append((int(c[1]), int(c[2]), c[5]))
    ref = pysam.FastaFile(a.ref_fa)
    with open(os.path.join(a.out_dir, "manifest.tsv"), "w") as o:
        o.write("ID\tsite\tlen\tsite_copies\ttarget\tn_pieces\n")
        for iid in ids:
            if iid not in seg:
                continue
            ch, pos = iid.split("-")[0], int(iid.split("-")[1])
            seq = "".join(seqs[iid]).upper()
            L = len(seq)
            sc = [x for x in loci.get(ch, []) if x[0] < pos <= x[1]]
            fl = max(a.min_flank, L)
            if sc:
                ts = max(1, min(x[0] for x in sc) + 1 - fl)
                te = min(ref.get_reference_length(ch), max(x[1] for x in sc) + fl)
            else:  # site not in an SD: the reference around the site
                ts, te = max(1, pos - fl), min(ref.get_reference_length(ch), pos + fl)
            with open(os.path.join(a.out_dir, f"{iid}.target.fa"), "w") as t:
                t.write(f">{ch}:{ts}-{te}\n{ref.fetch(ch, ts - 1, te).upper()}\n")
            with open(os.path.join(a.out_dir, f"{iid}.insertion.fa"), "w") as q:
                q.write(f">{iid}\n{seq}\n")
            # flank windows for the direct tandem test: insertion length + pad on each side of the site
            w = L + max(a.min_pad, int(a.pad_frac * L))
            rs, re_ = pos + 1, min(ref.get_reference_length(ch), pos + w)
            ls, le = max(1, pos - w + 1), pos
            with open(os.path.join(a.out_dir, f"{iid}.right.fa"), "w") as t:
                t.write(f">{ch}:{rs}-{re_}\n{ref.fetch(ch, rs - 1, re_).upper()}\n")
            with open(os.path.join(a.out_dir, f"{iid}.left.fa"), "w") as t:
                t.write(f">{ch}:{ls}-{le}\n{ref.fetch(ch, ls - 1, le).upper()}\n")
            n = 0
            with open(os.path.join(a.out_dir, f"{iid}.pieces.fa"), "w") as q:
                for r in sorted(seg[iid], key=lambda r: int(r["q_start"])):
                    qs, qe = int(r["q_start"]), int(r["q_end"])
                    if r["status"] == "aligned" and qe - qs + 1 >= a.min_piece:
                        n += 1
                        q.write(f">{iid}|{n}|{qs}|{qe}\n{seq[qs - 1:qe]}\n")
            o.write(f"{iid}\t{ch}:{pos}\t{L}\t{','.join(sorted({x[2] for x in sc})) or 'none'}\t{ch}:{ts}-{te}\t{n}\n")
    print(f"{len([i for i in ids if i in seg])} insertions with pieces written to {a.out_dir}")


if __name__ == "__main__":
    main()
