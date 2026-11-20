#!/usr/bin/env python3
"""Summarise BWA-MEM alignments (AlignQueriesBWA, bwa mem -a -Y, mapped records) of inserted sequences, one row per
insertion.
For every alignment: query span on the inserted sequence (clips from the CIGAR, flipped for reverse-strand records),
qcov = aligned query bases / insertion length, identity = matched bases / (M + I + D) with matched = M - (NM - I - D),
and the reference interval. The reported alignment: if any alignment covering >= --min-qcov lies within --adjacent bp
of the site (an insertion inside a tandem array has copies on both sides of the site, and the best-scoring one may be
a few units away), the best of those adjacent ones; otherwise the one with the most matched bases (qcov x identity;
ties: the one closest to the site). Its position relative to the insertion site (POS of the insertion's own contig):
  adjacent      the interval touches POS or ends within --adjacent bp of it (distance 0 when POS+1 is its first base,
                POS its last base, or POS lies inside it)
  spans site    POS lies inside the interval by more than --adjacent bp on both sides (insertion point inside the unit)
  <N> bp away   same contig, elsewhere
  other contig  another chromosome / unplaced contig
Call (best alignment covering >= --min-qcov of the insertion):
  tandem duplication        adjacent, or spans site (the insertion point lies inside the unit; relation = "spans site")
  interspersed duplication  same contig, more than --adjacent bp from the site
  interspersed (other contig)
otherwise "partial (<qcov>)" or "unmapped".
With --cytoband, origin_band = cytoband of the best copy's midpoint, and "acrocentric" is noted for p-arm
copies on chr13, chr14, chr15, chr21 and chr22.
Also reported: n_aln (all records), n_full (records with qcov >= --min-qcov), best mapq and the best record's
coordinates.
"""
import argparse
import gzip
import re

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def parse(c):
    ops = [(int(n), o) for n, o in CIG.findall(c[5])]
    lead = ops[0][0] if ops[0][1] in "SH" else 0
    tail = ops[-1][0] if ops[-1][1] in "SH" else 0
    m = sum(n for n, o in ops if o in "M=X")
    ins = sum(n for n, o in ops if o == "I")
    dl = sum(n for n, o in ops if o in "DN")
    qlen = sum(n for n, o in ops if o in "MIS=XH")
    nm = next((int(t[5:]) for t in c[11:] if t.startswith("NM:i:")), 0)
    rev = int(c[1]) & 16
    qs = tail if rev else lead
    qe = qlen - (lead if rev else tail)
    rs = int(c[3])
    re_ = rs + m + dl - 1
    ident = (m - (nm - ins - dl)) / max(1, m + ins + dl)
    return dict(chrom=c[2], rs=rs, re=re_, strand="-" if rev else "+", qs=qs, qe=qe, qlen=qlen, ident=ident,
                mapq=int(c[4]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sam", required=True, help="AlignQueriesBWA SAM (.gz)")
    ap.add_argument("--queries", required=True, help="Query FASTA (the sequences aligned)")
    ap.add_argument("--min-qcov", type=float, default=0.9, help="Min fraction of the insertion in one alignment [0.9]")
    ap.add_argument("--adjacent", type=int, default=5, help="Max site-to-copy distance for a tandem duplication, bp [5]")
    ap.add_argument("--cytoband", default=None, help="Optional UCSC cytoBand table (chrom start end name gieStain)")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    lens, name = {}, None
    for line in open(a.queries):
        if line.startswith(">"):
            name = line[1:].split()[0]
            lens[name] = 0
        else:
            lens[name] += len(line.strip())
    bands = {}
    if a.cytoband:
        for line in open(a.cytoband):
            c = line.rstrip("\n").split("\t")
            bands.setdefault(c[0], []).append((int(c[1]), int(c[2]), c[3]))

    def band_of(ch, p):
        nm = next((n for s, e, n in bands.get(ch, []) if s <= p < e), "")
        if not nm:
            return ""
        lab = ch.replace("chr", "") + nm
        return lab + (", acrocentric" if ch in ACRO and nm.startswith("p") else "")

    alns = {}
    with gzip.open(a.sam, "rt") as f:
        for line in f:
            if line.startswith("@"):
                continue
            c = line.rstrip("\n").split("\t")
            if int(c[1]) & 4:
                continue
            alns.setdefault(c[0], []).append(parse(c))
    with open(a.out, "w") as o:
        o.write("ID\tinsertion_site\taligned_location\tstrand\tlen\tcall\trelation\tdist_to_site\tbest_qcov\tbest_identity\t"
                "best_mapq\tbest_q\tn_aln\tn_full\torigin_band\n")
        for iid, L in lens.items():
            chrom, pos = iid.split("-")[0], int(iid.split("-")[1])
            al = alns.get(iid, [])
            if not al:
                o.write(f"{iid}\t{chrom}:{pos}\t\t\t{L}\tunmapped\t\t\t\t\t\t\t0\t0\t\n")
                continue

            def dist(x):
                if x["chrom"] != chrom:
                    return 10 ** 12
                return 0 if x["rs"] - 1 <= pos <= x["re"] else min(abs(x["rs"] - pos), abs(x["re"] - pos))

            for x in al:
                x["qcov"] = (x["qe"] - x["qs"]) / L
            adj = [x for x in al if x["qcov"] >= a.min_qcov and dist(x) <= a.adjacent]
            b = max(adj or al, key=lambda x: (round(x["qcov"] * x["ident"], 3), -dist(x)))
            n_full = sum(1 for x in al if x["qcov"] >= a.min_qcov)
            d = dist(b)
            if b["chrom"] != chrom:
                rel, d = "other contig", ""
            elif d <= a.adjacent:
                rel = "spans site" if b["rs"] - 1 < pos - a.adjacent and b["re"] > pos + a.adjacent else "adjacent"
            else:
                rel = f"{d} bp away"
            if b["qcov"] < a.min_qcov:
                call = f"partial ({b['qcov']:.2f})"
            elif rel == "adjacent":
                call = "tandem duplication"
            elif rel == "spans site":
                call = "tandem duplication"
            elif rel == "other contig":
                call = "interspersed (other contig)"
            else:
                call = "interspersed duplication"
            o.write(f"{iid}\t{chrom}:{pos}\t{b['chrom']}:{b['rs']}-{b['re']}\t{b['strand']}\t{L}\t{call}\t{rel}\t{d}\t"
                    f"{b['qcov']:.3f}\t{b['ident']:.3f}\t{b['mapq']}\t{b['qs'] + 1}-{b['qe']}\t{len(al)}\t{n_full}\t"
                    f"{band_of(b['chrom'], (b['rs'] + b['re']) // 2)}\n")


if __name__ == "__main__":
    main()
