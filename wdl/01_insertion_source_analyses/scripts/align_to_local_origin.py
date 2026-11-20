#!/usr/bin/env python3
"""Gapped local alignment of inserted sequences to their local origin (locate_local_source.py), to call interspersed
duplications.
For each ID, the inserted sequence (reverse complement when local_strand is "-") is aligned with Smith-Waterman
(match +2, mismatch -3, gap -5 per base) to the origin interval padded by --pad bp on each side. Reported: the aligned
reference interval, the aligned query span, identity over aligned columns, and aligned_frac = aligned query span /
insertion length.
Call: "interspersed duplication" when aligned_frac >= --min-frac and the aligned reference interval is more than
--adjacent bp from the insertion site (POS) and does not contain it; otherwise "not called" (reason given).
Insertions longer than --max-len are skipped.
"""
import argparse
import csv

import pysam


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def smith_waterman(q, t, match=2, mismatch=-3, gap=-5):
    n, m = len(q), len(t)
    prev = [0] * (m + 1)
    trace = []
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        cur = [0] * (m + 1)
        tr = bytearray(m + 1)
        qi = q[i - 1]
        for j in range(1, m + 1):
            d = prev[j - 1] + (match if qi == t[j - 1] else mismatch)
            u = prev[j] + gap
            left = cur[j - 1] + gap
            v, k = d, 1
            if u > v:
                v, k = u, 2
            if left > v:
                v, k = left, 3
            if v <= 0:
                v, k = 0, 0
            cur[j] = v
            tr[j] = k
            if v > best:
                best, bi, bj = v, i, j
        trace.append(tr)
        prev = cur
    i, j, cols, ident = bi, bj, 0, 0
    while i > 0 and j > 0:
        k = trace[i - 1][j]
        if k == 0:
            break
        cols += 1
        if k == 1:
            ident += q[i - 1] == t[j - 1]
            i, j = i - 1, j - 1
        elif k == 2:
            i -= 1
        else:
            j -= 1
    return best, i, bi, j, bj, ident / max(1, cols)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--local-source", required=True, help="locate_local_source.py output (ID, local_source, local_strand, ...)")
    ap.add_argument("--ids", required=True, help="File of insertion IDs")
    ap.add_argument("--remap", required=True, help="Truvari remap TSV, no header (CHROM POS REF ALT ID ...)")
    ap.add_argument("--ref-fa", required=True, help="Reference FASTA with .fai")
    ap.add_argument("--pad", type=int, default=200, help="Padding around the origin interval, bp [200]")
    ap.add_argument("--min-frac", type=float, default=0.9, help="Min aligned fraction of the insertion [0.9]")
    ap.add_argument("--adjacent", type=int, default=5, help="Origin within this distance of the site = tandem, not called [5]")
    ap.add_argument("--max-len", type=int, default=5000, help="Skip insertions longer than this [5000]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    ids = [x.strip() for x in open(a.ids) if x.strip()]
    src = {r["ID"]: r for r in csv.DictReader(open(a.local_source), delimiter="\t")}
    recs = {}
    with open(a.remap) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[4] in src:
                recs[c[4]] = (c[0], int(c[1]), c[3].upper()[1:])
    ref = pysam.FastaFile(a.ref_fa)
    with open(a.out, "w") as o:
        o.write("ID\tsite\tlen\torigin_aligned\tstrand\tq_start\tq_end\taligned_frac\tidentity\tdist_to_site\tcall\t"
                "inserted_seq\n")
        for iid in ids:
            ch, pos, ins = recs[iid]
            s = src[iid]
            if not s["local_source"]:
                o.write(f"{iid}\t{ch}:{pos}\t{len(ins)}\t\t\t\t\t\t\t\tno local source\t{ins}\n")
                continue
            if len(ins) > a.max_len:
                o.write(f"{iid}\t{ch}:{pos}\t{len(ins)}\t\t\t\t\t\t\t\tskipped (longer than {a.max_len} bp)\t{ins}\n")
                continue
            sch, rng = s["local_source"].split(":")
            ss, se = (int(x) for x in rng.split("-"))
            t0 = max(0, ss - 1 - a.pad)
            t = ref.fetch(sch, t0, se + a.pad).upper()
            q = rc(ins) if s["local_strand"] == "-" else ins
            _, qs, qe, ts, te, ident = smith_waterman(q, t)
            if s["local_strand"] == "-":
                qs, qe = len(q) - qe, len(q) - qs
            rs, re_ = t0 + ts + 1, t0 + te
            frac = (qe - qs) / len(ins)
            inside = sch == ch and rs - 1 <= pos <= re_
            dist = 0 if inside else min(abs(rs - pos), abs(re_ - pos))
            if sch == ch and dist <= a.adjacent:
                call = f"not called (tandem: origin within {a.adjacent} bp of the site)"
            elif frac < a.min_frac:
                call = f"not called (aligned < {a.min_frac:g})"
            else:
                call = "interspersed duplication"
            o.write(f"{iid}\t{ch}:{pos}\t{len(ins)}\t{sch}:{rs}-{re_}\t{s['local_strand']}\t{qs + 1}\t{qe}\t{frac:.3f}\t"
                    f"{ident:.3f}\t{dist}\t{call}\t{ins}\n")
            print(iid, f"{sch}:{rs}-{re_}", f"frac {frac:.3f}", f"ident {ident:.3f}", f"dist {dist}", call)


if __name__ == "__main__":
    main()
