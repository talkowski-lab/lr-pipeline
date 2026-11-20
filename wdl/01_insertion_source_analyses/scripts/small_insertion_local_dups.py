#!/usr/bin/env python3
"""Local duplication test for short insertions (e.g. 20-49 bp) at their insertion site.
For an insertion of length L after POS:
  tandem      the inserted sequence equals a tandem copy of the reference at the site: for some split k (0..L) it is
              ref[POS+1 .. POS+L-k] followed by ref[POS-k+1 .. POS] (k = 0: copy of the L bases right of the site;
              k = L: copy of the L bases left of it; in between: the insertion point lies inside the duplicated unit);
              identity = best fraction of equal bases over k (ungapped); the same on the reverse complement = inverted
  nearby      best gapped (Smith-Waterman, +2 / -3 / -5) alignment of the insertion, both strands, to the reference
              within --window bp of the site; coverage = aligned insertion bases / L, identity over aligned columns,
              distance from the site
Calls (identity >= --min-ident, coverage >= --min-cov for nearby): tandem duplication > inverted tandem duplication >
interspersed duplication (nearby copy more than --adjacent bp from the site) > copy at the site with indels (gapped
copy touching the site) > no local copy.
Input: a remap-style TSV (CHROM POS REF ALT ID ...; inserted sequence = ALT minus the 1-base anchor) and optional
--ids. Output: one row per insertion.
"""
import argparse

import pysam


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def tandem_ident(q, right, left):
    L = len(q)
    best = (0.0, 0)
    for k in range(L + 1):
        cand = right[:L - k] + left[L - k:] if k else right[:L]
        m = sum(1 for x, y in zip(q, cand) if x == y) / L
        if m > best[0]:
            best = (m, k)
    return best


def sw(q, t, match=2, mismatch=-3, gap=-5):
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
    i, j, cols, idn = bi, bj, 0, 0
    while i > 0 and j > 0:
        k = trace[i - 1][j]
        if k == 0:
            break
        cols += 1
        if k == 1:
            idn += q[i - 1] == t[j - 1]
            i, j = i - 1, j - 1
        elif k == 2:
            i -= 1
        else:
            j -= 1
    return i, bi, j, bj, idn / max(1, cols)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--remap", required=True, help="Remap-style TSV (CHROM POS REF ALT ID ...)")
    ap.add_argument("--ids", default=None, help="Optional file of IDs to test")
    ap.add_argument("--ref-fa", required=True, help="Reference FASTA with .fai")
    ap.add_argument("--window", type=int, default=1000, help="Reference searched on each side for nearby copies [1000]")
    ap.add_argument("--min-ident", type=float, default=0.9, help="Identity for a copy [0.9]")
    ap.add_argument("--min-cov", type=float, default=0.9, help="Insertion fraction in a nearby copy [0.9]")
    ap.add_argument("--adjacent", type=int, default=5, help="A copy within this distance of the site is at the site [5]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    keep = {x.strip() for x in open(a.ids) if x.strip()} if a.ids else None
    ref = pysam.FastaFile(a.ref_fa)
    with open(a.out, "w") as o:
        o.write("ID\tlen\tinserted_seq\ttandem_identity\ttandem_split\tinv_tandem_identity\tnearby_copy\tnearby_strand\t"
                "nearby_cov\tnearby_identity\tnearby_dist\tcall\n")
        for line in open(a.remap):
            c = line.rstrip("\n").split("\t")
            if keep is not None and c[4] not in keep:
                continue
            ch, pos, q = c[0], int(c[1]), c[3].upper()[1:]
            L = len(q)
            right = ref.fetch(ch, pos, pos + L).upper()
            left = ref.fetch(ch, max(0, pos - L), pos).upper()
            ti, tk = tandem_ident(q, right, left)
            ii, _ = tandem_ident(rc(q), right, left)
            w0 = max(0, pos - a.window)
            win = ref.fetch(ch, w0, pos + a.window).upper()
            best = None
            for strand, qq in (("+", q), ("-", rc(q))):
                qs, qe, ts, te, idn = sw(qq, win)
                cov = (qe - qs) / L
                key = cov * idn
                if best is None or key > best[0]:
                    rs, re_ = w0 + ts + 1, w0 + te
                    dist = 0 if rs - 1 <= pos <= re_ else min(abs(rs - pos), abs(re_ - pos))
                    best = (key, strand, cov, idn, f"{ch}:{rs}-{re_}", dist)
            _, st, cov, idn, reg, dist = best
            if ti >= a.min_ident:
                call = "tandem duplication" + ("" if tk in (0, L) else ", insertion point inside the unit")
            elif ii >= a.min_ident:
                call = "inverted tandem duplication"
            elif cov >= a.min_cov and idn >= a.min_ident and dist > a.adjacent:
                call = "interspersed duplication (nearby)" if st == "+" else "inverted interspersed duplication (nearby)"
            elif cov >= a.min_cov and idn >= a.min_ident:
                call = "copy at the site with indels"
            else:
                call = "no local copy"
            o.write(f"{c[4]}\t{L}\t{q}\t{ti:.3f}\t{tk}\t{ii:.3f}\t{reg}\t{st}\t{cov:.3f}\t{idn:.3f}\t{dist}\t{call}\n")


if __name__ == "__main__":
    main()
