#!/usr/bin/env python3
"""Check short insertions against T2T (hs1) from a BLAT PSL of three queries per insertion: <ID>|ins (inserted sequence),
<ID>|L and <ID>|R (hg38 flanks left / right of the site, e.g. 200 bp each).
The flanks place the site in T2T: best L hit and best R hit (coverage and identity >= --flank-min) on the same contig and
strand, within --max-gap of each other (the insertion's own chromosome first, then the best flank
alignments); the T2T site is the span between the end of L and the start of R (in L->R order).
Insertion hits count when coverage (matched + mismatched bases / query length) and identity (matches / aligned) are both
>= --min-cov / --min-ident.
Calls:
  present in T2T at the site           an insertion hit lies inside the T2T flank span and the span is >= --min-span x the
                                       insertion length (default 0.5, flank ends may be trimmed): T2T carries the
                                       inserted sequence where hg38 has none
  copy near the T2T site               an insertion hit within --window bp of the T2T site (a local copy in T2T)
  copy elsewhere in T2T                insertion hits only away from the site (or the site is not placed); nearest listed
  no T2T copy                          no insertion hit
Output: one row per insertion with the T2T site, flank span, number of insertion hits, best hit and call.
"""
import argparse
from collections import defaultdict


def parse_psl(path):
    hits, qsize = defaultdict(list), {}
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if len(c) < 21 or not c[0].isdigit():
            continue
        match, mis, rep = int(c[0]), int(c[1]), int(c[2])
        aligned = match + mis + rep
        qsize[c[9]] = int(c[10])
        hits[c[9]].append({"strand": c[8], "chrom": c[13], "ts": int(c[15]), "te": int(c[16]),
                           "cov": aligned / int(c[10]), "ident": (match + rep) / max(1, aligned)})
    return hits, qsize


def place_site(L, R, max_gap, chrom):
    """T2T site from flank hits: (chrom, start, end) of the span between L and R in L->R order. Pairs on the insertion's
    own chromosome first, then the best flank alignments, then the smallest gap."""
    best, bkey = None, None
    for lh in L:
        for rh in R:
            if lh["chrom"] != rh["chrom"] or lh["strand"] != rh["strand"]:
                continue
            s, e = (lh["te"], rh["ts"]) if lh["strand"] == "+" else (rh["te"], lh["ts"])
            if not -50 <= e - s <= max_gap:
                continue
            key = (lh["chrom"] != chrom, -(lh["cov"] * lh["ident"] + rh["cov"] * rh["ident"]), abs(e - s))
            if bkey is None or key < bkey:
                best, bkey = (lh["chrom"], s, e), key
    return best


def dist_to_site(h, site):
    if site is None or h["chrom"] != site[0]:
        return None
    lo, hi = min(site[1], site[2]), max(site[1], site[2])
    if h["ts"] <= hi and h["te"] >= lo:
        return 0
    return min(abs(h["ts"] - hi), abs(h["te"] - lo))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psl", required=True, help="BLAT PSL of the <ID>|ins / |L / |R queries against T2T")
    ap.add_argument("--ids", required=True, help="Insertion IDs to report (one per line)")
    ap.add_argument("--min-cov", type=float, default=0.9, help="Insertion hit coverage [0.9]")
    ap.add_argument("--min-ident", type=float, default=0.9, help="Insertion hit identity [0.9]")
    ap.add_argument("--flank-min", type=float, default=0.95, help="Flank hit coverage and identity [0.95]")
    ap.add_argument("--max-gap", type=int, default=10000, help="Largest T2T distance between the L and R flank hits [10000]")
    ap.add_argument("--min-span", type=float, default=0.5, help="T2T span / insertion length, present at site [0.5]")
    ap.add_argument("--window", type=int, default=1000, help="Distance from the T2T site for a nearby copy [1000]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    hits, qsize = parse_psl(a.psl)
    with open(a.out, "w") as o:
        o.write("ID\tins_len\tt2t_site\tt2t_span\tn_ins_hits\tbest_hit\tbest_cov\tbest_ident\tdist_to_site\tcall\n")
        for i in (x.strip() for x in open(a.ids) if x.strip()):
            n = qsize.get(f"{i}|ins")
            ins = [h for h in hits.get(f"{i}|ins", []) if h["cov"] >= a.min_cov and h["ident"] >= a.min_ident]
            L = [h for h in hits.get(f"{i}|L", []) if h["cov"] >= a.flank_min and h["ident"] >= a.flank_min]
            R = [h for h in hits.get(f"{i}|R", []) if h["cov"] >= a.flank_min and h["ident"] >= a.flank_min]
            site = place_site(L, R, a.max_gap, i.split("-")[0])
            span = site[2] - site[1] if site else None
            scored = sorted(((dist_to_site(h, site), h) for h in ins),
                            key=lambda t: (t[0] is None, t[0] or 0, -t[1]["cov"] * t[1]["ident"]))
            d, best = scored[0] if scored else (None, None)
            if best is None:
                call = "no T2T copy"
            elif d == 0 and n and span is not None and span >= a.min_span * n:
                call = "present in T2T at the site"
            elif d is not None and d <= a.window:
                call = "copy near the T2T site"
            else:
                call = "copy elsewhere in T2T"
            o.write("\t".join([i, str(n or ""), f"{site[0]}:{site[1]}-{site[2]}" if site else "not placed",
                               "" if span is None else str(span), str(len(ins)),
                               f"{best['chrom']}:{best['ts'] + 1}-{best['te']}({best['strand']})" if best else "",
                               f"{best['cov']:.3f}" if best else "", f"{best['ident']:.3f}" if best else "",
                               "" if d is None else str(d), call]) + "\n")


if __name__ == "__main__":
    main()
