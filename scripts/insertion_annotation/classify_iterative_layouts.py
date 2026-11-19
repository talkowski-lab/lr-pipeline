#!/usr/bin/env python3
"""Classify insertions from AlignQueriesIterative segments on two references (site reference, e.g. hg38, and a second
one, e.g. T2T-CHM13), for insertions whose aligned segments cover >= --min-cov of the inserted sequence on either.
Per reference: coverage, segment-length-weighted identity, loci (segments on one contig within --merge-bp; only
segments >= --min-locus-len bp at >= --min-ident identity define loci). Classes, first match wins:
  rotated tandem duplication   site reference, all + segments on the site contig, 5' group right of the site and 3'
                               group left of it, each allowed to overlap the site by --adjacent bp (gaps reported)
  tandem duplication (pieces)  site reference, all + segments on the site contig on one side, the nearest within
                               --adjacent bp of the site
  local repeat expansion       site reference, one locus touching or within --adjacent bp of the site, segments
                               overlapping each other on the reference (the same local sequence copied several times)
  chimeric insertion           >= 2 loci on either reference
  present in <second> only     covered on the second reference but not on the site reference
  low-identity match           covered, but weighted identity < --min-ident
  other single-locus           anything else (one locus away from the site)
"""
import argparse
import collections
import csv


def load(path):
    seg = collections.defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        seg[r["ID"]].append(r)
    return seg


def describe(ss, a):
    L = max(int(r["q_end"]) for r in ss)
    al = sorted((r for r in ss if r["status"] == "aligned"), key=lambda r: int(r["q_start"]))
    bp = [int(r["q_end"]) - int(r["q_start"]) + 1 for r in al]
    iv = []
    for r in al:
        c, rng = r["ref"].split(":")
        s, e = (int(x) for x in rng.split("-"))
        iv.append((c, s, e, r["strand"], float(r["identity"])))
    loci = []
    for c, s, e, _, idn in sorted((x for x, b in zip(iv, bp) if b >= a.min_locus_len and x[4] >= a.min_ident),
                                  key=lambda x: (x[0], x[1])):
        if loci and loci[-1][0] == c and s <= loci[-1][2] + a.merge_bp:
            loci[-1][2] = max(loci[-1][2], e)
        else:
            loci.append([c, s, e])
    return dict(L=L, bp=bp, cov=sum(bp) / L, wid=sum(x[4] * b for x, b in zip(iv, bp)) / max(1, sum(bp)), iv=iv, loci=loci,
                lay=" | ".join(f"q{r['q_start']}-{r['q_end']}>{r['ref']}({r['strand']},{float(r['identity']):.2f})" for r in al))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site-segments", required=True, help="Segments TSV on the reference of the insertion sites")
    ap.add_argument("--other-segments", required=True, help="Segments TSV on the second reference")
    ap.add_argument("--other-label", default="T2T", help="Name of the second reference [T2T]")
    ap.add_argument("--min-cov", type=float, default=0.9, help="Min aligned fraction of the insertion [0.9]")
    ap.add_argument("--min-ident", type=float, default=0.9, help="Identity for loci and for a good match [0.9]")
    ap.add_argument("--min-locus-len", type=int, default=50, help="Min segment length defining a locus [50]")
    ap.add_argument("--merge-bp", type=int, default=10000, help="Segments within this distance form one locus [10000]")
    ap.add_argument("--adjacent", type=int, default=5, help="Distance from the site counted as tandem [5]")
    ap.add_argument("--t2t-calls-out", default=None,
                    help="Optional name_t2t_matches.py input for the second-reference-only rows (largest segment = best hit)")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    hs, ts = load(a.site_segments), load(a.other_segments)
    out, t2t_rows = [], []
    for iid in sorted(set(hs) | set(ts)):
        h = describe(hs[iid], a) if iid in hs else None
        t = describe(ts[iid], a) if iid in ts else None
        hc, tc = (h["cov"] if h else 0), (t["cov"] if t else 0)
        if max(hc, tc) < a.min_cov:
            continue
        ch, pos = iid.split("-")[0], int(iid.split("-")[1])
        cls, ref, x, note = "", "", None, ""
        if hc >= a.min_cov:
            ref, x = "site", h
            iv = h["iv"]
            same = all(c == ch and st == "+" for c, _, _, st, _ in iv)
            sides = ["R" if s >= pos + 1 - a.adjacent else "L" if e <= pos + a.adjacent else "spans" for _, s, e, _, _ in iv]
            k = next((i for i in range(1, len(sides)) if sides[:i] == ["R"] * i and sides[i:] == ["L"] * (len(sides) - i)), None)
            dist = [0 if (c == ch and s - 1 <= pos <= e) else (min(abs(s - pos), abs(e - pos)) if c == ch else None)
                    for c, s, e, _, _ in iv]
            near = [d for d in dist if d is not None]
            if same and k is not None:
                gr = min(s for _, s, _, _, _ in iv[:k]) - (pos + 1)
                gl = pos - max(e for _, _, e, _, _ in iv[k:])
                cls, note = "rotated tandem duplication", f"gap right {gr}, gap left {gl}"
            elif same and set(sides) in ({"R"}, {"L"}) and min(near) <= a.adjacent:
                cls, note = "tandem duplication (pieces)", f"nearest piece {min(near)} bp from the site"
            elif len(h["loci"]) <= 1 and near and min(near) <= a.adjacent and len(iv) >= 2 and h["wid"] >= a.min_ident - 0.05:
                cls, note = "local repeat expansion", f"{len(iv)} overlapping copies of local sequence"
        if not cls and (len(h["loci"]) if h else 0) >= 2 and hc >= a.min_cov:
            cls, ref, x = "chimeric insertion", "site", h
        if not cls and (len(t["loci"]) if t else 0) >= 2 and tc >= a.min_cov and hc < a.min_cov:
            cls, ref, x = "chimeric insertion", a.other_label, t
        if not cls and hc < a.min_cov:
            ref, x = a.other_label, t
            cls = f"present in {a.other_label} only" if t["wid"] >= a.min_ident else "low-identity match"
        if not cls:
            ref, x = "site", h
            cls = "low-identity match" if h["wid"] < a.min_ident else "other single-locus"
        if cls == f"present in {a.other_label} only":
            big = max(zip(x["iv"], x["bp"]), key=lambda t: t[1])
            c, s0, e0, st, idn = big[0]
            t2t_rows.append((iid, x["L"], "present in T2T-CHM13 (iterative alignment)", f"{100 * big[1] / x['L']:.1f}",
                             f"{c}:{s0}-{e0}({st},{100 * idn:.1f}%)", f"{tc:.3f}", "", f"{hc:.3f}"))
        out.append((iid, x["L"], f"{hc:.3f}", f"{tc:.3f}", ref, cls, note, len(x["iv"]), len(x["loci"]), f"{x['wid']:.3f}",
                    ",".join(f"{c}:{s}-{e}" for c, s, e in x["loci"]), x["lay"]))
    with open(a.out, "w") as o:
        o.write("ID\tlen\tsite_ref_cov\tother_ref_cov\tref_used\tclass\tnote\tn_segments\tn_loci\tweighted_identity\t"
                "loci\tlayout\n")
        for r in sorted(out, key=lambda r: (r[5], r[0])):
            o.write("\t".join(str(v) for v in r) + "\n")
    if a.t2t_calls_out:
        with open(a.t2t_calls_out, "w") as o:
            o.write("ID\tlen\tt2t_call\tT2T_best_perc\tT2T_best_hit\tT2T_cov\thg38_best_perc\thg38_cov\n")
            o.writelines("\t".join(str(v) for v in r) + "\n" for r in t2t_rows)
    for k, v in collections.Counter(r[5] for r in out).most_common():
        print(v, k)


if __name__ == "__main__":
    main()
