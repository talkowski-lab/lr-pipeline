#!/usr/bin/env python3
"""Link insertions to segmental-duplication families: the insertion site and the insertion's origin (from any
alignment stage) lying in two copies of one SD family (build_sd_reference.py loci and copies).
For each insertion (chrom-POS-... ID) and each origin interval given with --origin:
  site copies    SD copies (loci) containing POS
  origin copies  SD copies containing >= --min-origin-frac of the origin interval
  link           the origin in one copy and the site in another copy of the same family: "<origin copy> inserted into
                 <site copy>" (copy = GD block name or family letter name; when both carry the same GD block name the
                 family letter names are added); the origin in the site's own copy: "copy within <name>"
  identity       of the SD pair joining the two copies when there is one (direct pair), else empty (same family through
                 other copies)
The best link per insertion (between-copy links first, then direct pairs, then identity) is reported with all others.
--origin is label=path,ID_column,region_column[,column:value] (region chrom:start-end, 1-based; anything after "("
is ignored; with column:value only rows having that value are used, e.g. status:aligned); repeatable.
Outputs <prefix>.tsv: ID, site_copies (up to 5), n_links, best_link, family, best_identity, best_source, best_origin,
all_links.
"""
import argparse
import bisect
import collections
import csv


def read_bed(path):
    with open(path) as f:
        hdr = f.readline().lstrip("#").rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in f]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", required=True, help="build_sd_reference.py <prefix>.loci.bed")
    ap.add_argument("--copies", required=True, help="build_sd_reference.py <prefix>.copies.bed (pair identities)")
    ap.add_argument("--origin", action="append", required=True, help="label=path,ID_column,region_column[,col:value]")
    ap.add_argument("--min-origin-frac", type=float, default=0.5, help="Origin fraction inside a copy [0.5]")
    ap.add_argument("--out-prefix", required=True, help="Output prefix")
    a = ap.parse_args()

    by = collections.defaultdict(list)
    for r in read_bed(a.loci):
        by[r["chrom"]].append((int(r["start"]), int(r["end"]), r["copy"], r["group"], r["name"]))
    idx = {ch: (sorted(v), [x[0] for x in sorted(v)], max(x[1] - x[0] for x in v)) for ch, v in by.items()}
    ident = {}
    for r in read_bed(a.copies):
        k = (r["copy"], r["partner_copy"])
        ident[k] = max(ident.get(k, 0.0), float(r["identity"]))

    def overlapping(ch, s, e):
        if ch not in idx:
            return []
        v, st, ml = idx[ch]
        out = []
        for x in v[bisect.bisect_left(st, s - ml):]:
            if x[0] >= e:
                break
            if x[1] > s:
                out.append(x)
        return out

    origins = collections.defaultdict(list)
    for spec in a.origin:
        lab, rest = spec.split("=", 1)
        path, idc, rc = rest.split(",")[:3]
        flt = rest.split(",")[3:]
        for r in csv.DictReader(open(path), delimiter="\t"):
            if flt and r.get(flt[0].split(":")[0]) != flt[0].split(":")[1]:
                continue
            reg = (r.get(rc) or "").split("(")[0].strip()
            if ":" not in reg:
                continue
            ch, rng = reg.split(":")
            s, e = (int(x) for x in rng.split("-"))
            origins[r[idc]].append((lab, ch, s - 1, e))
    rows = []
    for iid, ol in origins.items():
        ch, pos = iid.split("-")[0], int(iid.split("-")[1])
        site = [x for x in overlapping(ch, pos - 1, pos + 1) if x[0] < pos <= x[1]]
        links = {}
        for lab, oc, os_, oe in ol:
            for x in overlapping(oc, os_, oe):
                if (min(x[1], oe) - max(x[0], os_)) / max(1, oe - os_) < a.min_origin_frac:
                    continue
                for y in site:
                    if y[3] != x[3]:
                        continue
                    between = y[2] != x[2]
                    idn = ident.get((x[2], y[2]), 0.0)
                    if not between:
                        txt = f"copy within {y[4]}"
                    elif x[4] == y[4]:
                        txt = f"{x[4]} ({x[2]}) inserted into {y[4]} ({y[2]})"
                    else:
                        txt = f"{x[4]} inserted into {y[4]}"
                    key = (x[2], y[2])
                    if key not in links:
                        links[key] = (between, idn > 0, idn, txt, x[3], lab, f"{oc}:{os_ + 1}-{oe}")
        names = sorted({y[4] for y in site})
        ls = sorted(links.values(), key=lambda t: (not t[0], not t[1], -t[2]))
        b = ls[0] if ls else None
        rows.append((iid, ",".join(names[:5]) + (",..." if len(names) > 5 else ""), len(ls), b[3] if b else "",
                     b[4] if b else "", f"{b[2]:.4f}" if b and b[1] else "", b[5] if b else "", b[6] if b else "",
                     " | ".join(f"{t[3]} ({t[4]}, {f'{100 * t[2]:.1f}%' if t[1] else 'same family'}, {t[5]})" for t in ls[:10])))
    with open(a.out_prefix + ".tsv", "w") as o:
        o.write("ID\tsite_copies\tn_links\tbest_link\tfamily\tbest_identity\tbest_source\tbest_origin\tall_links\n")
        for r in sorted(rows, key=lambda r: (int(r[0].split("-")[1]), r[0])):
            o.write("\t".join(str(v) for v in r) + "\n")
    n = sum(1 for r in rows if r[3])
    nb = sum(1 for r in rows if r[3] and "inserted into" in r[3])
    print(f"{len(rows)} insertions with origins; {n} linked to an SD family ({nb} between copies, {n - nb} within one copy)")


if __name__ == "__main__":
    main()
