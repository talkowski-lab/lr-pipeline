#!/usr/bin/env python3
"""Genomic context of insertion sites (POS of chrom-POS-... IDs), for labelling insertions that sit in repetitive
regions. Per ID, comma-joined in a "context" column (empty when none applies):
  SD <block>                  inside a named segmental-duplication block (--sd-blocks, e.g. LCR22A-H)
  SD <identity>% (partner)    inside a UCSC genomicSuperDups copy (highest-identity pair; partner = other copy)
  subtelomeric                within --subtel-bp of a chromosome end
  pericentromeric             within --pericen-bp of the centromere (acen bands)
  acrocentric short arm       on the p arm of chr13, chr14, chr15, chr21, chr22
  TR <motifs>                 within --pad bp of tandem-repeat catalog loci (--tr-catalog, MOTIFS= in column 4; the
                              three shortest motifs)
  <class> <family>            within --pad bp of a RepeatMasker-style element of class Simple_repeat, Satellite or
                              Low_complexity (--repeats)
"""
import argparse
import bisect
import collections
import gzip

ACRO = {"chr13", "chr14", "chr15", "chr21", "chr22"}


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def index(rows):
    d = collections.defaultdict(list)
    for ch, s, e, v in rows:
        d[ch].append((s, e, v))
    out = {}
    for ch, v in d.items():
        v.sort()
        out[ch] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
    return out


def hits(idx, ch, p, pad):
    if ch not in idx:
        return []
    v, st, ml = idx[ch]
    i = bisect.bisect_left(st, p - pad - ml)
    res = []
    for s, e, val in v[i:]:
        if s > p + pad:
            break
        if e >= p - pad:
            res.append(val)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", required=True, help="File of insertion IDs (chrom-POS-...)")
    ap.add_argument("--cytoband", required=True, help="UCSC cytoBand table")
    ap.add_argument("--sds", required=True, help="UCSC genomicSuperDups.txt.gz")
    ap.add_argument("--sd-blocks", default=None, help="Optional BED of named SD blocks (chrom start end name)")
    ap.add_argument("--tr-catalog", default=None, help="Optional TR catalog BED(.gz) with MOTIFS= in column 4")
    ap.add_argument("--repeats", default=None, help="Optional RepeatMasker-style BED.gz (ref beg end name class family)")
    ap.add_argument("--subtel-bp", type=int, default=500000, help="Distance to a chromosome end for subtelomeric [500000]")
    ap.add_argument("--pericen-bp", type=int, default=3000000, help="Distance to the centromere for pericentromeric [3000000]")
    ap.add_argument("--pad", type=int, default=50, help="Distance to a TR / repeat element [50]")
    ap.add_argument("--out", required=True, help="Output TSV (ID, context)")
    a = ap.parse_args()

    ids = [x.strip() for x in open(a.ids) if x.strip()]
    chroms = {i.split("-")[0] for i in ids}
    clen, cen = {}, collections.defaultdict(list)
    for line in open(a.cytoband):
        ch, s, e, nm, stn = line.rstrip("\n").split("\t")[:5]
        clen[ch] = max(clen.get(ch, 0), int(e))
        if stn == "acen":
            cen[ch].append((int(s), int(e)))
    sd = []
    with gzip.open(a.sds, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[1] in chroms:
                sd.append((c[1], int(c[2]) + 1, int(c[3]), (float(c[26]), f"{c[7]}:{int(c[8]) + 1}-{c[9]}")))
    sdi = index(sd)
    blk = {}
    if a.sd_blocks:
        rows = [line.rstrip("\n").split("\t") for line in open(a.sd_blocks) if line.strip() and not line.startswith("#")]
        blk = index([(c[0], int(c[1]) + 1, int(c[2]), c[3]) for c in rows])
    tri = {}
    if a.tr_catalog:
        rows = []
        with opener(a.tr_catalog) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if c[0] in chroms:
                    m = next((x[7:] for x in c[3].split(";") if x.startswith("MOTIFS=")), "")
                    rows.append((c[0], int(c[1]) + 1, int(c[2]), m))
        tri = index(rows)
    rpi = {}
    if a.repeats:
        rows = []
        with opener(a.repeats) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if c[0] in chroms and c[4] in ("Simple_repeat", "Satellite", "Low_complexity"):
                    rows.append((c[0], int(c[1]) + 1, int(c[2]), c[4] if c[4] == c[5] else f"{c[4]} {c[5]}"))
        rpi = index(rows)
    with open(a.out, "w") as o:
        o.write("ID\tcontext\n")
        for iid in ids:
            ch, p = iid.split("-")[0], int(iid.split("-")[1])
            ctx = [f"SD {b}" for b in hits(blk, ch, p, 0)]
            s = hits(sdi, ch, p, 0)
            if s:
                idn, partner = max(s)
                ctx.append(f"SD {100 * idn:.1f}% ({partner})")
            if ch in clen and min(p, clen[ch] - p) <= a.subtel_bp:
                ctx.append("subtelomeric")
            if cen.get(ch):
                cs, ce = min(x[0] for x in cen[ch]), max(x[1] for x in cen[ch])
                if (0 if cs <= p <= ce else min(abs(p - cs), abs(p - ce))) <= a.pericen_bp:
                    ctx.append("pericentromeric")
                if ch in ACRO and p < cs:
                    ctx.append("acrocentric short arm")
            motifs = sorted({m for x in hits(tri, ch, p, a.pad) for m in x.split(",") if m}, key=lambda m: (len(m), m))
            if motifs:
                ctx.append("TR " + "/".join(motifs[:3]) + ("/..." if len(motifs) > 3 else ""))
            ctx += sorted(set(hits(rpi, ch, p, a.pad)))[:2]
            o.write(f"{iid}\t{', '.join(ctx)}\n")
    print(f"{len(ids)} sites")


if __name__ == "__main__":
    main()
