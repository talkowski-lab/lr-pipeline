#!/usr/bin/env python3
"""Last pipeline step: decide between tandem duplication, complex SV and chimeric insertion from the local realignment
(run_site_copy_realign.sh: iterative alignment of each insertion against the site's SD copy or the reference around
the site) and the genome-wide iterative alignment (best placement of every piece).
Per insertion:
  local      pieces aligned to the local target (coordinates converted to the chromosome), their coverage of the
             insertion and segment-length-weighted identity
  genome     coverage and weighted identity of the genome-wide pieces
  identity_drop = genome weighted identity - local weighted identity (over their own aligned bases)
Direct tandem test first (run_site_copy_realign.sh right.paf / left.paf): the insertion's prefix aligned on the +
strand starting at POS+1 and its suffix aligned on the + strand ending at POS; when together they cover >= --min-cov
with at most max(--site-tol, --max-gap-frac x length) bases between them and identity within --max-identity-drop of the
genome-wide placement, the insertion is a tandem duplication of unit (suffix start .. prefix end).
Then the piece-layout test (a tandem duplication of unit U inserted after POS reads the unit from POS+1, wrapping
from the unit end to its start as often as needed, and ends at POS): local pieces in insertion order, all + strand; the
first starts at
POS+1 and the last ends at POS (or at the unit edges when POS is at an edge), each within --site-tol bp; between
consecutive pieces the reference continues (next start = previous end + 1, with up to --site-tol bp shared bases or a
gap up to --max-gap-frac x unit) or wraps from the unit end to the unit start (same tolerances). U = span of the local
pieces; copies = insertion length / unit length.
Classes (first match wins):
  tandem duplication   tandem test passes, local coverage >= --min-cov, identity_drop <= --max-identity-drop
  complex SV           local coverage >= --min-cov and identity_drop <= --max-identity-drop, tandem test fails
                       (pieces out of order or inverted, all local)
  chimeric insertion   genome coverage >= --min-cov (pieces from other places); origins listed
  unresolved           neither covers --min-cov
"""
import argparse
import collections
import csv
import os


def load(path):
    seg = collections.defaultdict(list)
    if os.path.exists(path):
        for r in csv.DictReader(open(path), delimiter="\t"):
            seg[r["ID"]].append(r)
    return seg


def summary(rows, L, conv=None):
    al = []
    for r in sorted(rows, key=lambda r: int(r["q_start"])):
        if r["status"] != "aligned":
            continue
        c, rng = r["ref"].split(":", 1)
        if conv:  # target-relative "chrom:start-end" reference name
            tc, trng = r["ref"].rsplit(":", 1)
            s, e = (int(x) for x in trng.split("-"))
            base = conv[tc]
            c, s, e = base[0], base[1] + s - 1, base[1] + e - 1
        else:
            s, e = (int(x) for x in rng.split("-"))
        al.append(dict(qs=int(r["q_start"]), qe=int(r["q_end"]), chrom=c, s=s, e=e, strand=r["strand"],
                       ident=float(r["identity"]), mapq=r["mapq"]))
    bp = sum(x["qe"] - x["qs"] + 1 for x in al)
    wid = sum(x["ident"] * (x["qe"] - x["qs"] + 1) for x in al) / max(1, bp)
    return al, bp / max(1, L), wid


def tandem_test(al, pos, tol, gap_frac):
    if not al or any(x["strand"] != "+" for x in al) or len({x["chrom"] for x in al}) > 1:
        return False, "", "inverted or other contig" if al else "no local pieces"
    us, ue = min(x["s"] for x in al), max(x["e"] for x in al)
    U = ue - us + 1
    gap = max(tol, int(gap_frac * U))
    f, last = al[0], al[-1]
    start_ok = abs(f["s"] - (pos + 1)) <= tol or (pos >= ue - tol and abs(f["s"] - us) <= tol)
    end_ok = abs(last["e"] - pos) <= tol or (pos <= us + tol and abs(last["e"] - ue) <= tol)
    if not (start_ok and end_ok):
        return False, f"chr:{us}-{ue}", (f"first piece starts {f['s'] - pos - 1:+d} bp from POS+1, last ends "
                                         f"{last['e'] - pos:+d} bp from POS")
    for p, n in zip(al, al[1:]):
        d = n["s"] - (p["e"] + 1)
        wrap = abs(p["e"] - ue) <= gap and abs(n["s"] - us) <= gap
        if not (-tol <= d <= gap or wrap):
            return False, f"chr:{us}-{ue}", f"pieces not contiguous: {p['e']} -> {n['s']}"
    return True, f"{al[0]['chrom']}:{us}-{ue}", f"unit {U} bp"


def flank_test(d, iid, L, pos, tol):
    """Direct tandem test: the insertion's prefix aligned (+) from the first base right of POS (right.paf) and its
    suffix aligned (+) up to POS (left.paf). Returns (prefix_bp, prefix_identity, unit_end, suffix_start, suffix_identity,
    unit_start) with 0 / None when absent."""
    def read(side):
        out = []
        p = os.path.join(d, f"{iid}.{side}.paf")
        if os.path.exists(p):
            for line in open(p):
                c = line.rstrip("\n").split("\t")
                if c[4] == "+":
                    tname, tlen = c[5], int(c[6])
                    t0 = int(tname.rsplit(":", 1)[1].split("-")[0])
                    out.append((int(c[2]), int(c[3]), int(c[7]), int(c[8]), tlen, int(c[9]) / max(1, int(c[10])), t0))
        return out
    pre = [x for x in read("right") if x[0] <= tol and x[2] <= tol]
    suf = [x for x in read("left") if x[1] >= L - tol and x[3] >= x[4] - tol]
    a_ = max(pre, key=lambda x: x[1]) if pre else None
    b_ = min(suf, key=lambda x: x[0]) if suf else None
    return ((a_[1], a_[5], a_[6] + a_[3] - 1) if a_ else (0, None, None),
            (b_[0], b_[5], b_[6] + b_[2]) if b_ else (L, None, None))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True,
                    help="prepare_site_copy_realign.py manifest.tsv (its directory holds the local alignments)")
    ap.add_argument("--genome-segments", required=True, help="Genome-wide AlignQueriesIterative segments TSV (site reference)")
    ap.add_argument("--loci", default=None, help="Optional build_sd_reference.py loci.bed to name piece origins")
    ap.add_argument("--min-cov", type=float, default=0.9, help="Coverage of the insertion for a call [0.9]")
    ap.add_argument("--max-identity-drop", type=float, default=0.01,
                    help="Allowed weighted identity drop, local vs genome [0.01]")
    ap.add_argument("--site-tol", type=int, default=50, help="Shared bases / offset allowed at the site and junctions, bp [50]")
    ap.add_argument("--max-gap-frac", type=float, default=0.05, help="Gap allowed between consecutive pieces, x unit [0.05]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    d = os.path.dirname(os.path.abspath(a.manifest))
    gseg = load(a.genome_segments)
    names = collections.defaultdict(list)
    if a.loci:
        with open(a.loci) as f:
            f.readline()
            for line in f:
                c = line.rstrip("\n").split("\t")
                names[c[0]].append((int(c[1]), int(c[2]), c[5]))

    def sd_name(c, s, e):
        hit = [n for a0, b0, n in names.get(c, []) if min(b0, e) - max(a0, s - 1) >= 0.5 * (e - s + 1)]
        return hit[0] if hit else ""

    rows = []
    for m in csv.DictReader(open(a.manifest), delimiter="\t"):
        iid, L = m["ID"], int(m["len"])
        pos = int(m["site"].split(":")[1])
        tc, trng = m["target"].rsplit(":", 1)
        conv = {m["target"]: (tc, int(trng.split("-")[0]))}
        lseg = load(os.path.join(d, f"{iid}.local_segments.tsv"))
        lal, lcov, lwid = summary(lseg.get(iid, []), L, conv)
        gal, gcov, gwid = summary(gseg.get(iid, []), L)
        drop = gwid - lwid if lal else ""
        (pa, pid, uend), (sb, sid, ustart) = flank_test(d, iid, L, pos, a.site_tol)
        fcov = (min(L, pa + (L - sb))) / L if pa or sb < L else 0.0
        fid = ((pa * (pid or 0)) + ((L - sb) * (sid or 0))) / max(1, pa + (L - sb)) if (pid or sid) else 0.0
        gap = sb - pa  # > 0: insertion bases between the prefix and the suffix that match neither flank
        ok, unit, note = tandem_test(lal, pos, a.site_tol, a.max_gap_frac)
        drop_ok = lal and drop <= a.max_identity_drop
        if fcov >= a.min_cov and gap <= max(a.site_tol, a.max_gap_frac * L) and gwid - fid <= a.max_identity_drop:
            cls = "tandem duplication"
            if pa and sb < L and pa > sb:  # prefix and suffix overlap: shared bases at the junction, unit = insertion
                us, ue = pos - (L - sb) + 1, pos + sb
            else:
                us = ustart if ustart else pos + 1
                ue = uend if uend else pos
            unit = f"{m['site'].split(':')[0]}:{us}-{ue}"
            ok = True
            where = ("insertion point inside the unit" if pa and sb < L
                     else ("unit right of the site" if pa else "unit left of the site"))
            note = (f"direct test: bases 1-{pa} match from POS+1 ({100 * (pid or 0):.1f}%), bases {sb + 1}-{L} match up to "
                    f"POS ({100 * (sid or 0):.1f}%); {where}; unit {ue - us + 1} bp, {L / max(1, ue - us + 1):.2f} copies")
            lcov, lwid, drop = fcov, fid, gwid - fid
        elif ok and lcov >= a.min_cov and drop_ok:
            cls = "tandem duplication"
            ulen = int(unit.split(":")[1].split("-")[1]) - int(unit.split(":")[1].split("-")[0]) + 1
            note += f", {L / ulen:.2f} copies"
        elif lcov >= a.min_cov and drop_ok:
            cls = "complex SV"
        elif gcov >= a.min_cov:
            cls = "chimeric insertion"
        else:
            cls = "unresolved"

        def lay(al):
            return " | ".join(f"q{x['qs']}-{x['qe']}>{x['chrom']}:{x['s']}-{x['e']}({x['strand']},{100 * x['ident']:.1f}%"
                              f"{',' + sd_name(x['chrom'], x['s'], x['e']) if sd_name(x['chrom'], x['s'], x['e']) else ''})"
                              for x in al)
        rows.append((iid, m["site"], L, m["site_copies"], cls, f"{lcov:.3f}", f"{lwid:.4f}", f"{gcov:.3f}", f"{gwid:.4f}",
                     f"{drop:.4f}" if drop != "" else "", unit if ok else "", note, lay(lal), lay(gal)))
    with open(a.out, "w") as o:
        o.write("ID\tsite\tlen\tsite_sd\tclass\tlocal_cov\tlocal_identity\tgenome_cov\tgenome_identity\tidentity_drop\tunit\t"
                "note\tlocal_layout\tgenome_layout\n")
        for r in sorted(rows, key=lambda r: (r[4], int(r[0].split("-")[1]))):
            o.write("\t".join(str(v) for v in r) + "\n")
    print(collections.Counter(r[4] for r in rows))


if __name__ == "__main__":
    main()
