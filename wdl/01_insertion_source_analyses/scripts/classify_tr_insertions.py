#!/usr/bin/env python3
"""STR / VNTR stage of the insertion annotation (after MEIs).

1. MEIs (ins_call starting with "MEI " in the annotation table, i.e. the VCF SUB_FAMILY calls) are excluded.
2. Insertions whose site (POS +/- --site-pad) touches a TR catalog locus are compared with the local repeat: the
   fraction of the inserted sequence's k-mers found in the reference TR sequence (catalog interval +/- --tr-pad bp)
   = tr_containment, taking the best of the overlapping loci. Also reported: motif coverage = fraction of inserted bases
   covered by exact copies of the locus motifs (any rotation, either strand).
     tr_containment >= --min-contain or motif_cov >= --min-motif-cov -> "STR expansion" / "VNTR expansion" (longest motif)
     tr_containment >= --min-partial or motif_cov >= --min-motif-partial
        and the insertion is itself a tandem repeat (self_repeat >= --min-self) -> "<STR|VNTR> expansion (diverged units)"
        otherwise                                                            -> "in TR, partly matching"
     below: tandem-repeat insertion -> "in TR, different tandem repeat"; otherwise -> "in TR, non-repeat sequence"
   (motif coverage rescues expansions whose units diverge from the reference copies, which lowers k-mer containment)
3. Insertions outside the catalog: if their best-scoring alignment (annotate_hit_repeats.R) lands >= --min-tr-frac in
   a catalog STR/VNTR, covers >= --min-perc of the insertion at identity >= --min-ident, the insertion is
     "STR/VNTR copy, local" when that alignment is within --local-bp of the site (an expansion just outside the
     catalogued interval), otherwise "interspersed STR/VNTR insertion".
   Remaining sequences that are themselves tandem repeats (>= --min-self of k-mers recur within 100 bp, or >= --min-copies copies
     of a period at which >= --min-period-ident of bases equal the base one period later, for sequences <= 5 kb) are
     "tandem-repeat sequence, uncatalogued", typed STR / VNTR by that period (else the most common distance between
     recurring 8-mers; STR <= --str-max bp).
STR = longest motif <= --str-max bp, else VNTR. Python standard library + pysam.
"""
import argparse
import bisect
import collections
import csv
import gzip

import pysam


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def kmers(s, k):
    return [s[i:i + k] for i in range(len(s) - k + 1)]


def motif_cov(seq, motifs):
    """Fraction of bases in tandem runs of a locus motif (any rotation, either strand): exact copies that abut or
    overlap are merged into runs, and only runs >= max(3 copies, 10 bp) count, so short motifs (A, AT) do not
    score on scattered single matches."""
    cov = bytearray(len(seq))
    for m in motifs:
        if not m or len(m) > len(seq):
            continue
        forms = {m[i:] + m[:i] for i in range(len(m))}
        forms |= {rc(f) for f in forms}
        occ = []
        for f in forms:
            start = seq.find(f)
            while start != -1:
                occ.append((start, start + len(f)))
                start = seq.find(f, start + 1)
        occ.sort()
        runs = []
        for s0, e0 in occ:
            if runs and s0 <= runs[-1][1]:
                runs[-1][1] = max(runs[-1][1], e0)
            else:
                runs.append([s0, e0])
        for s0, e0 in runs:
            if e0 - s0 >= max(3 * len(m), 10):
                for j in range(s0, e0):
                    cov[j] = 1
    return sum(cov) / len(seq) if seq else 0.0


def self_repeat(seq, k=12, win=100):
    last, rep = {}, 0
    n = len(seq) - k + 1
    for i in range(max(0, n)):
        km = seq[i:i + k]
        if km in last and i - last[km] <= win:
            rep += 1
        last[km] = i
    return rep / n if n > 0 else 0.0


def dominant_period(seq, k=8, max_period=500):
    """Most common distance between consecutive occurrences of the same k-mer (<= max_period): the repeat unit length
    of a tandem-repeat sequence. 0 when nothing recurs."""
    last, dist = {}, collections.Counter()
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if km in last and i - last[km] <= max_period:
            dist[i - last[km]] += 1
        last[km] = i
    return dist.most_common(1)[0][0] if dist else 0


def period_identity(seq, max_period=200, max_len=5000, min_ident=0.8):
    """Smallest period p (1..max_period, len >= 2p) at which >= min_ident of bases equal the base p positions later;
    if none, the period with the highest identity. Returns (period, identity); (0, 0) for sequences longer than max_len
    (the k-mer self-repeat test covers those)."""
    n = len(seq)
    if n > max_len or n < 4:
        return 0, 0.0
    best = (0, 0.0)
    for p in range(1, min(max_period, n // 2) + 1):
        ident = sum(1 for i in range(n - p) if seq[i] == seq[i + p]) / (n - p)
        if ident >= min_ident:
            return p, round(ident, 3)
        if ident > best[1]:
            best = (p, round(ident, 3))
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotations", required=True,
                    help="Per-insertion table from build_insertion_annotations.R (.tsv.gz)")
    ap.add_argument("--remap", required=True, help="Truvari remap TSV, no header (CHROM POS REF ALT ID ...)")
    ap.add_argument("--tr-catalog", required=True, help="TR catalog BED.gz (chrom start end ID=...;MOTIFS=...)")
    ap.add_argument("--ref-fa", required=True, help="Reference FASTA with .fai")
    ap.add_argument("--hit-repeats", required=True, help="Per-hit table from annotate_hit_repeats.R (.tsv.gz)")
    ap.add_argument("--site-pad", type=int, default=10, help="Bases around POS that must touch a catalog locus [10]")
    ap.add_argument("--tr-pad", type=int, default=20,
                    help="Bases added to each side of the TR when building its k-mers [20]")
    ap.add_argument("--k", type=int, default=11, help="k-mer size for containment [11]")
    ap.add_argument("--min-contain", type=float, default=0.5, help="Containment for an STR/VNTR expansion [0.5]")
    ap.add_argument("--min-partial", type=float, default=0.2, help="Containment for 'partly matching' [0.2]")
    ap.add_argument("--min-motif-cov", type=float, default=0.8, help="Motif coverage for an STR/VNTR expansion [0.8]")
    ap.add_argument("--min-motif-partial", type=float, default=0.5, help="Motif coverage for 'partly matching' [0.5]")
    ap.add_argument("--str-max", type=int, default=6, help="Longest motif counted as STR [6]")
    ap.add_argument("--min-tr-frac", type=float, default=0.5,
                    help="Min TR fraction of the best alignment (outside catalog) [0.5]")
    ap.add_argument("--min-perc", type=float, default=80, help="Min %% of the insertion the best alignment covers [80]")
    ap.add_argument("--min-ident", type=float, default=0.9, help="Min identity of the best alignment [0.9]")
    ap.add_argument("--local-bp", type=int, default=1000,
                    help="Best alignment within insertion length + this of the site = local [1000]")
    ap.add_argument("--min-self", type=float, default=0.5, help="Self-repeat fraction for an uncatalogued tandem repeat [0.5]")
    ap.add_argument("--min-copies", type=float, default=3,
                    help="Min copies of the period in the insertion for the period-identity rule (2 copies = tandem dup) [3]")
    ap.add_argument("--min-period-ident", type=float, default=0.8,
                    help="Period identity for a tandem repeat (bases equal to the base one period later) [0.8]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    ann = {}
    with gzip.open(a.annotations, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            ann[r["ID"]] = r
    seqs = {}
    with open(a.remap) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            seqs[c[4]] = (c[0], int(c[1]), c[3].upper()[1:])
    chroms = {v[0] for v in seqs.values()}
    tr = collections.defaultdict(list)
    with gzip.open(a.tr_catalog, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[0] not in chroms:
                continue
            info = dict(x.split("=", 1) for x in c[3].split(";") if "=" in x)
            motifs = info.get("MOTIFS", "").split(",")
            tr[c[0]].append((int(c[1]), int(c[2]), info.get("ID", "."), motifs))
    for ch in tr:
        tr[ch].sort()
    starts = {ch: [x[0] for x in tr[ch]] for ch in tr}
    best = {}
    with gzip.open(a.hit_repeats, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["ID"] not in best or float(r["score"]) > float(best[r["ID"]]["score"]):
                best[r["ID"]] = r
    ref = pysam.FastaFile(a.ref_fa)

    cols = ["ID", "len", "remap_class", "ins_call", "tr_stage", "tr_id", "tr_interval", "tr_motifs", "tr_type",
            "tr_containment", "motif_cov", "self_repeat", "period", "period_identity", "best_hit", "best_hit_perc",
            "best_hit_identity", "best_hit_tr_frac", "best_hit_dist"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for iid, (chrom, pos, seq) in seqs.items():
            r = ann.get(iid, {})
            call = r.get("ins_call", "")
            row = dict(ID=iid, len=len(seq), remap_class=r.get("remap_class", ""), ins_call=call)
            row["self_repeat"] = round(self_repeat(seq), 3)
            pi_p, pi_id = period_identity(seq)
            row["period"] = pi_p if pi_id >= a.min_period_ident else dominant_period(seq)
            row["period_identity"] = pi_id
            if call.startswith("MEI "):
                row["tr_stage"] = "MEI (excluded)"
            else:
                lst = tr.get(chrom, [])
                i = bisect.bisect_right(starts.get(chrom, []), pos + a.site_pad)
                hits = [x for x in lst[max(0, i - 200):i] if x[1] >= pos - a.site_pad]
                cand = []
                ks = set(kmers(seq, a.k))
                for s, e, tid, motifs in hits:
                    tseq = ref.fetch(chrom, max(0, s - a.tr_pad), e + a.tr_pad).upper()
                    tk = set(kmers(tseq, a.k))
                    cont = len(ks & tk) / len(ks) if ks else 0.0
                    cand.append((cont, s, e, tid, motifs))
                if cand:
                    cont, s, e, tid, motifs = max(cand)
                    ttype = "STR" if max(len(m) for m in motifs) <= a.str_max else "VNTR"
                    row.update(tr_id=tid, tr_interval=f"{chrom}:{s + 1}-{e}", tr_motifs=",".join(motifs), tr_type=ttype,
                               tr_containment=round(cont, 3), motif_cov=round(motif_cov(seq, motifs), 3))
                    mc = row["motif_cov"]
                    if cont >= a.min_contain or mc >= a.min_motif_cov:
                        row["tr_stage"] = f"{ttype} expansion"
                    elif cont >= a.min_partial or mc >= a.min_motif_partial:
                        # partly matching the local repeat; if the insertion is itself a tandem repeat, its units have
                        # diverged from the reference copies
                        if row["self_repeat"] >= a.min_self:
                            row["tr_stage"] = f"{ttype} expansion (diverged units)"
                        else:
                            row["tr_stage"] = "in TR, partly matching"
                    elif row["self_repeat"] >= a.min_self:
                        row["tr_stage"] = "in TR, different tandem repeat"
                    else:
                        row["tr_stage"] = "in TR, non-repeat sequence"
                else:
                    b = best.get(iid)
                    if b:
                        trf = max(float(b["f_STR"]), float(b["f_VNTR"]))
                        row.update(best_hit=b["coords"], best_hit_perc=b["perc"], best_hit_identity=b["identity"],
                                   best_hit_tr_frac=round(trf, 3))
                        bc = b["coords"].split(":")[0]
                        bs, be = (int(x) for x in b["coords"].split(":")[1].split("-"))
                        if bc != chrom:
                            dist = -1
                        else:
                            dist = 0 if bs <= pos <= be else min(abs(bs - pos), abs(be - pos))
                        row["best_hit_dist"] = dist
                        if trf >= a.min_tr_frac and float(b["perc"]) >= a.min_perc and float(b["identity"]) >= a.min_ident:
                            local = dist >= 0 and dist <= len(seq) + a.local_bp
                            kind = "STR" if float(b["f_STR"]) >= float(b["f_VNTR"]) else "VNTR"
                            if local:
                                row["tr_stage"] = f"{kind} copy, local (outside catalog interval)"
                            else:
                                row["tr_stage"] = f"interspersed {kind} insertion"
                    if "tr_stage" not in row:
                        periodic = row["period_identity"] >= a.min_period_ident and len(seq) >= a.min_copies * row["period"] > 0
                        is_tr = row["self_repeat"] >= a.min_self or periodic
                        row["tr_stage"] = "tandem-repeat sequence, uncatalogued" if is_tr else "not TR"
                        if row["tr_stage"] != "not TR":
                            row["tr_type"] = "STR" if 0 < row["period"] <= a.str_max else "VNTR"
            o.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
