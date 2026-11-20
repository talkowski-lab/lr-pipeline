#!/usr/bin/env python3
"""MEI stage of the insertion annotation, from RepeatMasker run on the inserted sequences.

Per insertion:
  rm_chain        ordered RepeatMasker elements along the inserted sequence, name(+/-), simple repeats and
                  low-complexity pieces < 30 bp left out
  te_family       dominant mobile-element group by bp: Alu (incl. FLAM / FRAM monomers), L1, SVA, HERVK
                  (LTR5* / HERVK*), other_TE; te_frac = its bp / (insertion length - tail), tail = 3' poly-A
                  for + strand elements, 5' poly-T (with up to 30 bp before it) for - strand elements
  te_subfamily    RepeatMasker name of the largest element of that group
  te_identity     1 - RepeatMasker divergence (Kimura-free % div / 100) of that largest element to its consensus
  cons_start/end, cons_len, completeness ("full-length" / "5'-truncated" / "3'-truncated" / "internal fragment")
                  and kept_parts (approximate consensus regions, below) for the dominant element, merging the
                  fragments of the same RepeatMasker element ID
  l1_5prime_inversion  L1 pieces of one element in opposite orientations (twin priming)
  tsd_bp, polyA   target-site duplication from the reference flanks (insertions are left-aligned, so the TSD is the
                  start of the inserted sequence = start of the right flank, or its end = end of the left flank)
                  and a poly-A (or poly-T at the 5' end for minus-strand elements)
  site_in_ref     family of a reference element (RepeatMasker-style BED) at the insertion site, if any
Consensus regions (approximate): Alu left arm 1-42%, right arm 47-100% of the consensus; L1 5'UTR 0-15%, ORF1 15-32%,
ORF2 33-96%, 3'UTR 96-100%; SVA 5' hexamer/Alu-like 1-435, VNTR 436-870, SINE-R 871-end of the RepeatMasker
consensus (VNTR collapsed); HERV-K: LTR5* = LTR, HERVK-int = internal.
mei_stage (first match wins):
  VCF MEI                     "VCF MEI <family>" (re-evaluated: RepeatMasker family, subfamily, completeness, parts);
                              flagged when the site is inside a same-family reference element and the sequence is
                              VNTR-like (possible VNTR variant rather than a new insertion)
  VNTR variant within SVA     site inside a reference SVA and the insertion GC-rich (>= --sva-vntr-gc) or a TR expansion
  TR variant within L1        site inside a reference L1, the insertion a TR expansion (classify_tr_insertions.py)
                              and >= 50% L1 sequence (an L1-derived tandem repeat)
  new <family> insertion      dominant MEI family covers >= --min-te-frac, a TSD of --min-tsd to 29 bp (a match
                              reaching the 30 bp search limit is a tandem duplication, not a TSD) or a poly-A, and the
                              site is not inside a same-family reference element
  Alu chain / nested Alus     >= 2 Alu-family elements (incl. FLAM / FRAM) covering >= --min-chain-frac
  <family> sequence without insertion hallmarks
                              dominant MEI family covers >= --min-te-frac but no TSD / poly-A or the site is inside a
                              same-family reference element (duplications / copies of existing elements)
  (empty)                     left for the TR / duplication stages
"""
import argparse
import bisect
import collections
import csv
import gzip
import re

MEI = ("Alu", "L1", "SVA", "HERVK")


def group_of(name, cls):
    if cls.startswith("SINE/Alu"):
        return "Alu"
    if cls.startswith("LINE/L1"):
        return "L1"
    if cls.startswith("Retroposon/SVA"):
        return "SVA"
    if cls.startswith("LTR/ERVK") and re.match(r"^(LTR5|HERVK)", name):
        return "HERVK"
    if cls.split("/")[0] in ("SINE", "LINE", "LTR", "DNA", "Retroposon", "RC"):
        return "other_TE"
    return "other"


def parts(fam, name, s, e, L):
    if fam == "Alu":
        if name.startswith("FLAM"):
            return ["left monomer"]
        if name.startswith("FRAM"):
            return ["right monomer"]
        reg = [("left arm", 0, 0.42 * L), ("right arm", 0.47 * L, L)]
    elif fam == "L1":
        reg = [("5'UTR", 0, 0.15 * L), ("ORF1", 0.15 * L, 0.32 * L), ("ORF2", 0.33 * L, 0.96 * L), ("3'UTR", 0.96 * L, L)]
    elif fam == "SVA":
        reg = [("5' hexamer/Alu-like", 0, 435), ("VNTR", 435, 870), ("SINE-R", 870, L)]
    elif fam == "HERVK":
        return ["LTR"] if name.startswith("LTR5") else ["internal (gag-pro-pol-env)"]
    else:
        return []
    return [n for n, a, b in reg if min(e, b) - max(s, a) >= min(30, 0.3 * (b - a))]


def completeness(fam, s, e, L, qspan):
    tol5 = {"Alu": 15, "L1": 100, "SVA": 40}.get(fam, 30)
    tol3 = {"Alu": 20, "L1": 60, "SVA": 40}.get(fam, 30)
    full5, full3 = s <= tol5, e >= L - tol3
    if fam == "L1" and full5 and full3 and qspan < 5500:
        full5 = False
    if full5 and full3:
        return "full-length"
    if full3:
        return "5'-truncated"
    if full5:
        return "3'-truncated"
    return "internal fragment"


def read_rm(path):
    els = collections.defaultdict(list)
    with open(path) as f:
        for i, line in enumerate(f):
            c = line.split()
            if i < 3 or len(c) < 15:
                continue
            q, qs, qe = c[4], int(c[5]), int(c[6])
            strand = "+" if c[8] == "+" else "-"
            name, cls = c[9], c[10]
            if strand == "+":
                cs, ce, left = int(c[11]), int(c[12]), int(c[13].strip("()"))
            else:
                left, ce, cs = int(c[11].strip("()")), int(c[12]), int(c[13])
            els[q].append(dict(qs=qs, qe=qe, strand=strand, name=name, cls=cls, cs=cs, ce=ce, L=ce + left,
                               rmid=c[14], div=float(c[1]), grp=group_of(name, cls)))
    for q in els:
        els[q].sort(key=lambda d: d["qs"])
    return els


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rm-out", required=True, help="RepeatMasker .out of the inserted sequences")
    ap.add_argument("--sequences", required=True, help="Inserted sequences FASTA (extract_insertion_sequences.py)")
    ap.add_argument("--flanks", required=True, help="Flank table (extract_insertion_sequences.py)")
    ap.add_argument("--annotations", required=True,
                    help="build_insertion_annotations.R output (.tsv.gz; ins_call, remap_class)")
    ap.add_argument("--tr-stage", required=True, help="classify_tr_insertions.py output")
    ap.add_argument("--repeats", required=True, help="RepeatMasker-style reference BED.gz (ref beg end name class family ...)")
    ap.add_argument("--site-pad", type=int, default=10, help="Bases around POS that must touch a reference element [10]")
    ap.add_argument("--min-te-frac", type=float, default=0.8, help="Dominant MEI family fraction for new / copy calls [0.8]")
    ap.add_argument("--min-chain-frac", type=float, default=0.7, help="Alu-family fraction for an Alu chain [0.7]")
    ap.add_argument("--min-tsd", type=int, default=7, help="Min TSD length counted as a hallmark [7]")
    ap.add_argument("--sva-vntr-gc", type=float, default=0.75, help="GC of an SVA VNTR variant [0.75]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    els = read_rm(a.rm_out)
    seqs = read_fasta(a.sequences)
    flanks = {r["ID"]: r for r in csv.DictReader(open(a.flanks), delimiter="\t")}
    with gzip.open(a.annotations, "rt") as f:
        ann = {r["ID"]: r for r in csv.DictReader(f, delimiter="\t")}
    trs = {r["ID"]: r for r in csv.DictReader(open(a.tr_stage), delimiter="\t")}
    chroms = {r["CHROM"] for r in flanks.values()}
    ref = collections.defaultdict(list)
    with gzip.open(a.repeats, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] in chroms and c[5] in ("Alu", "L1", "SVA", "ERVK"):
                ref[c[0]].append((int(c[1]), int(c[2]), {"ERVK": "HERVK"}.get(c[5], c[5])))
    for ch in ref:
        ref[ch].sort()
    rstarts = {ch: [x[0] for x in ref[ch]] for ch in ref}

    cols = ["ID", "len", "ins_call", "remap_class", "tr_stage", "mei_stage", "te_family", "te_frac", "te_subfamily", "te_strand",
            "te_identity",
            "cons_start", "cons_end", "cons_len", "completeness", "kept_parts", "l1_5prime_inversion", "tsd_bp", "tsd_seq",
            "polyA",
            "site_in_ref", "gc", "n_alu_like", "rm_chain"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for iid, fl in flanks.items():
            seq = seqs[iid]
            n = len(seq)
            r = ann.get(iid, {})
            row = dict(ID=iid, len=n, ins_call=r.get("ins_call", ""), remap_class=r.get("remap_class", ""),
                       tr_stage=trs.get(iid, {}).get("tr_stage", ""))
            row["gc"] = round((seq.count("G") + seq.count("C")) / max(1, n), 3)
            e = els.get(iid, [])
            chain = [d for d in e if not (d["cls"] in ("Simple_repeat", "Low_complexity") and d["qe"] - d["qs"] < 30)]
            row["rm_chain"] = "|".join(f"{d['name']}({d['strand']})" for d in chain)
            bp = collections.Counter()
            for d in e:
                bp[d["grp"]] += d["qe"] - d["qs"] + 1
            # poly-A / poly-T
            pa = re.search(r"A{10,}[ACGT]{0,30}$", seq) is not None
            pt = re.match(r"^[ACGT]{0,30}T{10,}", seq) is not None
            mei_bp = {g: bp[g] for g in MEI if bp[g] > 0}
            fam = max(mei_bp, key=mei_bp.get) if mei_bp else ("other_TE" if bp["other_TE"] else "")
            # the tail is a 3' poly-A for + strand elements and a 5' poly-T for - strand elements
            fam_strand = max((d for d in e if d["grp"] == fam), key=lambda d: d["qe"] - d["qs"])["strand"] if fam in MEI else "+"
            tail = re.search(r"A{8,}[ACGT]{0,30}$", seq) if fam_strand == "+" else re.match(r"^[ACGT]{0,30}?T{8,}", seq)
            denom = n - (len(tail.group()) if tail else 0)
            row["te_family"] = fam
            row["te_frac"] = round(min(1.0, bp[fam] / max(1, denom)), 3) if fam else 0
            row["n_alu_like"] = sum(1 for d in e if d["grp"] == "Alu")
            if fam in MEI:
                fe = [d for d in e if d["grp"] == fam]
                big = max(fe, key=lambda d: d["qe"] - d["qs"])
                same = [d for d in fe if d["rmid"] == big["rmid"]]
                s, en, L = min(d["cs"] for d in same), max(d["ce"] for d in same), big["L"]
                qspan = sum(d["qe"] - d["qs"] + 1 for d in same)
                row.update(te_subfamily=big["name"], te_strand=big["strand"], te_identity=round(1 - big["div"] / 100, 3),
                           cons_start=s, cons_end=en, cons_len=L,
                           completeness=completeness(fam, s, en, L, qspan),
                           kept_parts=",".join(parts(fam, big["name"], s, en, L)))
                strands = {d["strand"] for d in same}
                row["l1_5prime_inversion"] = fam == "L1" and len(strands) == 2
                pa_dir = pa if big["strand"] == "+" else pt
            else:
                pa_dir = pa or pt
            row["polyA"] = pa_dir
            ins, lf, rf = seq, fl["left_flank"], fl["right_flank"]
            ks = [k for k in range(5, 31) if ins[:k] == rf[:k] or ins[-k:] == lf[-k:]]
            k = max(ks) if ks else 0
            row["tsd_bp"] = k
            row["tsd_seq"] = (ins[:k] if ins[:k] == rf[:k] else ins[-k:]) if k else ""
            ch, pos = fl["CHROM"], int(fl["POS"])
            i = bisect.bisect_right(rstarts.get(ch, []), pos + a.site_pad)
            fams = sorted({x[2] for x in ref.get(ch, [])[max(0, i - 50):i] if x[1] >= pos - a.site_pad})
            row["site_in_ref"] = ",".join(fams)

            tr = row["tr_stage"]
            tr_exp = "expansion" in tr
            stage = ""
            if row["ins_call"].startswith("MEI "):
                vf = row["ins_call"][4:]
                stage = f"VCF MEI {vf}"
                if vf in fams and (row["gc"] >= a.sva_vntr_gc or tr_exp):
                    stage += " (site inside a reference " + vf + "; VNTR-like: check)"
            elif "SVA" in fams and (row["gc"] >= a.sva_vntr_gc or tr_exp):
                stage = "VNTR variant within SVA"
            elif "L1" in fams and tr_exp and fam == "L1" and row["te_frac"] >= 0.5:
                stage = "TR variant within L1"
            elif fam in MEI and row["te_frac"] >= a.min_te_frac and (a.min_tsd <= k < 30 or pa_dir) and fam not in fams:
                stage = f"new {fam} insertion (not in VCF)"
            elif row["n_alu_like"] >= 2 and bp["Alu"] / max(1, denom) >= a.min_chain_frac:
                stage = "Alu chain / nested Alus"
            elif fam in MEI and row["te_frac"] >= a.min_te_frac:
                stage = f"{fam} sequence without insertion hallmarks"
            row["mei_stage"] = stage
            o.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
