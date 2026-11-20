#!/usr/bin/env python3
"""Final annotation of short insertions (e.g. 20-49 bp, below the remap minimum), same columns as build_final_annotation.R
(ID, len, final_group, final_class, final_detail). Stages, first call wins:
  1. TR stage (classify_tr_insertions.py): STR / VNTR expansion -> STR / VNTR; uncatalogued tandem repeat -> STR_novel /
     VNTR_novel; in TR, partly matching / different tandem repeat -> other TR (group "STR / VNTR")
  2. local copy in hg38 (small_insertion_local_dups.py), group "duplication": tandem duplication (copy at the site, origin
     from the split; ", insertion point inside the unit" kept in the detail), inverted tandem duplication, copy at the site
     with indels -> tandem duplication, nearby copy more than 5 bp away -> interspersed duplication (strand in the detail)
  3. T2T (small_insertion_t2t_check.py), group "present in T2T": T2T_<RefSeq Curated genes on hs1 overlapping the best
     T2T hit, or T2T_<chrom>:<start>-<end>> (LOC* and -AS genes left out); the detail says whether the T2T copy is at the
     T2T site, near it or elsewhere
     A best T2T hit on chrM is a NUMT instead (group and class NUMT, chrM region, strand and identity in the detail)
  4. the rest (no local hg38 copy, no T2T copy): unresolved
"""
import argparse
import collections
import csv
import gzip
import re


def read(path):
    return {r["ID"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def t2t_label(hit, genes, max_genes):
    m = re.match(r"^([^:]+):(\d+)-(\d+)", hit)
    ch, s, e = m.group(1), int(m.group(2)), int(m.group(3))
    names = []
    for gs, ge, n in genes.get(ch, []):
        if ge > s - 1 and gs < e and n and not n.startswith("LOC") and "-AS" not in n and n not in names:
            names.append(n)
    return "T2T_" + ("/".join(names[:max_genes]) if names else f"{ch}:{s}-{e}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tr-stage", required=True, help="classify_tr_insertions.py output (all insertions)")
    ap.add_argument("--remap", required=True, help="Remap-style TSV (CHROM POS REF ALT ID): sites for the tandem origins")
    ap.add_argument("--local-dups", required=True, help="small_insertion_local_dups.py output")
    ap.add_argument("--t2t", required=True, help="small_insertion_t2t_check.py output")
    ap.add_argument("--genes-table", required=True,
                    help="T2T gene table (genePred-extended .txt.gz: chrom col 3, txStart 5, txEnd 6, name2 13)")
    ap.add_argument("--max-genes", type=int, default=3, help="Max gene names in a T2T label [3]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    tr, ld, t2 = read(a.tr_stage), read(a.local_dups), read(a.t2t)
    site = {}
    for line in open(a.remap):
        c = line.split("\t")
        site[c[4]] = (c[0], int(c[1]))
    genes = collections.defaultdict(list)
    with gzip.open(a.genes_table, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            genes[c[2]].append((int(c[4]), int(c[5]), c[12]))

    n = collections.Counter()
    with open(a.out, "w") as o:
        o.write("ID\tlen\tfinal_group\tfinal_class\tfinal_detail\n")
        for i, r in tr.items():
            st = r["tr_stage"]
            if st.startswith("STR expansion") or st.startswith("VNTR expansion"):
                grp, cls = "STR / VNTR", st.split(" ")[0]
                det = f"{st}; {r['tr_id']}"
            elif st == "tandem-repeat sequence, uncatalogued":
                grp, cls = "STR / VNTR", (r["tr_type"] or "VNTR") + "_novel"
                det = f"not in the TR catalog; period {r['period']} bp"
            elif st in ("in TR, partly matching", "in TR, different tandem repeat"):
                grp, cls, det = "STR / VNTR", "other TR", f"{st}; {r['tr_id']}"
            else:
                d = ld[i]
                L = int(d["len"])
                ch, pos = site[i]
                t = t2.get(i)
                if d["call"].startswith("tandem duplication"):
                    k = int(d["tandem_split"])
                    grp, cls = "duplication", "tandem duplication"
                    det = (f"origin {ch}:{pos - k + 1}-{pos + L - k} (+), 100% of the insertion at "
                           f"{100 * float(d['tandem_identity']):.1f}% identity, 0 bp from the site"
                           + (", insertion point inside the unit" if 0 < k < L else ""))
                elif d["call"] == "inverted tandem duplication":
                    grp, cls = "duplication", "inverted tandem duplication"
                    det = f"origin at the site (-), {100 * float(d['inv_tandem_identity']):.1f}% identity"
                elif d["call"] == "copy at the site with indels":
                    grp, cls = "duplication", "tandem duplication"
                    det = (f"origin {d['nearby_copy']} ({d['nearby_strand']}), {100 * float(d['nearby_cov']):.0f}% of the "
                           f"insertion at {100 * float(d['nearby_identity']):.1f}% identity with indels, 0 bp from the site")
                elif d["call"].endswith("interspersed duplication (nearby)"):
                    grp, cls = "duplication", "interspersed duplication"
                    det = (f"origin {d['nearby_copy']} ({d['nearby_strand']}), {100 * float(d['nearby_cov']):.0f}% of the "
                           f"insertion at {100 * float(d['nearby_identity']):.1f}% identity, {d['nearby_dist']} bp from the site")
                elif t is not None and t["best_hit"].startswith("chrM:"):
                    grp, cls = "NUMT", "NUMT"
                    m = re.match(r"^(chrM:\d+-\d+)\((.)\)", t["best_hit"])
                    det = (f"origin {m.group(1)} ({m.group(2)}), {100 * float(t['best_cov']):.0f}% of the insertion at "
                           f"{100 * float(t['best_ident']):.1f}% identity; T2T-CHM13 BLAT")
                elif t is not None and t["call"] != "no T2T copy":
                    grp, cls = "present in T2T", t2t_label(t["best_hit"], genes, a.max_genes)
                    where = {"present in T2T at the site": "at the T2T site",
                             "copy near the T2T site": f"{t['dist_to_site']} bp from the T2T site",
                             "copy elsewhere in T2T": "away from the T2T site" if t["t2t_site"] != "not placed"
                             else "T2T site not placed"}[t["call"]]
                    det = (f"present in T2T-CHM13: {t['best_hit']}, {100 * float(t['best_cov']):.0f}% of the insertion at "
                           f"{100 * float(t['best_ident']):.1f}% identity, {where}; {t['n_ins_hits']} T2T hits; "
                           "no local hg38 copy")
                else:
                    grp, cls, det = "unresolved", "unresolved", "no local hg38 copy; no T2T copy"
                if st == "in TR, non-repeat sequence":
                    det += f"; in TR {r['tr_id']}"
            n[(grp, cls if grp != "present in T2T" else "T2T_*")] += 1
            o.write(f"{i}\t{r['len']}\t{grp}\t{cls}\t{det}\n")
    for k, v in sorted(n.items(), key=lambda x: -x[1]):
        print(v, *k, sep="\t")


if __name__ == "__main__":
    main()
