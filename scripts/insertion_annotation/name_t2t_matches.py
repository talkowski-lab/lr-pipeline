#!/usr/bin/env python3
"""Label T2T-matched insertions (t2t_match_calls.R) as T2T_<name>: the RefSeq Curated genes on T2T-CHM13 (UCSC hs1,
public API, coordinates only) overlapping the best T2T alignment (LOC* and antisense -AS genes left out, up to
--max-genes names joined by "/"), or T2T_<chrom>:<start>-<end> when no gene overlaps. Insertions that are not
"present in T2T-CHM13" get no label. With --genes-table a local copy of the gene table is used instead of the API.
"""
import argparse
import collections
import csv
import gzip
import json
import re
import time
import urllib.request


def genes_at(genome, track, chrom, start0, end):
    u = f"https://api.genome.ucsc.edu/getData/track?genome={genome};track={track};chrom={chrom};start={start0};end={end}"
    for attempt in range(3):
        try:
            req = urllib.request.Request(u, headers={"User-Agent": "insertion-annotation-script"})
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.load(r)
            time.sleep(0.3)
            return [g.get("name2", "") for g in d.get(track, [])]
        except Exception:
            time.sleep(2 * (attempt + 1))
    return []


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", required=True, help="t2t_match_calls.R <prefix>.t2t_calls.tsv")
    ap.add_argument("--genome", default="hs1", help="UCSC genome of the T2T assembly [hs1]")
    ap.add_argument("--track", default="ncbiRefSeqCurated", help="UCSC gene track [ncbiRefSeqCurated]")
    ap.add_argument("--max-genes", type=int, default=3, help="Max gene names in a label [3]")
    ap.add_argument("--genes-table", default=None,
                    help="Optional local gene table (genePred-extended .txt.gz: chrom col 3, txStart 5, txEnd 6, name2 13)")
    ap.add_argument("--out", required=True, help="Output TSV (input columns + t2t_genes, t2t_label)")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.calls), delimiter="\t"))
    local = None
    if a.genes_table:
        local = collections.defaultdict(list)
        with gzip.open(a.genes_table, "rt") as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                local[c[2]].append((int(c[4]), int(c[5]), c[12]))
    with open(a.out, "w") as o:
        cols = list(rows[0].keys()) + ["t2t_genes", "t2t_label"]
        o.write("\t".join(cols) + "\n")
        for r in rows:
            r["t2t_genes"], r["t2t_label"] = "", ""
            if r["t2t_call"].startswith("present in T2T"):
                m = re.match(r"^([^:]+):(\d+)-(\d+)", r["T2T_best_hit"])
                ch, s, e = m.group(1), int(m.group(2)), int(m.group(3))
                names = []
                hits = ([n for gs, ge, n in local.get(ch, []) if ge > s - 1 and gs < e] if local is not None
                        else genes_at(a.genome, a.track, ch, s - 1, e))
                for g in hits:
                    if g and not g.startswith("LOC") and "-AS" not in g and g not in names:
                        names.append(g)
                r["t2t_genes"] = ",".join(names)
                r["t2t_label"] = "T2T_" + ("/".join(names[:a.max_genes]) if names else f"{ch}:{s}-{e}")
            o.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
