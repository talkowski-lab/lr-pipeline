#!/usr/bin/env python3
"""Look for the source of inserted sequences in the reference around the insertion site (works without alignments).

For each ID, the window POS +/- max(--min-window, --window-x x insertion length) of the reference is scanned for the
inserted sequence's k-mers on both strands. Matched k-mers are chained into the largest cluster (consecutive
reference positions <= --gap apart, same strand); the cluster gives
  local_cov     fraction of the inserted bases covered by k-mers of that cluster
  local_source  reference interval of the cluster, strand
  local_dist    0 if the interval touches the insertion point (within --adjacent bp), else the distance to it
local_call:
  tandem duplication      local_cov >= --min-cov, + strand, source adjacent to the site
  inverted duplication    local_cov >= --min-cov, - strand, adjacent
  nearby duplication      local_cov >= --min-cov, not adjacent (dispersed within the window)
  partly local            --min-partial <= local_cov < --min-cov
  no local source         below; combined with the remap class when one is given (e.g. interspersed / novel)
Python standard library + pysam.
"""
import argparse
import collections
import csv
import gzip

import pysam


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


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
    ap.add_argument("--ids", required=True, help="File with one insertion ID per line")
    ap.add_argument("--sequences", required=True, help="Inserted sequences FASTA")
    ap.add_argument("--flanks", required=True, help="Flank table (extract_insertion_sequences.py; CHROM, POS)")
    ap.add_argument("--ref-fa", required=True, help="Reference FASTA with .fai")
    ap.add_argument("--annotations", default=None,
                    help="Optional build_insertion_annotations.R output (remap_class, remap_coords)")
    ap.add_argument("--k", type=int, default=13, help="k-mer size [13]")
    ap.add_argument("--min-window", type=int, default=5000, help="Minimum window on each side, bp [5000]")
    ap.add_argument("--window-x", type=float, default=3.0, help="Window = this x insertion length if larger [3]")
    ap.add_argument("--gap", type=int, default=200, help="Max reference gap between chained k-mers [200]")
    ap.add_argument("--adjacent", type=int, default=5, help="Source within this distance of POS is adjacent [5]")
    ap.add_argument("--min-cov", type=float, default=0.8, help="Coverage for a local duplication [0.8]")
    ap.add_argument("--min-partial", type=float, default=0.3, help="Coverage for 'partly local' [0.3]")
    ap.add_argument("--out", required=True, help="Output TSV")
    a = ap.parse_args()

    ids = [x.strip() for x in open(a.ids) if x.strip()]
    seqs = read_fasta(a.sequences)
    fl = {r["ID"]: r for r in csv.DictReader(open(a.flanks), delimiter="\t")}
    ann = {}
    if a.annotations:
        op = gzip.open if a.annotations.endswith(".gz") else open
        with op(a.annotations, "rt") as f:
            ann = {r["ID"]: r for r in csv.DictReader(f, delimiter="\t")}
    ref = pysam.FastaFile(a.ref_fa)
    k = a.k
    cols = ["ID", "len", "remap_class", "remap_coords", "local_call", "local_cov", "local_source", "local_strand", "local_dist",
            "source_len"]
    with open(a.out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for iid in ids:
            seq = seqs[iid]
            n = len(seq)
            ch, pos = fl[iid]["CHROM"], int(fl[iid]["POS"])
            w = max(a.min_window, int(a.window_x * n))
            ws = max(0, pos - w)
            win = ref.fetch(ch, ws, pos + w).upper()
            best = (0.0, None)
            for strand, q in (("+", seq), ("-", rc(seq))):
                qidx = collections.defaultdict(list)
                for qi in range(len(q) - k + 1):
                    qidx[q[qi:qi + k]].append(qi)
                # stream k-mer hits in reference order (no hit list kept); one coverage array per cluster
                cov, c0, last = None, None, None
                for rp in range(len(win) - k + 1):
                    qis = qidx.get(win[rp:rp + k])
                    if not qis:
                        continue
                    if cov is not None and rp - last > a.gap:
                        frac = sum(cov) / n
                        if frac > best[0]:
                            best = (frac, (strand, c0 + ws, last + ws + k))
                        cov = None
                    if cov is None:
                        cov, c0 = bytearray(len(q)), rp
                    last = rp
                    for qi in qis:
                        cov[qi:qi + k] = b"\x01" * k
                if cov is not None:
                    frac = sum(cov) / n
                    if frac > best[0]:
                        best = (frac, (strand, c0 + ws, last + ws + k))
            frac, info = best
            r = ann.get(iid, {})
            row = dict(ID=iid, len=n, remap_class=r.get("remap_class", ""), remap_coords=r.get("remap_coords", ""),
                       local_cov=round(frac, 3))
            if info:
                strand, s0, e0 = info
                dist = 0 if s0 - a.adjacent <= pos <= e0 + a.adjacent else min(abs(s0 - pos), abs(e0 - pos))
                row.update(local_source=f"{ch}:{s0 + 1}-{e0}", local_strand=strand, local_dist=dist, source_len=e0 - s0)
            if frac >= a.min_cov:
                if row["local_dist"] == 0:
                    row["local_call"] = "tandem duplication" if row["local_strand"] == "+" else "inverted duplication"
                else:
                    row["local_call"] = "nearby duplication"
            elif frac >= a.min_partial:
                row["local_call"] = "partly local"
            else:
                row["local_call"] = "no local source" + (f" (remap {row['remap_class']})" if row["remap_class"] else "")
            o.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


if __name__ == "__main__":
    main()
