#!/usr/bin/env python3
"""Iterative alignment of query sequences with minimap2 (same logic as wdl/tools/AlignQueriesIterative.wdl): each piece's
primary alignment is cut out and the unaligned parts on either side (>= --min-piece-len) are realigned, up to
--max-rounds; pieces longer than --sr-max-len use the asm20 preset, shorter ones the sr preset. Needs minimap2 on PATH.
Output columns: ID, q_start, q_end, ref, strand, identity, mapq, aligned_len, preset, round, status.
"""
import argparse
import re
import subprocess

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--queries", required=True, help="Query FASTA (inserted sequences)")
ap.add_argument("--asm20-index", required=True, help="minimap2 -x asm20 index (.mmi) or FASTA of the target")
ap.add_argument("--sr-index", required=True, help="minimap2 -x sr index (.mmi) or FASTA of the target")
ap.add_argument("--min-piece-len", type=int, default=30, help="Pieces shorter than this are not realigned [30]")
ap.add_argument("--max-rounds", type=int, default=20, help="Maximum rounds [20]")
ap.add_argument("--sr-max-len", type=int, default=500, help="Pieces up to this length use the sr preset [500]")
ap.add_argument("--threads", default="2", help="minimap2 threads [2]")
ap.add_argument("--out", required=True, help="Output segments TSV")
a = ap.parse_args()
queries_fa, min_len, max_rounds, sr_max, threads, out = a.queries, a.min_piece_len, a.max_rounds, a.sr_max_len, a.threads, a.out
CIG = re.compile(r"(\d+)([MIDNSHP=X])")

seqs, name = {}, None
for line in open(queries_fa):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]
        seqs[name] = []
    elif line:
        seqs[name].append(line.upper())
seqs = {k: "".join(v) for k, v in seqs.items()}


def align(pieces, preset, index):
    best = {}
    if not pieces:
        return best
    tmp = out + ".pieces.tmp.fa"
    with open(tmp, "w") as o:
        for q, s, e in pieces:
            o.write(f">{q}|{s}|{e}\n{seqs[q][s:e]}\n")
    sam = subprocess.run(["minimap2", "-a", "-x", preset, "-t", threads, index, tmp],
                         check=True, capture_output=True, text=True).stdout
    for line in sam.splitlines():
        if line.startswith("@"):
            continue
        c = line.split("\t")
        flag = int(c[1])
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        ops = [(int(n), op) for n, op in CIG.findall(c[5])]
        lead = ops[0][0] if ops[0][1] in "SH" else 0
        tail = ops[-1][0] if ops[-1][1] in "SH" else 0
        m = sum(n for n, op in ops if op in "M=X")
        ins = sum(n for n, op in ops if op == "I")
        dl = sum(n for n, op in ops if op in "DN")
        qlen = sum(n for n, op in ops if op in "MIS=XH")
        nm = next((int(t[5:]) for t in c[11:] if t.startswith("NM:i:")), 0)
        rev = bool(flag & 16)
        qs, qe = (tail, qlen - lead) if rev else (lead, qlen - tail)
        rs = int(c[3])
        best[c[0]] = dict(qs=qs, qe=qe, chrom=c[2], rs=rs, re=rs + m + dl - 1, strand="-" if rev else "+",
                          ident=(m - (nm - ins - dl)) / max(1, m + ins + dl), mapq=int(c[4]), preset=preset)
    return best


# pieces still to align: (query ID, start, end), 0-based half-open on the inserted sequence
todo = [(q, 0, len(s)) for q, s in seqs.items() if len(s) >= min_len]
segments = [(q, 0, len(s), "", "", "", "", "", "", 0, "shorter than min_piece_len")
            for q, s in seqs.items() if len(s) < min_len]
for rnd in range(1, max_rounds + 1):
    if not todo:
        break
    best = align([p for p in todo if p[2] - p[1] > sr_max], "asm20", a.asm20_index)
    best.update(align([p for p in todo if p[2] - p[1] <= sr_max], "sr", a.sr_index))
    nxt = []
    for q, s, e in todo:
        b = best.get(f"{q}|{s}|{e}")
        if b is None:
            segments.append((q, s, e, "", "", "", "", "", "", rnd, "unaligned"))
            continue
        a0, a1 = s + b["qs"], s + b["qe"]
        segments.append((q, a0, a1, f"{b['chrom']}:{b['rs']}-{b['re']}", b["strand"], f"{b['ident']:.3f}", b["mapq"],
                         a1 - a0, b["preset"], rnd, "aligned"))
        for p0, p1 in ((s, a0), (a1, e)):
            if p1 - p0 >= min_len:
                nxt.append((q, p0, p1))
            elif p1 > p0:
                segments.append((q, p0, p1, "", "", "", "", "", "", rnd, "shorter than min_piece_len"))
    todo = nxt
for q, s, e in todo:
    segments.append((q, s, e, "", "", "", "", "", "", max_rounds, "not aligned within max_rounds"))
with open(out, "w") as o:
    o.write("ID\tq_start\tq_end\tref\tstrand\tidentity\tmapq\taligned_len\tpreset\tround\tstatus\n")
    for r in sorted(segments, key=lambda r: (r[0], r[1])):
        o.write(f"{r[0]}\t{r[1] + 1}\t{r[2]}\t" + "\t".join(str(v) for v in r[3:]) + "\n")
