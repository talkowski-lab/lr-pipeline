#!/usr/bin/env python3
"""Recurrence (dot) plot data for one long read against the reference region it aligns to.

1. Fetch the read from a BAM/CRAM (local or gs://) by name: within --region if given (fast, indexed), otherwise by
   scanning the whole file (fine for a small BAM, e.g. one already subset to the read). Its supplementary alignments
   are found through the SA tag, wherever they are. The full read sequence comes from the alignment that is not
   hard-clipped (normally the primary).
2. Chop to --region: the read's alignments overlapping --region are walked base by base through their CIGARs to find
   the read positions aligned to the region's first and last reference bases; the read segment between them (+/- --pad,
   union over alignments) is windowed, and the reference segment is --region +/- --pad. So for a long contig only the
   part aligned to the region is analysed. Without --region: all alignments on the primary's chromosome and their full
   spans. --read-range / --ref-region override either segment. Read positions in the outputs stay in original read
   coordinates. Sequence from --ref (bgzipped or plain FASTA with .fai;
   may also hold just a segment named chrom:start-end, e.g. from samtools faidx, that contains the reference segment).
3. Cut the read into consecutive --window bp windows (--step defaults to the window, so windows tile the read; a
   smaller --step makes them overlap). For each window find every place in the reference segment where it matches
   with edit distance <= --max-ed-frac x window (forward and reverse complement). Seed and verify, exact: the window
   is split into max_ed + 1 pieces; any match with <= max_ed edits contains at least one piece unedited (pigeonhole),
   so exact piece hits (sorted k-mer index of the reference) give every candidate, and a semi-global edit-distance DP
   is run only around them. Matches within --min-sep bp of a better one on the same strand are dropped, so each copy
   of a repeated sequence gives one hit per window.

Outputs:
  <out>.hits.tsv    read_start read_end strand ref_start ref_end edit_distance identity (ref 1-based, inclusive)
  <out>.aln.tsv     the read's alignment blocks: alignment (primary/supplementary), strand, read_start read_end chrom
                    ref_start ref_end (original read coordinates; CIGAR split at indels >= 50 bp; blocks >= 50 bp)
  <out>.info.tsv    read name, read length, reference segment, settings
Plot with plot_read_ref_recurrence.R.

Usage:
  read_ref_recurrence.py --bam sample.bam --read-name m84011_..../ccs [--region chr16:29400000-30300000] \\
      --ref Homo_sapiens_assembly38.fasta --out read1 [--window 100] [--step 100] [--max-ed-frac 0.1] [--pad 5000]
      [--ref-region chr16:29300000-30400000] [--min-sep 50] [--requester-pays-project PROJECT]
For gs:// inputs the OAuth token comes from $GCS_OAUTH_TOKEN or `gcloud auth print-access-token`, the requester-pays
project from --requester-pays-project, $GCS_REQUESTER_PAYS_PROJECT or $GOOGLE_PROJECT.
"""
import argparse
import bisect
import multiprocessing as mp
import os
import subprocess
import sys

import numpy as np
import pysam

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def encode(s):
    return np.frombuffer(s.upper().encode(), dtype=np.uint8)


def semiglobal_end_distances(pat, ref):
    """Edit distance of `pat` (length m) aligned so that it ends at each reference position j (any start): last DP row.
    Horizontal moves are resolved with a running minimum: D[j] = min_k<=j (V[k] + j - k)."""
    m, n = len(pat), len(ref)
    idx = np.arange(n + 1, dtype=np.int32)
    prev = np.zeros(n + 1, dtype=np.int32)  # row 0: free start anywhere in the reference
    for i in range(1, m + 1):
        cost = (ref != pat[i - 1]).astype(np.int32)
        diag = np.empty(n + 1, dtype=np.int32)
        diag[0] = i
        diag[1:] = prev[:-1] + cost
        v = np.minimum(diag, prev + 1)
        v[0] = i
        cur = np.minimum.accumulate(v - idx) + idx
        prev = cur
    return prev[1:]  # distance of the best alignment ending at reference position j (0-based)


BASE = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    BASE[_b] = _i


def kmer_codes(seq, k):
    """2-bit code of every k-mer of an encoded sequence; -1 where the k-mer contains a non-ACGT base."""
    b = BASE[seq]
    n = len(b) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.int64)
    code = np.zeros(n, dtype=np.int64)
    bad = np.zeros(n, dtype=bool)
    for t in range(k):
        code = code * 4 + np.minimum(b[t:t + n], 3)
        bad |= b[t:t + n] == 4
    code[bad] = -1
    return code


class RefIndex:
    def __init__(self, ref, k):
        self.ref, self.k = ref, k
        code = kmer_codes(ref, k)
        self.order = np.argsort(code, kind="stable")
        self.codes = code[self.order]

    def positions(self, piece):
        c = kmer_codes(piece, self.k)
        if not len(c) or c[0] < 0:
            return np.empty(0, dtype=np.int64)
        lo, hi = np.searchsorted(self.codes, c[0], "left"), np.searchsorted(self.codes, c[0], "right")
        return self.order[lo:hi]


def dedup(hits, min_sep):
    """Keep the best hit (lowest distance) and drop hits whose end is within min_sep of an already kept one."""
    kept = []
    for j, d in sorted(hits, key=lambda x: (x[1], x[0])):
        i = bisect.bisect_left(kept, j)
        if (i == len(kept) or kept[i] - j > min_sep) and (i == 0 or j - kept[i - 1] > min_sep):
            kept.insert(i, j)
    best = dict(hits)  # ends are unique: candidate intervals are merged before the DP
    return [(j, best[j]) for j in kept]


def hits_for(pat, idx, max_ed, min_sep):
    """All reference ends j (0-based) where pat matches with edit distance <= max_ed, via exact seeds + local DP."""
    m, k, ref = len(pat), idx.k, idx.ref
    cand = [idx.positions(pat[o:o + k]) - o for o in range(0, (max_ed + 1) * k, k) if o + k <= m]
    cand = np.unique(np.concatenate(cand)) if cand else np.empty(0, dtype=np.int64)
    if not len(cand):
        return []
    starts = np.maximum(cand - max_ed - 1, 0)
    ends = np.minimum(cand + m + max_ed + 1, len(ref))
    hits = []
    s0, e0 = starts[0], ends[0]
    for s1, e1 in list(zip(starts[1:], ends[1:])) + [(None, None)]:
        if s1 is not None and s1 <= e0:
            e0 = max(e0, e1)
            continue
        d = semiglobal_end_distances(pat, ref[s0:e0])
        for j in np.where(d <= max_ed)[0]:
            hits.append((int(s0 + j), int(d[j])))
        if s1 is not None:
            s0, e0 = s1, e1
    return dedup(hits, min_sep)


def gs_env(paths, project):
    if any(p.startswith("gs://") for p in paths):
        if not os.environ.get("GCS_OAUTH_TOKEN"):
            os.environ["GCS_OAUTH_TOKEN"] = subprocess.check_output(["gcloud", "auth", "print-access-token"], text=True).strip()
        rp = project or os.environ.get("GCS_REQUESTER_PAYS_PROJECT") or os.environ.get("GOOGLE_PROJECT")
        if rp:
            os.environ["GCS_REQUESTER_PAYS_PROJECT"] = rp


def fetch_ref(path, chrom, start, end):
    """Sequence of chrom:start-end (1-based, inclusive) from a FASTA holding either whole chromosomes or segments named
    'chrom:s-e' (samtools faidx style) that contain the interval."""
    fa = pysam.FastaFile(path)
    if chrom in fa.references:
        return fa.fetch(chrom, start - 1, end).upper()
    for name in fa.references:
        if name.startswith(chrom + ":"):
            _, s, e = parse_region(name)
            if s <= start and end <= e:
                return fa.fetch(name, start - s, end - s + 1).upper()
    sys.exit(f"{chrom}:{start}-{end} is not in {path} (neither as a chromosome nor inside a chrom:start-end segment)")


def parse_region(r):
    chrom, rng = r.rsplit(":", 1)
    s, e = (int(x.replace(",", "")) for x in rng.split("-"))
    return chrom, s, e


def fetch_read(bam, read_name, region):
    """Primary + supplementary alignments of the read (SA tag followed), the alignment carrying the full sequence, and
    the read sequence in its original orientation."""
    if region:
        chrom, rs, re_ = parse_region(region)
        it = bam.fetch(chrom, rs - 1, re_)
    else:
        print("no --region: scanning the whole BAM for the read", file=sys.stderr)
        it = bam.fetch(until_eof=True)
    alns = [r for r in it if r.query_name == read_name and not r.is_secondary and not r.is_unmapped]
    if not alns:
        sys.exit(f"read {read_name} not found" + (f" in {region}" if region else ""))
    for r in list(alns):  # supplementary alignments elsewhere, from the SA tag (chrom,pos,strand,CIGAR,mapq,NM;)
        if not r.has_tag("SA"):
            continue
        for e in r.get_tag("SA").strip(";").split(";"):
            c, pos = e.split(",")[:2]
            pos = int(pos) - 1
            if any(x.reference_name == c and x.reference_start == pos for x in alns):
                continue
            alns += [x for x in bam.fetch(c, pos, pos + 1)
                     if x.query_name == read_name and x.reference_start == pos and not x.is_secondary]
    alns = list({(r.reference_name, r.reference_start, r.is_supplementary): r for r in alns}.values())
    full = [r for r in alns if r.query_sequence and r.cigartuples and not any(op == 5 for op, _ in r.cigartuples)]
    if not full:
        sys.exit("no alignment of the read carries the full sequence (all hard-clipped)")
    prim = next((r for r in full if not r.is_supplementary), full[0])
    seq = prim.query_sequence if not prim.is_reverse else revcomp(prim.query_sequence)  # original read orientation
    return alns, prim, seq


def read_coords(r, qs, qe):
    """Convert query_sequence coordinates [qs, qe) of an alignment to original read coordinates."""
    hard_left = r.cigartuples[0][1] if r.cigartuples[0][0] == 5 else 0
    q0, q1 = qs + hard_left, qe + hard_left
    if r.is_reverse:
        rlen = r.infer_read_length()
        q0, q1 = rlen - q1, rlen - q0
    return q0, q1


def alignment_blocks(alns):
    """(alignment, strand, read_start, read_end, chrom, ref_start, ref_end), 1-based inclusive."""
    blocks = []
    for r in alns:
        for qs, qe, ts, te in matched_blocks(r):
            q0, q1 = read_coords(r, qs, qe)
            blocks.append(("supplementary" if r.is_supplementary else "primary", "-" if r.is_reverse else "+",
                           q0 + 1, q1, r.reference_name, ts + 1, te))
    return blocks


def target_alignments(alns, prim, region):
    if region:
        chrom, rs, re_ = parse_region(region)
        target = [r for r in alns if r.reference_name == chrom and r.reference_end > rs - 1 and r.reference_start < re_]
    else:
        target = [r for r in alns if r.reference_name == prim.reference_name]
    if not target:
        sys.exit("no alignment of the read overlaps --region")
    return target


def read_segment(target, region, pad, qlen, read_range=None):
    """0-based [start, end) on the read: the part aligned to --region (CIGAR walk), else the target alignments' span,
    +/- pad; or --read-range."""
    if read_range:
        q_lo, q_hi = (int(x.replace(",", "")) for x in read_range.split("-"))
        q_lo -= 1
    else:
        spans = []
        for r in target:
            if region:  # read positions aligned to the region's first and last reference bases
                _, rs, re_ = parse_region(region)
                t0, t1 = max(rs - 1, r.reference_start), min(re_, r.reference_end) - 1
                spans.append(read_coords(r, ref_to_query(r, t0), ref_to_query(r, t1) + 1))
            else:
                spans.append(read_coords(r, r.query_alignment_start, r.query_alignment_end))
            print(f"  {r.reference_name}:{r.reference_start + 1:,}-{r.reference_end:,} ({'-' if r.is_reverse else '+'}) "
                  f"-> read {spans[-1][0] + 1:,}-{spans[-1][1]:,}", file=sys.stderr)
        q_lo, q_hi = min(x[0] for x in spans) - pad, max(x[1] for x in spans) + pad
    return max(0, q_lo), min(qlen, q_hi)


def run_windows(seq, q_lo, q_hi, ref, rstart, window, step, max_ed_frac, min_sep, threads):
    """Hits of every window of seq[q_lo:q_hi] in ref (which starts at 1-based rstart)."""
    max_ed = int(round(max_ed_frac * window))
    k = window // (max_ed + 1)
    if k < 8:
        sys.exit(f"seed length {k} (window // (max edit distance + 1)) is too short; "
                 "use a larger --window or a smaller --max-ed-frac")
    G.update(idx=RefIndex(ref, k), seq=seq, window=window, max_ed=max_ed, min_sep=min_sep, rstart=rstart)
    starts = list(range(q_lo, max(q_lo + 1, q_hi - window + 1), step))
    chunks = [starts[i:i + 200] for i in range(0, len(starts), 200)]
    rows = []
    if threads > 1:
        with mp.get_context("fork").Pool(threads) as pool:
            for n, part in enumerate(pool.imap(scan_windows, chunks), 1):
                rows += part
                print(f"  {min(n * 200, len(starts))}/{len(starts)} windows", file=sys.stderr)
    else:
        for n, ch in enumerate(chunks, 1):
            rows += scan_windows(ch)
            print(f"  {min(n * 200, len(starts))}/{len(starts)} windows", file=sys.stderr)
    return rows, max_ed, len(starts)


def orient_to_ref(rows, blocks, q_lo, q_hi, qlen, max_copies=5):
    """If most single-copy hits (windows with <= max_copies hits) are on the reverse strand, report everything on the
    reverse complement of the read: read position p -> qlen - p + 1 and strands flipped, so the read's main diagonal runs
    forward. Returns rows, blocks, (seg_start, seg_end) 1-based, and whether the read was flipped."""
    n = {}
    for r in rows:
        n[r[0]] = n.get(r[0], 0) + 1
    minus = sum(1 for r in rows if n[r[0]] <= max_copies and r[2] == "-")
    plus = sum(1 for r in rows if n[r[0]] <= max_copies and r[2] == "+")
    if minus <= plus:
        return rows, blocks, (q_lo + 1, q_hi), False
    flip = {"+": "-", "-": "+"}
    rows = [(qlen - r[1] + 1, qlen - r[0] + 1, flip[r[2]]) + tuple(r[3:]) for r in rows]
    blocks = [(b[0], flip[b[1]], qlen - b[3] + 1, qlen - b[2] + 1) + tuple(b[4:]) for b in blocks]
    return sorted(rows), blocks, (qlen - q_hi + 1, qlen - q_lo), True


def write_outputs(out, rows, blocks, info):
    with open(out + ".hits.tsv", "w") as fh:
        fh.write("read_start\tread_end\tstrand\tref_start\tref_end\tedit_distance\tidentity\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    with open(out + ".aln.tsv", "w") as fh:
        fh.write("alignment\tstrand\tread_start\tread_end\tchrom\tref_start\tref_end\n")
        for b in blocks:
            fh.write("\t".join(map(str, b)) + "\n")
    with open(out + ".info.tsv", "w") as fh:
        fh.write("key\tvalue\n")
        for key, v in info:
            fh.write(f"{key}\t{v}\n")
    print(f"{len(rows)} hits -> {out}.hits.tsv", file=sys.stderr)


def common_args(p):
    p.add_argument("--bam", required=True)
    p.add_argument("--read-name", required=True)
    p.add_argument("--region", help="chrom:start-end overlapping the read's alignment (speeds up finding it; without it "
                   "the whole BAM is scanned)")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--window", type=int, default=100)
    p.add_argument("--step", type=int, help="bp between window starts (default: the window size, no overlap)")
    p.add_argument("--max-ed-frac", type=float, default=0.1, help="max edit distance as a fraction of the window")
    p.add_argument("--pad", type=int, default=5000, help="bp added on each side of the region, on the reference and "
                   "on the read")
    p.add_argument("--min-sep", type=int, default=50, help="hits closer than this (bp) to a better hit are dropped")
    p.add_argument("--read-range", help="start-end (1-based) on the read to window, instead of the part aligned to "
                   "--region")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--orient-to-ref", action="store_true", help="if the read matches the reference mostly on the reverse "
                   "strand (e.g. a reverse-strand supplementary alignment), report it reverse-complemented so that its "
                   "diagonal runs forward; read positions are then on the reverse-complemented read")
    p.add_argument("--max-ref", type=int, default=5000000, help="refuse reference segments longer than this")
    p.add_argument("--requester-pays-project")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    common_args(p)
    p.add_argument("--ref", required=True, help="reference FASTA (with .fai), local or gs://")
    p.add_argument("--ref-region", help="use this reference segment instead of --region +/- --pad")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    a.step = a.step or a.window
    gs_env([a.bam, a.ref], a.requester_pays_project)

    bam = pysam.AlignmentFile(a.bam, reference_filename=a.ref if a.bam.endswith(".cram") else None)
    alns, prim, seq = fetch_read(bam, a.read_name, a.region)
    qlen = len(seq)
    blocks = alignment_blocks(alns)
    target = target_alignments(alns, prim, a.region)
    if a.ref_region:
        rchrom, rstart, rend = parse_region(a.ref_region)
    elif a.region:
        rchrom, rs, re_ = parse_region(a.region)
        rstart, rend = max(1, rs - a.pad), re_ + a.pad
    else:
        rchrom = target[0].reference_name
        rstart = max(1, min(r.reference_start for r in target) + 1 - a.pad)
        rend = max(r.reference_end for r in target) + a.pad
    q_lo, q_hi = read_segment(target, a.region, a.pad, qlen, a.read_range)
    print(f"{len(target)} target alignment(s); read segment {q_lo + 1:,}-{q_hi:,} of {qlen:,}", file=sys.stderr)
    if rend - rstart + 1 > a.max_ref:
        sys.exit(f"reference segment {rend - rstart + 1:,} bp is longer than --max-ref {a.max_ref:,}; use --ref-region")
    ref = encode(fetch_ref(a.ref, rchrom, rstart, rend))
    print(f"read {a.read_name}: {qlen:,} bp, {len(alns)} alignments; reference {rchrom}:{rstart:,}-{rend:,} "
          f"({len(ref):,} bp)", file=sys.stderr)

    rows, max_ed, n_win = run_windows(seq, q_lo, q_hi, ref, rstart, a.window, a.step, a.max_ed_frac, a.min_sep, a.threads)
    seg, flipped = (q_lo + 1, q_hi), False
    if a.orient_to_ref:
        rows, blocks, seg, flipped = orient_to_ref(rows, blocks, q_lo, q_hi, qlen)
    write_outputs(a.out, rows, blocks, [
        ("read_name", a.read_name), ("read_length", qlen), ("reference", os.path.basename(a.ref)), ("chrom", rchrom),
        ("ref_start", rstart), ("ref_end", rend), ("window", a.window), ("step", a.step),
        ("read_seg_start", seg[0]), ("read_seg_end", seg[1]), ("max_edit_distance", max_ed),
        ("read_orientation", "reverse_complemented" if flipped else "original"),
        ("primary", f"{prim.reference_name}:{prim.reference_start + 1}-{prim.reference_end} "
                    f"({'-' if prim.is_reverse else '+'})")])


G = {}


def scan_windows(starts):
    out = []
    seq, w, idx = G["seq"], G["window"], G["idx"]
    for w0 in starts:
        win = seq[w0:w0 + w]
        for strand, pat in (("+", win), ("-", revcomp(win))):
            for j, ed in hits_for(encode(pat), idx, G["max_ed"], G["min_sep"]):
                out.append((w0 + 1, w0 + len(win), strand, G["rstart"] + max(0, j - len(win) + 1), G["rstart"] + j, ed,
                            round(1 - ed / len(win), 4)))
    return out


def ref_to_query(r, pos):
    """Query index (into query_sequence, 0-based) aligned to reference position pos (0-based), from the CIGAR; a
    position inside a deletion maps to the next aligned query base."""
    q, t = r.query_alignment_start, r.reference_start
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            if t + n > pos:
                return q + max(0, pos - t)
            q += n
            t += n
        elif op == 1:
            q += n
        elif op in (2, 3):
            if t + n > pos:
                return q
            t += n
    return q


def matched_blocks(r, min_gap=50, min_len=50):
    """Alignment blocks from the CIGAR, split only at indels/skips >= min_gap bp (query 0-based start/end excl., ref
    0-based/excl.); blocks shorter than min_len are dropped."""
    out = []
    q, t = r.query_alignment_start, r.reference_start
    qs, ts = q, t
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            q += n
            t += n
        elif op == 1:
            if n >= min_gap:
                out.append((qs, q, ts, t))
                qs, ts = q + n, t
            q += n
        elif op in (2, 3):
            if n >= min_gap:
                out.append((qs, q, ts, t))
                qs, ts = q, t + n
            t += n
    out.append((qs, q, ts, t))
    return [b for b in out if b[1] - b[0] >= min_len]


if __name__ == "__main__":
    main()
