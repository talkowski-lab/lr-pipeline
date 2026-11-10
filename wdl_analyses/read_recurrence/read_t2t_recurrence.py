#!/usr/bin/env python3
"""Recurrence (dot) plot data for one read/contig aligned to hg38, compared against T2T-CHM13 instead of hg38.

The read and the part of it to analyse come from the hg38 BAM exactly as in read_ref_recurrence.py (read fetched in
--region, supplementary alignments through the SA tag, CIGAR walk to the read segment aligned to --region, +/- --pad).
Then:
1. Lift --region to T2T with a UCSC chain file (hg38 -> CHM13, e.g. hg38-chm13v2.over.chain.gz or
   hg38ToHs1.over.chain.gz): points every --lift-step bp across the region are lifted through the best-scoring chain
   that covers each; the chain carrying most points defines the T2T region (its lifted span), +/- --pad. Points in
   chain gaps (sequence absent from the other assembly) are skipped; the lifted fraction is reported. Inside large SDs
   the lift can be partial or split across chains; check the printed summary, or give --t2t-region directly.
2. Lift the read's hg38 alignment blocks (cut into --lift-step bp pieces) through the chain, for the grey guide lines.
3. Match every --window bp window of the read segment against the T2T segment (edit distance <= --max-ed-frac x window,
   both strands), as in read_ref_recurrence.py.

Outputs (same format as read_ref_recurrence.py, reference coordinates on T2T; plot with plot_read_ref_recurrence.R):
  <out>.hits.tsv, <out>.aln.tsv (lifted blocks), <out>.info.tsv (also hg38 region, lifted fraction)
  <out>.lift.tsv    every lifted region point: hg38_pos chain_id t2t_chrom t2t_pos t2t_strand

Needs read_ref_recurrence.py in the same directory.

Usage:
  read_t2t_recurrence.py --bam sample.hg38.bam --read-name h1tg003138l --region chr16:29324175-29638676 \\
      --chain hg38-chm13v2.over.chain.gz --t2t-ref chm13v2.0.fa.bgz --out h1tg003138l.t2t \\
      [--t2t-region chr16:...] [--hg38-ref hg38.fa (CRAM only)] [--window 100] [--step 100] [--max-ed-frac 0.1]
      [--pad 5000] [--lift-step 1000] [--threads 8]
"""
import argparse
import gzip
import os
import sys
from collections import Counter

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_ref_recurrence as rr  # noqa: E402


class Chain:
    def __init__(self, header):
        f = header.split()
        self.score = float(f[1])
        self.tname, self.tstart, self.tend = f[2], int(f[5]), int(f[6])
        self.qname, self.qsize, self.qstrand = f[7], int(f[8]), f[9]
        self.id = f[12] if len(f) > 12 else f"{f[2]}:{f[5]}"
        self.tb, self.qb, self.sz = [], [], []
        self._t, self._q = self.tstart, int(f[10])

    def add(self, fields):
        size = int(fields[0])
        self.tb.append(self._t)
        self.qb.append(self._q)
        self.sz.append(size)
        if len(fields) == 3:
            self._t += size + int(fields[1])
            self._q += size + int(fields[2])

    def finish(self):
        self.tb, self.qb, self.sz = np.array(self.tb), np.array(self.qb), np.array(self.sz)

    def lift(self, p):
        """0-based target position -> 0-based query position on the query's + strand, or None in a gap."""
        if not self.tstart <= p < self.tend:
            return None
        i = np.searchsorted(self.tb, p, "right") - 1
        if i < 0 or p >= self.tb[i] + self.sz[i]:
            return None
        q = int(self.qb[i] + p - self.tb[i])
        return self.qsize - 1 - q if self.qstrand == "-" else q


def load_chains(path, chrom):
    chains, cur, skip = [], None, True
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith("chain"):
                if cur:
                    cur.finish()
                cur = Chain(line)
                if cur.tname == chrom:
                    chains.append(cur)
                else:
                    cur = None
                    skip = True
                    continue
                skip = False
            elif line.strip() and not line.startswith("#") and not skip:
                cur.add(line.split())
    if cur:
        cur.finish()
    if not chains:
        sys.exit(f"no chain in {path} starts on {chrom} (is it an hg38 -> T2T chain?)")
    return sorted(chains, key=lambda c: -c.score)


def lift_point(chains, p):
    for c in chains:
        q = c.lift(p)
        if q is not None:
            return c, q
    return None, None


def lift_region(chains, region, lift_step, lift_tsv=None):
    """Lift points every lift_step bp across region; the chain carrying most points gives the T2T span.
    Returns (t2t_chrom, start0, end0, lifted fraction); writes every lifted point to lift_tsv."""
    chrom, rs, re_ = rr.parse_region(region)
    pts = sorted(set(range(rs - 1, re_, lift_step)) | {re_ - 1})
    lifted = []
    for pt in pts:
        c, q = lift_point(chains, pt)
        if c is not None:
            lifted.append((pt, c, q))
    if not lifted:
        sys.exit(f"no point of {region} lifts through the chain file; give --t2t-region")
    if lift_tsv:
        with open(lift_tsv, "w") as fh:
            fh.write("hg38_pos\tchain_id\tt2t_chrom\tt2t_pos\tt2t_strand\n")
            for pt, c, q in lifted:
                fh.write(f"{pt + 1}\t{c.id}\t{c.qname}\t{q + 1}\t{c.qstrand}\n")
    use = Counter(c.id for _, c, _ in lifted)
    print(f"lift {region}: {len(lifted)}/{len(pts)} points lifted; by chain: "
          + ", ".join(f"{cid} {n}" for cid, n in use.most_common(5)), file=sys.stderr)
    best = use.most_common(1)[0][0]
    on_best = [(pt, c, q) for pt, c, q in lifted if c.id == best]
    tchrom = on_best[0][1].qname
    t_lo, t_hi = min(q for _, _, q in on_best), max(q for _, _, q in on_best)
    print(f"  main chain {best}: hg38 {on_best[0][0] + 1:,}-{on_best[-1][0] + 1:,} -> T2T {tchrom}:{t_lo + 1:,}-"
          f"{t_hi + 1:,} ({on_best[0][1].qstrand})", file=sys.stderr)
    return tchrom, t_lo, t_hi, round(len(lifted) / len(pts), 3)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    rr.common_args(p)
    p.add_argument("--chain", required=True, help="UCSC chain file hg38 -> T2T-CHM13 (.gz ok)")
    p.add_argument("--t2t-ref", required=True, help="T2T-CHM13 FASTA (bgzipped or plain, with .fai)")
    p.add_argument("--t2t-region", help="T2T segment to compare against, instead of the lifted --region +/- --pad")
    p.add_argument("--hg38-ref", help="hg38 FASTA, only needed to read a CRAM")
    p.add_argument("--lift-step", type=int, default=1000, help="bp between lifted points / size of lifted block pieces")
    if len(sys.argv) == 1:
        p.print_help()
        sys.exit(1)
    a = p.parse_args()
    if not a.region:
        sys.exit("--region (hg38) is required: it is lifted to T2T")
    a.step = a.step or a.window
    rr.gs_env([a.bam, a.t2t_ref], a.requester_pays_project)

    bam = pysam.AlignmentFile(a.bam, reference_filename=a.hg38_ref)
    alns, prim, seq = rr.fetch_read(bam, a.read_name, a.region)
    qlen = len(seq)
    target = rr.target_alignments(alns, prim, a.region)
    q_lo, q_hi = rr.read_segment(target, a.region, a.pad, qlen, a.read_range)
    print(f"{len(target)} target alignment(s); read segment {q_lo + 1:,}-{q_hi:,} of {qlen:,}", file=sys.stderr)

    # 1. Lift the region
    chrom, rs, re_ = rr.parse_region(a.region)
    chains = load_chains(a.chain, chrom)
    tchrom, t_lo, t_hi, frac = lift_region(chains, a.region, a.lift_step, a.out + ".lift.tsv")
    if a.t2t_region:
        tchrom, tstart, tend = rr.parse_region(a.t2t_region)
    else:
        tstart, tend = max(1, t_lo + 1 - a.pad), t_hi + 1 + a.pad
    if tend - tstart + 1 > a.max_ref:
        sys.exit(f"T2T segment {tend - tstart + 1:,} bp is longer than --max-ref {a.max_ref:,}; use --t2t-region")

    # 2. Lift the alignment blocks (pieces of --lift-step bp; kept when both ends lift through one chain)
    blocks = []
    for kind, strand, r0, r1, bchrom, t0, t1 in rr.alignment_blocks(alns):
        if bchrom != chrom:
            continue
        for s in range(t0 - 1, t1, a.lift_step):
            e = min(s + a.lift_step, t1) - 1
            c1, x = lift_point(chains, s)
            c2, y = lift_point(chains, e)
            if c1 is None or c1 is not c2 or abs(abs(y - x) - (e - s)) > 0.2 * (e - s) + 10:
                continue
            # read coordinates of the piece (linear within the block)
            f0, f1 = (s - (t0 - 1)) / (t1 - t0 + 1), (e + 1 - (t0 - 1)) / (t1 - t0 + 1)
            if strand == "+":
                p0, p1 = r0 + f0 * (r1 - r0 + 1), r0 + f1 * (r1 - r0 + 1) - 1
            else:
                p0, p1 = r1 + 1 - f1 * (r1 - r0 + 1), r1 - f0 * (r1 - r0 + 1)
            st = strand if c1.qstrand == "+" else ("-" if strand == "+" else "+")
            blocks.append((kind, st, int(p0), int(p1), c1.qname, min(x, y) + 1, max(x, y) + 1))

    # 3. Match against T2T
    ref = rr.encode(rr.fetch_ref(a.t2t_ref, tchrom, tstart, tend))
    print(f"read {a.read_name}: {qlen:,} bp; T2T reference {tchrom}:{tstart:,}-{tend:,} ({len(ref):,} bp)",
          file=sys.stderr)
    rows, max_ed, n_win = rr.run_windows(seq, q_lo, q_hi, ref, tstart, a.window, a.step, a.max_ed_frac, a.min_sep,
                                         a.threads)
    seg, flipped = (q_lo + 1, q_hi), False
    if a.orient_to_ref:
        rows, blocks, seg, flipped = rr.orient_to_ref(rows, blocks, q_lo, q_hi, qlen)
    rr.write_outputs(a.out, rows, blocks, [
        ("read_name", a.read_name), ("read_length", qlen), ("reference", "T2T-CHM13"), ("chrom", tchrom),
        ("ref_start", tstart), ("ref_end", tend), ("window", a.window), ("step", a.step),
        ("read_seg_start", seg[0]), ("read_seg_end", seg[1]), ("max_edit_distance", max_ed),
        ("read_orientation", "reverse_complemented" if flipped else "original"),
        ("hg38_region", a.region), ("lifted_fraction", frac),
        ("primary", f"hg38 {prim.reference_name}:{prim.reference_start + 1}-{prim.reference_end} "
                    f"({'-' if prim.is_reverse else '+'})")])


if __name__ == "__main__":
    main()
