#!/usr/bin/env python3
"""Cohort QC for CpG methylation BEDs from pb-cpg-tools (PacBio) or modkit (ONT).

All per-sample data are projected onto a fixed reference CpG index (every CG
dinucleotide in the selected reference contigs), so per-site statistics can be
accumulated additively across thousands of samples and compared across platforms.

Subcommands, in workflow order:
    build-index  reference FASTA (+ optional CpG-island BED) -> CpG index (.npz)
    sample       one methylation BED -> dense per-site arrays, QC row, histograms
    accumulate   a batch of sample .npz files -> additive per-site partial sums
    finalize     partial sums -> per-site stats, binned summaries, PCA, QC flags
    plot         summary tables -> per-platform figures
    compare      two platforms' outputs -> merged tables and comparison figures

Per-site variance is reported three ways because depth confounds it:
    raw          SD of beta across samples (sites with cov >= min_coverage)
    excess       observed variance minus mean expected binomial variance
    matched      SD after hypergeometric downsampling of every sample to the same
                 depth (matched_coverage), removing depth differences entirely
"""

import argparse
import gzip
import hashlib
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

CONTEXTS = ["island", "shore", "shelf", "open_sea"]
NO_CONTEXT = 255
SHORE_BP = 2000
SHELF_BP = 4000

BETA_BINS = 50
COV_MAX = 200
NEIGHBOR_COV_EDGES = [0, 15, 20, 30, 50, 10**9]
SD_BIN_WIDTH = 0.005
SD_BINS = 100
MEAN_BINS_COARSE = 10
GRID_BINS = 50
EXCESS_EDGES = np.round(np.arange(-0.05, 0.1505, 0.001), 4)
PLATFORM_FORMATS = {"PacBio": "pb_cpg_tools", "ONT": "modkit_bedmethyl"}
READ_CHUNK = 5_000_000


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def is_autosome(contig):
    return re.fullmatch(r"(chr)?\d+", contig) is not None


def contig_kind(contig):
    if is_autosome(contig):
        return "autosome"
    if re.fullmatch(r"(chr)?X", contig):
        return "X"
    if re.fullmatch(r"(chr)?Y", contig):
        return "Y"
    return "other"


# ---------------------------------------------------------------------------
# CpG index
# ---------------------------------------------------------------------------


class CpgIndex:
    def __init__(self, path):
        d = np.load(path, allow_pickle=False)
        self.contigs = [str(c) for c in d["contigs"]]
        self.offsets = d["offsets"].astype(np.int64)
        self.pos = d["pos"]
        self.context = d["context"]
        self.has_context = bool(d["has_context"])
        self.n = int(self.offsets[-1])
        self.contig_id = np.repeat(
            np.arange(len(self.contigs), dtype=np.int16), np.diff(self.offsets)
        )
        kinds = np.array([contig_kind(c) for c in self.contigs])
        self.kind = kinds[self.contig_id]

    def contig_slice(self, i):
        return slice(self.offsets[i], self.offsets[i + 1])

    def mask(self, kind):
        return self.kind == kind

    def context_labels(self):
        if not self.has_context:
            return np.full(self.n, "NA", dtype=object)
        names = np.array(CONTEXTS + ["NA"], dtype=object)
        codes = np.where(self.context == NO_CONTEXT, len(CONTEXTS), self.context)
        return names[codes]


def read_fai(path):
    fai = {}
    with open(path) as f:
        for line in f:
            name, length, offset, linebases, linewidth = line.split("\t")[:5]
            fai[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return fai


def iter_gzip_contigs(fa_path, contigs):
    """Stream a (b)gzipped FASTA once, yielding (name, seq) for requested contigs."""
    wanted, seen, name, parts = set(contigs), set(), None, []
    with gzip.open(fa_path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if name in wanted:
                    seen.add(name)
                    yield name, np.frombuffer(b"".join(parts).upper(), dtype=np.uint8)
                if seen == wanted:
                    return
                name, parts = line[1:].split()[0].decode(), []
            elif name in wanted:
                parts.append(line.rstrip(b"\r\n"))
    if name in wanted:
        yield name, np.frombuffer(b"".join(parts).upper(), dtype=np.uint8)


def read_contig(fa_path, fai_entry):
    length, offset, linebases, linewidth = fai_entry
    n_lines = (length + linebases - 1) // linebases
    nbytes = length + n_lines * (linewidth - linebases)
    with open(fa_path, "rb") as f:
        f.seek(offset)
        raw = f.read(nbytes)
    seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:length]
    return np.frombuffer(seq.upper(), dtype=np.uint8)


def find_cpgs(seq):
    return np.nonzero((seq[:-1] == ord("C")) & (seq[1:] == ord("G")))[0].astype(
        np.int32
    )


def merge_intervals(starts, ends):
    order = np.argsort(starts)
    starts, ends = starts[order], ends[order]
    merged_s, merged_e = [], []
    for s, e in zip(starts, ends):
        if merged_e and s <= merged_e[-1]:
            merged_e[-1] = max(merged_e[-1], e)
        else:
            merged_s.append(s)
            merged_e.append(e)
    return np.array(merged_s, dtype=np.int64), np.array(merged_e, dtype=np.int64)


def classify_context(pos, starts, ends):
    """Island / shore (<=2kb) / shelf (2-4kb) / open sea, from merged islands."""
    ctx = np.full(len(pos), CONTEXTS.index("open_sea"), dtype=np.uint8)
    if len(starts) == 0:
        return ctx
    i = np.searchsorted(starts, pos, side="right") - 1
    dist = np.full(len(pos), np.iinfo(np.int64).max, dtype=np.int64)
    left = i >= 0
    li = np.clip(i, 0, None)
    d_left = np.where(pos < ends[li], 0, pos - ends[li] + 1)
    dist = np.where(left, d_left, dist)
    right = i + 1 < len(starts)
    ri = np.clip(i + 1, None, len(starts) - 1)
    d_right = starts[ri] - pos
    dist = np.where(right, np.minimum(dist, d_right), dist)
    ctx[dist <= SHELF_BP] = CONTEXTS.index("shelf")
    ctx[dist <= SHORE_BP] = CONTEXTS.index("shore")
    ctx[dist == 0] = CONTEXTS.index("island")
    return ctx


def cmd_build_index(args):
    fai = read_fai(args.ref_fai)
    islands = None
    if args.cpg_islands_bed:
        opener = gzip.open if args.cpg_islands_bed.endswith(".gz") else open
        rows = []
        with opener(args.cpg_islands_bed, "rt") as f:
            for line in f:
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                fields = line.split("\t")
                rows.append((fields[0], int(fields[1]), int(fields[2])))
        islands = pd.DataFrame(rows, columns=["chrom", "start", "end"])

    missing = [c for c in args.contigs if c not in fai]
    if missing:
        raise SystemExit(f"Contigs not found in {args.ref_fai}: {missing}")
    if args.ref_fa.endswith(".gz"):
        seqs = iter_gzip_contigs(args.ref_fa, args.contigs)
    else:
        seqs = ((c, read_contig(args.ref_fa, fai[c])) for c in args.contigs)

    found = {}
    for contig, seq in seqs:
        pos = find_cpgs(seq)
        if islands is not None:
            sub = islands[islands.chrom == contig]
            s, e = merge_intervals(sub.start.to_numpy(), sub.end.to_numpy())
            ctx = classify_context(pos.astype(np.int64), s, e)
        else:
            ctx = np.full(len(pos), NO_CONTEXT, dtype=np.uint8)
        found[contig] = (pos, ctx)
        log(f"{contig}: {len(pos):,} CpGs")

    contigs, offsets, all_pos, all_ctx = [], [0], [], []
    for contig in args.contigs:
        pos, ctx = found[contig]
        contigs.append(contig)
        offsets.append(offsets[-1] + len(pos))
        all_pos.append(pos)
        all_ctx.append(ctx)

    np.savez(
        args.output,
        contigs=np.array(contigs),
        offsets=np.array(offsets, dtype=np.int64),
        pos=np.concatenate(all_pos),
        context=np.concatenate(all_ctx),
        has_context=np.array(islands is not None),
    )
    log(f"Total: {offsets[-1]:,} CpGs on {len(contigs)} contigs")


# ---------------------------------------------------------------------------
# Per-sample parsing
# ---------------------------------------------------------------------------


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def detect_format(path):
    with open_text(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if (
                len(fields) >= 18
                and fields[5] in ("+", "-", ".")
                and len(fields[3]) <= 5
            ):
                return "modkit_bedmethyl"
            if (
                len(fields) >= 8
                and not fields[4].lstrip("-").replace(".", "").isdigit()
            ):
                return "pb_cpg_tools"
            raise SystemExit(
                f"Unrecognized methylation BED format; first record: {line!r}"
            )
    raise SystemExit(f"No records found in {path}")


def iter_records(path, fmt):
    """Yield chunks of (chrom, pos_plus_strand_C, cov, mod_count, beta, strand_mode)."""
    if fmt == "pb_cpg_tools":
        cols = [0, 1, 3, 4, 5, 6]
        names = ["chrom", "start", "mod_score", "type", "cov", "est_mod_count"]
    else:
        cols = [0, 1, 3, 5, 9, 11]
        names = ["chrom", "start", "code", "strand", "cov", "nmod"]
    reader = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        comment="#",
        usecols=cols,
        dtype={0: str, 4: str} if fmt == "pb_cpg_tools" else {0: str, 3: str, 5: str},
        chunksize=READ_CHUNK,
    )
    for chunk in reader:
        chunk.columns = names
        if fmt == "pb_cpg_tools":
            if (chunk["type"] != "Total").any() and (chunk["type"] == "Total").any():
                chunk = chunk[chunk["type"] == "Total"]
            cov = chunk["cov"].to_numpy(dtype=np.int64)
            # mod_score is the model pileup score (raw beta in "model" mode).
            # est_mod_count is the per-read call count, which is not recoverable from
            # mod_score * cov; it backs "counts" mode and the downsampling.
            beta = chunk["mod_score"].to_numpy(dtype=np.float64) / 100.0
            mod = chunk["est_mod_count"].to_numpy(dtype=np.int64)
            pos = chunk["start"].to_numpy(dtype=np.int64)
            yield chunk["chrom"].to_numpy(), pos, cov, mod, beta, "combined"
        else:
            chunk = chunk[chunk["code"] == "m"]
            strand = chunk["strand"].to_numpy()
            pos = chunk["start"].to_numpy(dtype=np.int64)
            minus = strand == "-"
            pos = np.where(minus, pos - 1, pos)
            cov = chunk["cov"].to_numpy(dtype=np.int64)
            mod = chunk["nmod"].to_numpy(dtype=np.int64)
            with np.errstate(divide="ignore", invalid="ignore"):
                beta = np.where(cov > 0, mod / np.maximum(cov, 1), np.nan)
            mode = "stranded" if (minus.any() or (strand == "+").any()) else "combined"
            yield chunk["chrom"].to_numpy(), pos, cov, mod, beta, mode


def load_records(path, fmt, index):
    contig_lookup = {c: i for i, c in enumerate(index.contigs)}
    parts = {k: [] for k in ("cid", "pos", "cov", "mod", "beta")}
    strand_modes = set()
    n_total = 0
    for chrom, pos, cov, mod, beta, mode in iter_records(path, fmt):
        n_total += len(chrom)
        strand_modes.add(mode)
        cid = pd.Series(chrom).map(contig_lookup).to_numpy(dtype=np.float64)
        keep = ~np.isnan(cid)
        parts["cid"].append(cid[keep].astype(np.int64))
        parts["pos"].append(pos[keep])
        parts["cov"].append(cov[keep])
        parts["mod"].append(mod[keep])
        parts["beta"].append(beta[keep])
    rec = {k: np.concatenate(v) if v else np.array([]) for k, v in parts.items()}
    strand_mode = "stranded" if "stranded" in strand_modes else "combined"
    return rec, n_total, strand_mode


def map_to_index(rec, index):
    """Return (global index, matched mask, offset, match fraction)."""
    best = None
    for off in (0, -1, 1):
        gidx = np.zeros(len(rec["pos"]), dtype=np.int64)
        matched = np.zeros(len(rec["pos"]), dtype=bool)
        for c in np.unique(rec["cid"]):
            sel = rec["cid"] == c
            ref = index.pos[index.contig_slice(int(c))]
            q = rec["pos"][sel] + off
            j = np.searchsorted(ref, q)
            jc = np.clip(j, 0, max(len(ref) - 1, 0))
            hit = (
                (j < len(ref)) & (ref[jc] == q) if len(ref) else np.zeros(len(q), bool)
            )
            gidx[sel] = index.offsets[int(c)] + jc
            matched[sel] = hit
        frac = matched.mean() if len(matched) else 0.0
        if best is None or frac > best[3]:
            best = (gidx, matched, off, frac)
        if frac > 0.99:
            break
    return best


def seed_for(sample_id):
    return int.from_bytes(hashlib.sha256(sample_id.encode()).digest()[:8], "little")


def downsample(cov, mod, target, rng):
    """Hypergeometric draw of exactly `target` reads at sites with cov >= target."""
    out = np.full(len(cov), np.nan, dtype=np.float32)
    ok = cov >= target
    if ok.any():
        mod_ok = np.clip(mod[ok], 0, cov[ok])
        k = rng.hypergeometric(mod_ok, cov[ok] - mod_ok, target)
        out[ok] = k / target
    return out


def neighbor_pairs(index, max_dist):
    """Index i such that CpG i and i+1 are on the same contig within max_dist bp."""
    same = index.contig_id[:-1] == index.contig_id[1:]
    close = (index.pos[1:].astype(np.int64) - index.pos[:-1]) <= max_dist
    auto = index.kind[:-1] == "autosome"
    return np.nonzero(same & close & auto)[0]


def beta_fracs(b):
    if len(b) == 0:
        return np.nan, np.nan, np.nan
    return (
        float((b < 0.2).mean()),
        float(((b >= 0.2) & (b <= 0.8)).mean()),
        float((b > 0.8).mean()),
    )


def safe_mean(x):
    return float(np.mean(x)) if len(x) else np.nan


def safe_median(x):
    return float(np.median(x)) if len(x) else np.nan


def cmd_sample(args):
    index = CpgIndex(args.cpg_index)
    fmt = detect_format(args.bed)
    expected = PLATFORM_FORMATS[args.platform]
    if fmt != expected:
        raise SystemExit(
            f"{args.bed} looks like {fmt}, "
            f"but platform {args.platform} expects {expected}"
        )
    # ONT bedMethyl only has read counts; PacBio offers the model score or counts.
    beta_source = args.pacbio_beta_source if fmt == "pb_cpg_tools" else "counts"
    log(f"{args.sample_id}: parsing {fmt}, beta from {beta_source}")
    rec, n_total, strand_mode = load_records(args.bed, fmt, index)
    gidx, matched, offset, match_frac = map_to_index(rec, index)
    log(f"{args.sample_id}: offset {offset}, reference CpG match {match_frac:.4f}")
    if match_frac < args.min_ref_cpg_match_frac:
        raise SystemExit(
            f"Only {match_frac:.3f} of records on index contigs match reference CpGs "
            f"(best offset {offset}); expected >= {args.min_ref_cpg_match_frac}. "
            "Check the reference build and that modkit was run with --cpg."
        )

    g = gidx[matched]
    cov_full = np.bincount(g, weights=rec["cov"][matched], minlength=index.n)
    mod_full = np.bincount(g, weights=rec["mod"][matched], minlength=index.n)
    cov_full = np.rint(cov_full).astype(np.int64)
    mod_full = np.rint(mod_full).astype(np.int64)
    if beta_source == "model":
        bw = np.bincount(
            g,
            weights=np.nan_to_num(rec["beta"][matched]) * rec["cov"][matched],
            minlength=index.n,
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.where(cov_full > 0, bw / np.maximum(cov_full, 1), np.nan)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = np.where(cov_full > 0, mod_full / np.maximum(cov_full, 1), np.nan)
    beta = np.clip(beta, 0, 1).astype(np.float32)
    rng = np.random.default_rng(seed_for(args.sample_id))
    beta_m = downsample(cov_full, mod_full, args.matched_coverage, rng)
    cov16 = np.minimum(cov_full, np.iinfo(np.uint16).max).astype(np.uint16)

    np.savez_compressed(
        f"{args.prefix}.methyl.npz",
        sample_id=np.array(args.sample_id),
        platform=np.array(args.platform),
        n_sites=np.array(index.n),
        cov=cov16,
        beta=beta.astype(np.float16),
        beta_matched=beta_m.astype(np.float16),
    )

    auto = index.mask("autosome")
    good = auto & (cov_full >= args.min_coverage)
    b_good = beta[good]
    bm_auto = beta_m[auto]
    bm_good = bm_auto[~np.isnan(bm_auto)]
    cov_auto = cov_full[auto]
    covered = cov_auto[cov_auto > 0]
    mean_cov_auto = safe_mean(covered)

    pairs = neighbor_pairs(index, args.max_neighbor_distance)
    c1, c2 = cov_full[pairs], cov_full[pairs + 1]
    both = (c1 >= args.min_coverage) & (c2 >= args.min_coverage)
    d = np.abs(beta[pairs][both] - beta[pairs + 1][both])
    pm1, pm2 = beta_m[pairs], beta_m[pairs + 1]
    both_m = ~np.isnan(pm1) & ~np.isnan(pm2)
    dm = np.abs(pm1[both_m] - pm2[both_m])
    r = (
        float(np.corrcoef(beta[pairs][both], beta[pairs + 1][both])[0, 1])
        if both.sum() > 2
        else np.nan
    )

    low, mid, high = beta_fracs(b_good)
    low_m, mid_m, high_m = beta_fracs(bm_good)
    qc = {
        "sample_id": args.sample_id,
        "platform": args.platform,
        "input_format": fmt,
        "beta_source": beta_source,
        "strand_mode": strand_mode,
        "position_offset": offset,
        "n_records": n_total,
        "n_records_on_index_contigs": int(len(gidx)),
        "ref_cpg_match_frac": match_frac,
        "n_cpgs_covered": int((cov_full > 0).sum()),
        "frac_ref_cpgs_cov1": float((cov_auto >= 1).mean()),
        "frac_ref_cpgs_cov5": float((cov_auto >= 5).mean()),
        "frac_ref_cpgs_cov10": float((cov_auto >= 10).mean()),
        "frac_ref_cpgs_cov20": float((cov_auto >= 20).mean()),
        "frac_ref_cpgs_min_cov": float((cov_auto >= args.min_coverage).mean()),
        "mean_cov": mean_cov_auto,
        "median_cov": safe_median(covered),
        "mean_cov_all_ref_cpgs": float(cov_auto.mean()),
        "n_sites_min_cov": int(good.sum()),
        "mean_beta": safe_mean(b_good),
        "median_beta": safe_median(b_good),
        "frac_low": low,
        "frac_intermediate": mid,
        "frac_high": high,
        "n_sites_matched": int(len(bm_good)),
        "mean_beta_matched": safe_mean(bm_good),
        "frac_low_matched": low_m,
        "frac_intermediate_matched": mid_m,
        "frac_high_matched": high_m,
        "neighbor_absdiff": safe_mean(d),
        "neighbor_absdiff_matched": safe_mean(dm),
        "neighbor_pearson_r": r,
        "n_neighbor_pairs": int(both.sum()),
    }
    for kind in ("X", "Y"):
        m = index.mask(kind)
        if m.any() and mean_cov_auto and not np.isnan(mean_cov_auto):
            ck = cov_full[m]
            qc[f"chr{kind}_cov_ratio"] = (
                safe_mean(ck[ck > 0]) / mean_cov_auto if (ck > 0).any() else 0.0
            )
            bk = beta[m & (cov_full >= args.min_coverage)]
            qc[f"chr{kind}_mean_beta"] = safe_mean(bk)
        else:
            qc[f"chr{kind}_cov_ratio"] = np.nan
            qc[f"chr{kind}_mean_beta"] = np.nan
    for i, name in enumerate(CONTEXTS):
        if index.has_context:
            qc[f"mean_beta_{name}"] = safe_mean(beta[good & (index.context == i)])
        else:
            qc[f"mean_beta_{name}"] = np.nan
    pd.DataFrame([qc]).to_csv(f"{args.prefix}.qc.tsv", sep="\t", index=False)

    rows = []
    edges = np.linspace(0, 1, BETA_BINS + 1)
    for name, vals in (("beta", b_good), ("beta_matched", bm_good)):
        counts, _ = np.histogram(vals, bins=edges)
        for lo, hi, c in zip(edges[:-1], edges[1:], counts):
            rows.append((name, lo, hi, int(c), np.nan))
    cov_counts = np.bincount(np.minimum(cov_auto, COV_MAX + 1), minlength=COV_MAX + 2)
    for v, c in enumerate(cov_counts):
        rows.append(("coverage", v, v + 1 if v <= COV_MAX else np.inf, int(c), np.nan))
    cmin = np.minimum(c1, c2)[both]
    for lo, hi in zip(NEIGHBOR_COV_EDGES[:-1], NEIGHBOR_COV_EDGES[1:]):
        lo_eff = max(lo, args.min_coverage)
        if lo_eff >= hi:
            continue
        sel = (cmin >= lo_eff) & (cmin < hi)
        rows.append(
            (
                "neighbor_absdiff_by_cov",
                lo_eff,
                hi if hi < 10**9 else np.inf,
                int(sel.sum()),
                safe_mean(d[sel]),
            )
        )
    hists = pd.DataFrame(rows, columns=["hist", "bin_lo", "bin_hi", "count", "value"])
    hists.insert(0, "platform", args.platform)
    hists.insert(0, "sample_id", args.sample_id)
    hists.to_csv(f"{args.prefix}.hists.tsv", sep="\t", index=False)
    log(f"{args.sample_id}: done")


# ---------------------------------------------------------------------------
# Accumulation and finalization
# ---------------------------------------------------------------------------

SUM_KEYS = ["n", "s", "ss", "sb", "sc", "n_any", "sc_all", "nm", "sm", "ssm"]


def new_sums(n):
    sums = {k: np.zeros(n, dtype=np.float64) for k in SUM_KEYS}
    for k in ("n", "n_any", "nm"):
        sums[k] = np.zeros(n, dtype=np.uint32)
    return sums


def add_sample(sums, cov, beta, beta_m, min_cov):
    cov = cov.astype(np.float64)
    m = cov >= min_cov
    b = beta[m].astype(np.float64)
    c = cov[m]
    sums["n"][m] += 1
    sums["s"][m] += b
    sums["ss"][m] += b * b
    # Unbiased estimate of the binomial sampling variance p(1-p)/c given b.
    sums["sb"][m] += b * (1 - b) / np.maximum(c - 1, 1)
    sums["sc"][m] += c
    sums["n_any"][cov > 0] += 1
    sums["sc_all"] += cov
    mm = ~np.isnan(beta_m)
    bm = beta_m[mm].astype(np.float64)
    sums["nm"][mm] += 1
    sums["sm"][mm] += bm
    sums["ssm"][mm] += bm * bm


def cmd_accumulate(args):
    index = CpgIndex(args.cpg_index)
    sums = new_sums(index.n)
    pca_sites = np.arange(0, index.n, args.pca_site_stride)
    ids, platforms, block = [], [], []
    for path in args.samples:
        d = np.load(path, allow_pickle=False)
        if int(d["n_sites"]) != index.n:
            raise SystemExit(f"{path} was built against a different CpG index")
        cov, beta = d["cov"], d["beta"]
        add_sample(sums, cov, beta, d["beta_matched"], args.min_coverage)
        pb = beta[pca_sites].astype(np.float16)
        pb[cov[pca_sites] < args.min_coverage] = np.nan
        block.append(pb)
        ids.append(str(d["sample_id"]))
        platforms.append(str(d["platform"]))
        log(f"accumulated {ids[-1]}")
    np.savez(
        args.output,
        n_samples=np.array(len(ids)),
        sample_ids=np.array(ids),
        platforms=np.array(platforms),
        pca_sites=pca_sites,
        pca_beta=np.vstack(block),
        **sums,
    )


def per_site_stats(sums):
    n = sums["n"].astype(np.float64)
    nm = sums["nm"].astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = sums["s"] / n
        var = (sums["ss"] - sums["s"] ** 2 / n) / (n - 1)
        var = np.where(n >= 2, np.maximum(var, 0), np.nan)
        binom_var = np.where(n >= 2, sums["sb"] / n, np.nan)
        mean_m = sums["sm"] / nm
        var_m = (sums["ssm"] - sums["sm"] ** 2 / nm) / (nm - 1)
        var_m = np.where(nm >= 2, np.maximum(var_m, 0), np.nan)
        mean_cov = sums["sc"] / n
    return {
        "n": sums["n"],
        "mean_beta": mean,
        "sd_beta": np.sqrt(var),
        "mean_cov": mean_cov,
        "binom_sd": np.sqrt(binom_var),
        "excess_var": var - binom_var,
        "n_matched": sums["nm"],
        "mean_beta_matched": mean_m,
        "sd_beta_matched": np.sqrt(var_m),
    }


def write_per_site(path, index, stats, n_samples, sums):
    keep = np.nonzero(sums["n_any"] > 0)[0]
    contexts = index.context_labels()
    chunk = 2_000_000
    with open(path[:-3], "w") as out:
        out.write(
            "#chrom\tstart\tend\tcontext\tn_samples_covered\tmean_cov_all\tn\t"
            "mean_beta\tsd_beta\tmean_cov\tbinom_sd\texcess_var\tn_matched\t"
            "mean_beta_matched\tsd_beta_matched\n"
        )
    for i in range(0, len(keep), chunk):
        k = keep[i : i + chunk]
        pos = index.pos[k].astype(np.int64)
        df = pd.DataFrame(
            {
                "chrom": np.array(index.contigs)[index.contig_id[k]],
                "start": pos,
                "end": pos + 2,
                "context": contexts[k],
                "n_samples_covered": sums["n_any"][k],
                "mean_cov_all": sums["sc_all"][k] / n_samples,
                **{key: stats[key][k] for key in stats},
            }
        )
        df.to_csv(
            path[:-3],
            sep="\t",
            index=False,
            header=False,
            mode="a",
            float_format="%.5g",
            na_rep="NA",
        )
    subprocess.run(["bgzip", "-f", path[:-3]], check=True)
    subprocess.run(["tabix", "-f", "-p", "bed", path], check=True)


def context_groups(index, base_mask):
    groups = [("all", base_mask)]
    if index.has_context:
        for i, name in enumerate(CONTEXTS):
            groups.append((name, base_mask & (index.context == i)))
    return groups


def mean_bin_labels(mean, nbins):
    b = np.clip((mean * nbins).astype(np.int64), 0, nbins - 1)
    return b


def summarize_sites(index, stats, min_samples, matched_coverage):
    auto = index.mask("autosome")
    sd_edges = np.arange(SD_BINS + 1) * SD_BIN_WIDTH
    grid_edges = np.linspace(0, 1, GRID_BINS + 1)
    sd_grid_edges = np.linspace(0, 0.5, GRID_BINS + 1)
    hist_rows, by_mean_rows, grid_rows = [], [], []
    for metric, nkey, mkey, skey in (
        ("raw", "n", "mean_beta", "sd_beta"),
        ("matched", "n_matched", "mean_beta_matched", "sd_beta_matched"),
    ):
        ok = auto & (stats[nkey] >= min_samples) & ~np.isnan(stats[skey])
        for ctx, mask in context_groups(index, ok):
            mean, sd = stats[mkey][mask], stats[skey][mask]
            binom = stats["binom_sd"][mask]
            excess = stats["excess_var"][mask]
            mb = mean_bin_labels(mean, MEAN_BINS_COARSE)
            for b in range(-1, MEAN_BINS_COARSE):
                sel = np.ones(len(mean), bool) if b < 0 else mb == b
                lo = 0.0 if b < 0 else b / MEAN_BINS_COARSE
                hi = 1.0 if b < 0 else (b + 1) / MEAN_BINS_COARSE
                label = "all" if b < 0 else f"{lo:.1f}-{hi:.1f}"
                s = sd[sel]
                counts, _ = np.histogram(
                    np.minimum(s, sd_edges[-1] - 1e-9), bins=sd_edges
                )
                for elo, ehi, c in zip(sd_edges[:-1], sd_edges[1:], counts):
                    hist_rows.append((metric, ctx, label, lo, hi, elo, ehi, int(c)))
                q = (
                    np.quantile(s, [0.05, 0.25, 0.5, 0.75, 0.95])
                    if len(s)
                    else [np.nan] * 5
                )
                by_mean_rows.append(
                    (
                        metric,
                        ctx,
                        label,
                        lo,
                        hi,
                        int(len(s)),
                        *q,
                        safe_mean(s),
                        safe_mean(binom[sel]) if metric == "raw" else np.nan,
                        safe_mean(excess[sel]) if metric == "raw" else np.nan,
                        float((s > 0.1).mean()) if len(s) else np.nan,
                        float((s > 0.2).mean()) if len(s) else np.nan,
                        matched_coverage,
                    )
                )
            if ctx == "all":
                h2, _, _ = np.histogram2d(
                    np.clip(mean, 0, 1 - 1e-9),
                    np.minimum(sd, 0.5 - 1e-9),
                    bins=[grid_edges, sd_grid_edges],
                )
                for i in range(GRID_BINS):
                    for j in range(GRID_BINS):
                        if h2[i, j]:
                            grid_rows.append(
                                (
                                    metric,
                                    grid_edges[i],
                                    grid_edges[i + 1],
                                    sd_grid_edges[j],
                                    sd_grid_edges[j + 1],
                                    int(h2[i, j]),
                                )
                            )
    sd_hist = pd.DataFrame(
        hist_rows,
        columns=[
            "metric",
            "context",
            "mean_bin",
            "mean_lo",
            "mean_hi",
            "sd_lo",
            "sd_hi",
            "count",
        ],
    )
    by_mean = pd.DataFrame(
        by_mean_rows,
        columns=[
            "metric",
            "context",
            "mean_bin",
            "mean_lo",
            "mean_hi",
            "n_sites",
            "sd_q05",
            "sd_q25",
            "sd_median",
            "sd_q75",
            "sd_q95",
            "sd_mean",
            "binom_sd_mean",
            "excess_var_mean",
            "frac_sd_gt_0.1",
            "frac_sd_gt_0.2",
            "matched_coverage",
        ],
    )
    grid = pd.DataFrame(
        grid_rows, columns=["metric", "mean_lo", "mean_hi", "sd_lo", "sd_hi", "count"]
    )

    ok = auto & (stats["n"] >= min_samples) & ~np.isnan(stats["excess_var"])
    ex_rows = []
    for ctx, mask in context_groups(index, ok):
        ex = np.clip(
            stats["excess_var"][mask], EXCESS_EDGES[0], EXCESS_EDGES[-1] - 1e-9
        )
        counts, _ = np.histogram(ex, bins=EXCESS_EDGES)
        for lo, hi, c in zip(EXCESS_EDGES[:-1], EXCESS_EDGES[1:], counts):
            ex_rows.append((ctx, lo, hi, int(c)))
    excess = pd.DataFrame(ex_rows, columns=["context", "bin_lo", "bin_hi", "count"])
    return sd_hist, by_mean, grid, excess


def run_pca(beta, min_frac_present=0.8, n_pcs=10):
    x = beta.astype(np.float32)
    present = ~np.isnan(x)
    keep = present.mean(axis=0) >= min_frac_present
    x, present = x[:, keep], present[:, keep]
    if x.shape[0] < 3 or x.shape[1] < 3:
        return None, None, int(keep.sum())
    col_mean = np.nanmean(x, axis=0)
    x = np.where(present, x, col_mean) - col_mean
    gram = x @ x.T
    vals, vecs = np.linalg.eigh(gram.astype(np.float64))
    order = np.argsort(vals)[::-1]
    vals, vecs = np.clip(vals[order], 0, None), vecs[:, order]
    k = min(n_pcs, x.shape[0] - 1)
    coords = vecs[:, :k] * np.sqrt(vals[:k])
    var_frac = vals[:k] / vals.sum() if vals.sum() > 0 else np.zeros(k)
    return coords, var_frac, int(keep.sum())


def write_pca(prefix, ids, platforms, coords, var_frac, n_sites):
    if coords is None:
        pd.DataFrame(columns=["sample_id", "platform"]).to_csv(
            f"{prefix}.pca_coords.tsv", sep="\t", index=False
        )
        pd.DataFrame(columns=["pc", "variance_fraction", "n_sites"]).to_csv(
            f"{prefix}.pca_variance.tsv", sep="\t", index=False
        )
        return
    df = pd.DataFrame(coords, columns=[f"PC{i + 1}" for i in range(coords.shape[1])])
    df.insert(0, "platform", platforms)
    df.insert(0, "sample_id", ids)
    df.to_csv(f"{prefix}.pca_coords.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(
        {
            "pc": [f"PC{i + 1}" for i in range(len(var_frac))],
            "variance_fraction": var_frac,
            "n_sites": n_sites,
        }
    ).to_csv(f"{prefix}.pca_variance.tsv", sep="\t", index=False)


OUTLIER_METRICS = [
    "mean_cov",
    "frac_ref_cpgs_min_cov",
    "mean_beta",
    "frac_intermediate",
    "frac_intermediate_matched",
    "neighbor_absdiff",
    "neighbor_absdiff_matched",
    "ref_cpg_match_frac",
]


def flag_outliers(qc, threshold):
    flagged = []
    for m in OUTLIER_METRICS:
        x = qc[m].astype(float)
        med = np.nanmedian(x)
        mad = 1.4826 * np.nanmedian(np.abs(x - med))
        z = (x - med) / mad if mad > 0 else pd.Series(0.0, index=x.index)
        qc[f"{m}_robust_z"] = z
        flagged.append(np.where(np.abs(z) > threshold, m, ""))
    flagged = np.array(flagged).T
    qc["outlier_metrics"] = [",".join(f for f in row if f) for row in flagged]
    qc["is_outlier"] = qc["outlier_metrics"] != ""
    return qc


def cmd_finalize(args):
    index = CpgIndex(args.cpg_index)
    sums = None
    ids, platforms, pca_blocks, pca_sites = [], [], [], None
    for path in args.partials:
        d = np.load(path, allow_pickle=False)
        if sums is None:
            sums = {k: d[k].copy() for k in SUM_KEYS}
            pca_sites = d["pca_sites"]
        else:
            for k in SUM_KEYS:
                sums[k] += d[k]
        ids.extend(str(x) for x in d["sample_ids"])
        platforms.extend(str(x) for x in d["platforms"])
        pca_blocks.append(d["pca_beta"])
        log(f"merged {path}")
    n_samples = len(ids)
    if len(set(ids)) != n_samples:
        raise SystemExit("Duplicate sample IDs across inputs")

    stats = per_site_stats(sums)
    write_per_site(
        f"{args.prefix}.per_site_stats.tsv.gz", index, stats, n_samples, sums
    )
    sd_hist, by_mean, grid, excess = summarize_sites(
        index, stats, args.min_samples_per_site, args.matched_coverage
    )
    for df, name in (
        (sd_hist, "site_sd_hist"),
        (by_mean, "site_sd_by_mean"),
        (grid, "site_mean_sd_grid"),
        (excess, "site_excess_var_hist"),
    ):
        df.insert(0, "platform", args.platform)
        df.to_csv(
            f"{args.prefix}.{name}.tsv", sep="\t", index=False, float_format="%.6g"
        )

    pca_beta = np.vstack(pca_blocks)
    np.savez_compressed(
        f"{args.prefix}.pca_matrix.npz",
        sample_ids=np.array(ids),
        platforms=np.array(platforms),
        chrom=np.array(index.contigs)[index.contig_id[pca_sites]],
        pos=index.pos[pca_sites],
        beta=pca_beta,
    )
    auto_cols = index.kind[pca_sites] == "autosome"
    coords, var_frac, n_used = run_pca(pca_beta[:, auto_cols])
    write_pca(args.prefix, ids, platforms, coords, var_frac, n_used)

    qc = pd.read_csv(args.sample_qc, sep="\t")
    flag_outliers(qc, args.outlier_mad).to_csv(
        f"{args.prefix}.sample_qc.tsv", sep="\t", index=False, float_format="%.6g"
    )

    cohort = {
        "platform": args.platform,
        "beta_source": ",".join(sorted(qc["beta_source"].astype(str).unique())),
        "n_samples": n_samples,
        "n_ref_cpgs": index.n,
        "n_sites_any_coverage": int((sums["n_any"] > 0).sum()),
        "min_coverage": args.min_coverage,
        "matched_coverage": args.matched_coverage,
        "min_samples_per_site": args.min_samples_per_site,
        "n_outlier_samples": int(qc["is_outlier"].sum()),
    }
    raw_all = by_mean[
        (by_mean.metric == "raw")
        & (by_mean.context == "all")
        & (by_mean.mean_bin == "all")
    ]
    m_all = by_mean[
        (by_mean.metric == "matched")
        & (by_mean.context == "all")
        & (by_mean.mean_bin == "all")
    ]
    if len(raw_all):
        cohort["n_autosomal_sites_raw"] = int(raw_all.n_sites.iloc[0])
        cohort["median_site_sd_raw"] = float(raw_all.sd_median.iloc[0])
        cohort["mean_binom_sd"] = float(raw_all.binom_sd_mean.iloc[0])
        cohort["mean_excess_var"] = float(raw_all.excess_var_mean.iloc[0])
    if len(m_all):
        cohort["n_autosomal_sites_matched"] = int(m_all.n_sites.iloc[0])
        cohort["median_site_sd_matched"] = float(m_all.sd_median.iloc[0])
    pd.DataFrame([cohort]).to_csv(
        f"{args.prefix}.cohort_summary.tsv", sep="\t", index=False, float_format="%.6g"
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

# Categorical slots, fixed order (blue, orange, aqua, yellow, magenta).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
PLATFORM_COLORS = {"PacBio": SERIES[0], "ONT": SERIES[1]}
TEXT_SECONDARY = "#52514e"
GRID_COLOR = "#e4e3df"
OUTLIER_COLOR = "#0b0b0b"


def setup_matplotlib():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 150,
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": TEXT_SECONDARY,
            "axes.labelcolor": "#0b0b0b",
            "axes.grid": True,
            "grid.color": GRID_COLOR,
            "grid.linewidth": 0.6,
            "xtick.color": TEXT_SECONDARY,
            "ytick.color": TEXT_SECONDARY,
            "lines.linewidth": 2,
            "legend.frameon": False,
        }
    )
    return plt


def sequential_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "seq_blue", ["#f3f7fd", "#2a78d6", "#0c2f5c"]
    )


def savefig(fig, path_base, written):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{path_base}.{ext}")
        written.append(f"{path_base}.{ext}")


def context_colors(contexts):
    order = ["all"] + CONTEXTS
    return {
        c: SERIES[order.index(c)] if c in order else TEXT_SECONDARY for c in contexts
    }


def plot_sample_qc(plt, qc, out, written, color_by="platform"):
    metrics = [
        ("mean_cov", "Mean coverage (covered CpGs)"),
        ("frac_ref_cpgs_min_cov", "Frac. ref CpGs >= min cov"),
        ("mean_beta", "Mean beta"),
        ("frac_intermediate", "Frac. intermediate (0.2-0.8)"),
        ("frac_intermediate_matched", "Frac. intermediate, matched depth"),
        ("neighbor_absdiff", "Neighbour |delta beta|"),
        ("neighbor_absdiff_matched", "Neighbour |delta beta|, matched"),
        ("chrX_cov_ratio", "chrX / autosome coverage"),
    ]
    groups = list(dict.fromkeys(qc[color_by]))
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.5))
    rng = np.random.default_rng(0)
    for ax, (m, title) in zip(axes.flat, metrics):
        if m not in qc or qc[m].isna().all():
            ax.set_visible(False)
            continue
        for gi, g in enumerate(groups):
            sub = qc[qc[color_by] == g]
            vals = sub[m].dropna().to_numpy()
            if len(vals) > 1:
                vp = ax.violinplot(vals, positions=[gi], showextrema=False, widths=0.8)
                for body in vp["bodies"]:
                    body.set_facecolor(PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)]))
                    body.set_alpha(0.35)
            x = gi + rng.uniform(-0.15, 0.15, len(sub))
            out_mask = sub.get(
                "is_outlier", pd.Series(False, index=sub.index)
            ).to_numpy(bool)
            col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
            ax.scatter(
                x[~out_mask],
                sub[m].to_numpy()[~out_mask],
                s=6,
                color=col,
                alpha=0.6,
                linewidths=0,
            )
            if out_mask.any():
                ax.scatter(
                    x[out_mask],
                    sub[m].to_numpy()[out_mask],
                    s=22,
                    marker="x",
                    color=OUTLIER_COLOR,
                    linewidths=1,
                    label="outlier (>4 MAD)",
                )
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([f"{g}\nn={(qc[color_by] == g).sum()}" for g in groups])
        ax.set_title(title, fontsize=9)
    handles, labels = [], []
    for ax in axes.flat:
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    if handles:
        fig.legend(handles, labels, loc="upper right")
    fig.suptitle("Per-sample QC metrics")
    savefig(fig, out, written)
    plt.close(fig)


def hist_matrix(hists, name):
    h = hists[hists["hist"] == name]
    piv = h.pivot_table(
        index="sample_id", columns="bin_lo", values="count", aggfunc="sum"
    )
    return piv


def plot_beta_distributions(plt, hists, out, written, group_col="platform"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, name, title in zip(
        axes, ("beta", "beta_matched"), ("Raw (cov >= min coverage)", "Matched depth")
    ):
        sub = hists[hists["hist"] == name]
        for gi, (g, gdf) in enumerate(sub.groupby(group_col, sort=False)):
            piv = hist_matrix(gdf, name)
            if piv.empty:
                continue
            dens = piv.div(piv.sum(axis=1).replace(0, np.nan), axis=0)
            centers = np.array(piv.columns, dtype=float) + 0.5 / BETA_BINS
            col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
            show = dens.sample(min(200, len(dens)), random_state=0)
            for _, row in show.iterrows():
                ax.plot(centers, row.to_numpy(), color=col, alpha=0.08, linewidth=0.7)
            ax.plot(
                centers,
                dens.median(axis=0).to_numpy(),
                color=col,
                linewidth=2,
                label=f"{g} median (n={len(dens)})",
            )
        ax.set_title(title)
        ax.set_xlabel("Methylation beta")
        ax.legend()
    axes[0].set_ylabel("Fraction of CpGs")
    fig.suptitle("Per-sample methylation distributions (thin lines: up to 200 samples)")
    savefig(fig, out, written)
    plt.close(fig)


def plot_coverage(plt, hists, qc, out, written, group_col="platform"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    sub = hists[(hists["hist"] == "coverage") & (hists["bin_lo"] <= COV_MAX)]
    for gi, (g, gdf) in enumerate(sub.groupby(group_col, sort=False)):
        piv = hist_matrix(gdf, "coverage")
        dens = piv.div(piv.sum(axis=1).replace(0, np.nan), axis=0)
        x = np.array(piv.columns, dtype=float)
        col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
        q25, q50, q75 = (dens.quantile(q, axis=0).to_numpy() for q in (0.25, 0.5, 0.75))
        axes[0].fill_between(x, q25, q75, color=col, alpha=0.25, linewidth=0)
        axes[0].plot(x, q50, color=col, label=f"{g} median, IQR band")
    axes[0].set_xlabel("Coverage at reference CpG (autosomes)")
    axes[0].set_ylabel("Fraction of reference CpGs")
    axes[0].set_xlim(0, 100)
    axes[0].legend()
    for gi, (g, gdf) in enumerate(qc.groupby(group_col, sort=False)):
        col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
        axes[1].hist(gdf["mean_cov"].dropna(), bins=40, color=col, alpha=0.6, label=g)
    axes[1].set_xlabel("Per-sample mean coverage (covered CpGs)")
    axes[1].set_ylabel("Samples")
    axes[1].legend()
    fig.suptitle("Coverage")
    savefig(fig, out, written)
    plt.close(fig)


def sd_density(sd_hist, metric, context, mean_bin="all"):
    h = sd_hist[
        (sd_hist.metric == metric)
        & (sd_hist.context == context)
        & (sd_hist.mean_bin == mean_bin)
    ]
    total = h["count"].sum()
    centers = (h.sd_lo + h.sd_hi).to_numpy() / 2
    if total == 0:
        return centers, np.zeros(len(h)), np.zeros(len(h))
    frac = h["count"].to_numpy() / total
    return centers, frac / SD_BIN_WIDTH, np.cumsum(frac)


def plot_site_sd(plt, sd_hist, out, written):
    contexts = list(dict.fromkeys(sd_hist.context))
    colors = context_colors(contexts)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, metric in zip(axes, ("raw", "matched")):
        for c in contexts:
            x, dens, _ = sd_density(sd_hist, metric, c)
            ax.plot(x, dens, color=colors[c], label=c)
        ax.set_title(f"{metric} beta")
        ax.set_xlabel("Per-site SD of beta across samples")
        ax.set_xlim(0, 0.35)
        ax.legend(title="CpG context")
    axes[0].set_ylabel("Density (autosomal sites)")
    fig.suptitle("Per-CpG methylation variation across samples")
    savefig(fig, out, written)
    plt.close(fig)


def plot_mean_vs_sd(plt, grid, by_mean, out, written):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    cmap = sequential_cmap()
    from matplotlib.colors import LogNorm

    for ax, metric in zip(axes, ("raw", "matched")):
        g = grid[grid.metric == metric]
        if g.empty:
            ax.set_visible(False)
            continue
        mat = np.zeros((GRID_BINS, GRID_BINS))
        i = np.clip(
            np.rint(g.mean_lo.to_numpy() * GRID_BINS).astype(int), 0, GRID_BINS - 1
        )
        j = np.clip(
            np.rint(g.sd_lo.to_numpy() * 2 * GRID_BINS).astype(int), 0, GRID_BINS - 1
        )
        np.add.at(mat, (i, j), g["count"].to_numpy())
        mat = np.ma.masked_equal(mat, 0)
        im = ax.imshow(
            mat.T,
            origin="lower",
            extent=[0, 1, 0, 0.5],
            aspect="auto",
            cmap=cmap,
            norm=LogNorm(),
        )
        fig.colorbar(im, ax=ax, label="CpG sites")
        bm = by_mean[
            (by_mean.metric == metric)
            & (by_mean.context == "all")
            & (by_mean.mean_bin != "all")
        ]
        mid = (bm.mean_lo + bm.mean_hi) / 2
        ax.plot(
            mid,
            bm.sd_median,
            color=SERIES[1],
            marker="o",
            markersize=4,
            label="median SD",
        )
        p = np.linspace(0.005, 0.995, 200)
        if metric == "raw":
            ax.plot(
                mid,
                bm.binom_sd_mean,
                color="#0b0b0b",
                linestyle="--",
                linewidth=1.5,
                label="expected binomial SD (mean)",
            )
        else:
            mc = float(bm.matched_coverage.iloc[0]) if len(bm) else np.nan
            ax.plot(
                p,
                np.sqrt(p * (1 - p) / mc),
                color="#0b0b0b",
                linestyle="--",
                linewidth=1.5,
                label=f"binomial SD at {mc:g}x",
            )
        ax.set_title(f"{metric} beta")
        ax.set_xlabel("Per-site mean beta")
        ax.legend(loc="upper left")
    axes[0].set_ylabel("Per-site SD of beta")
    fig.suptitle("Mean-variance relationship (autosomal CpGs)")
    savefig(fig, out, written)
    plt.close(fig)


def plot_excess(plt, excess, out, written):
    contexts = list(dict.fromkeys(excess.context))
    colors = context_colors(contexts)
    fig, ax = plt.subplots(figsize=(7, 4))
    for c in contexts:
        h = excess[excess.context == c]
        tot = h["count"].sum()
        if tot:
            ax.plot(
                (h.bin_lo + h.bin_hi) / 2, h["count"] / tot, color=colors[c], label=c
            )
    ax.axvline(0, color=TEXT_SECONDARY, linewidth=1)
    ax.set_xlabel("Excess variance (observed - expected binomial)")
    ax.set_ylabel("Fraction of sites")
    ax.legend(title="CpG context")
    ax.set_title("Variance beyond read-sampling noise (clipped to plot range)")
    savefig(fig, out, written)
    plt.close(fig)


def plot_neighbor(plt, hists, qc, out, written, group_col="platform"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    nb = hists[hists["hist"] == "neighbor_absdiff_by_cov"].dropna(subset=["value"])
    for gi, (g, gdf) in enumerate(nb.groupby(group_col, sort=False)):
        col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
        agg = gdf.groupby("bin_lo")["value"].quantile([0.25, 0.5, 0.75]).unstack()
        x = np.arange(len(agg)) + (gi - 0.5) * 0.15
        axes[0].errorbar(
            x,
            agg[0.5],
            yerr=[agg[0.5] - agg[0.25], agg[0.75] - agg[0.5]],
            fmt="o",
            color=col,
            label=f"{g} median, IQR",
            markersize=5,
        )
        axes[0].set_xticks(np.arange(len(agg)))
        axes[0].set_xticklabels([f">={int(v)}x" for v in agg.index])
    axes[0].set_xlabel("Min coverage of the CpG pair")
    axes[0].set_ylabel("Mean |delta beta|, adjacent CpGs")
    axes[0].legend()
    for gi, (g, gdf) in enumerate(qc.groupby(group_col, sort=False)):
        col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
        axes[1].scatter(
            gdf.mean_cov,
            gdf.neighbor_absdiff,
            s=8,
            color=col,
            alpha=0.6,
            label=f"{g} raw",
            linewidths=0,
        )
        axes[1].scatter(
            gdf.mean_cov,
            gdf.neighbor_absdiff_matched,
            s=10,
            marker="^",
            facecolors="none",
            edgecolors=col,
            alpha=0.6,
            label=f"{g} matched depth",
        )
    axes[1].set_xlabel("Per-sample mean coverage")
    axes[1].set_ylabel("Mean |delta beta|, adjacent CpGs")
    axes[1].legend()
    fig.suptitle("Within-sample neighbour discordance (lower = smoother/cleaner)")
    savefig(fig, out, written)
    plt.close(fig)


def plot_intermediate_vs_cov(plt, qc, out, written, group_col="platform"):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, m, title in zip(
        axes,
        ("frac_intermediate", "frac_intermediate_matched"),
        ("Raw", "Matched depth"),
    ):
        for gi, (g, gdf) in enumerate(qc.groupby(group_col, sort=False)):
            col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
            ax.scatter(
                gdf.mean_cov, gdf[m], s=8, color=col, alpha=0.6, label=g, linewidths=0
            )
        ax.set_title(title)
        ax.set_xlabel("Per-sample mean coverage")
        ax.legend()
    axes[0].set_ylabel("Fraction of CpGs with 0.2 <= beta <= 0.8")
    fig.suptitle("Intermediate methylation vs depth")
    savefig(fig, out, written)
    plt.close(fig)


def plot_pca(plt, coords, var, qc, out, written, color_by="platform"):
    if coords is None or coords.empty or "PC2" not in coords:
        return
    df = coords.merge(
        qc[["sample_id", "mean_cov", "is_outlier"]], on="sample_id", how="left"
    )
    pairs = [("PC1", "PC2")] + ([("PC3", "PC4")] if "PC4" in df else [])
    fig, axes = plt.subplots(1, len(pairs) + 1, figsize=(5 * (len(pairs) + 1), 4.3))
    vf = dict(zip(var.pc, var.variance_fraction)) if len(var) else {}
    for ax, (a, b) in zip(axes, pairs):
        for gi, (g, gdf) in enumerate(df.groupby(color_by, sort=False)):
            col = PLATFORM_COLORS.get(g, SERIES[gi % len(SERIES)])
            ax.scatter(
                gdf[a], gdf[b], s=10, color=col, alpha=0.7, label=g, linewidths=0
            )
        out_df = df[df.is_outlier.fillna(False).astype(bool)]
        if len(out_df):
            ax.scatter(
                out_df[a],
                out_df[b],
                s=30,
                marker="x",
                color=OUTLIER_COLOR,
                linewidths=1,
                label="QC outlier",
            )
        ax.set_xlabel(f"{a} ({100 * vf.get(a, np.nan):.1f}%)")
        ax.set_ylabel(f"{b} ({100 * vf.get(b, np.nan):.1f}%)")
        ax.legend()
    sc = axes[-1].scatter(
        df.PC1, df.PC2, c=df.mean_cov, s=10, cmap=sequential_cmap(), linewidths=0
    )
    fig.colorbar(sc, ax=axes[-1], label="Mean coverage")
    axes[-1].set_xlabel("PC1")
    axes[-1].set_ylabel("PC2")
    axes[-1].set_title("Coloured by depth")
    fig.suptitle("PCA of samples (autosomal CpG subset, mean-imputed)")
    savefig(fig, out, written)
    plt.close(fig)


def read_tsv(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def cmd_plot(args):
    plt = setup_matplotlib()
    qc = read_tsv(args.sample_qc)
    hists = read_tsv(args.sample_hists)
    sd_hist = read_tsv(args.site_sd_hist)
    by_mean = read_tsv(args.site_sd_by_mean)
    grid = read_tsv(args.site_mean_sd_grid)
    excess = read_tsv(args.site_excess_var_hist)
    coords = read_tsv(args.pca_coords)
    var = read_tsv(args.pca_variance)
    os.makedirs(args.outdir, exist_ok=True)
    p = os.path.join(args.outdir, args.prefix)
    written = []
    plot_sample_qc(plt, qc, f"{p}.sample_qc_metrics", written)
    plot_beta_distributions(plt, hists, f"{p}.beta_distributions", written)
    plot_coverage(plt, hists, qc, f"{p}.coverage", written)
    plot_site_sd(plt, sd_hist, f"{p}.site_sd_distribution", written)
    plot_mean_vs_sd(plt, grid, by_mean, f"{p}.site_mean_vs_sd", written)
    plot_excess(plt, excess, f"{p}.site_excess_variance", written)
    plot_neighbor(plt, hists, qc, f"{p}.neighbor_discordance", written)
    plot_intermediate_vs_cov(plt, qc, f"{p}.intermediate_vs_coverage", written)
    plot_pca(plt, coords, var, qc, f"{p}.pca", written)
    log("\n".join(written))


# ---------------------------------------------------------------------------
# Cross-platform comparison
# ---------------------------------------------------------------------------


def concat_platforms(paths, labels):
    frames = []
    for path, label in zip(paths, labels):
        df = read_tsv(path)
        df["platform"] = label
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_site_table(path, label, min_samples):
    cols = [
        "#chrom",
        "start",
        "context",
        "n",
        "mean_beta",
        "sd_beta",
        "n_matched",
        "mean_beta_matched",
        "sd_beta_matched",
    ]
    frames = []
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=cols,
        chunksize=4_000_000,
        dtype={"#chrom": "category", "context": "category"},
        na_values="NA",
    ):
        keep = chunk["n"] >= min_samples
        frames.append(chunk[keep])
    df = pd.concat(frames, ignore_index=True).rename(columns={"#chrom": "chrom"})
    for c in ("mean_beta", "sd_beta", "mean_beta_matched", "sd_beta_matched"):
        df[c] = df[c].astype(np.float32)
    df = df.drop(columns=["n"]).rename(
        columns={
            c: f"{c}_{label}"
            for c in (
                "mean_beta",
                "sd_beta",
                "n_matched",
                "mean_beta_matched",
                "sd_beta_matched",
            )
        }
    )
    return df


def cmd_compare(args):
    plt = setup_matplotlib()
    labels = [args.label_a, args.label_b]
    p = os.path.join(args.outdir, args.prefix)
    os.makedirs(args.outdir, exist_ok=True)
    written = []

    merged = {}
    for name, paths in (
        ("sample_qc", args.sample_qc),
        ("sample_hists", args.sample_hists),
        ("site_sd_hist", args.site_sd_hist),
        ("site_sd_by_mean", args.site_sd_by_mean),
        ("site_mean_sd_grid", args.site_mean_sd_grid),
        ("site_excess_var_hist", args.site_excess_var_hist),
        ("cohort_summary", args.cohort_summary),
    ):
        merged[name] = concat_platforms(paths, labels)
        merged[name].to_csv(
            f"{p}.merged_{name}.tsv", sep="\t", index=False, float_format="%.6g"
        )
    mc = merged["cohort_summary"]["matched_coverage"].unique()
    if len(mc) > 1:
        raise SystemExit(f"Runs used different matched_coverage values: {mc}")

    # Joint per-site comparison on sites both platforms cover.
    a = load_site_table(args.per_site_stats[0], labels[0], args.min_samples_per_site)
    b = load_site_table(args.per_site_stats[1], labels[1], args.min_samples_per_site)
    chroms = sorted(set(a.chrom.cat.categories) | set(b.chrom.cat.categories))
    a["chrom"] = a["chrom"].cat.set_categories(chroms)
    b["chrom"] = b["chrom"].cat.set_categories(chroms)
    joint = a.merge(b.drop(columns=["context"]), on=["chrom", "start"], how="inner")
    del a, b
    joint = joint[joint.chrom.isin([c for c in chroms if is_autosome(c)])]
    for m in ("sd_beta", "sd_beta_matched", "mean_beta"):
        joint[f"delta_{m}"] = joint[f"{m}_{labels[0]}"] - joint[f"{m}_{labels[1]}"]
    joint.insert(2, "end", joint["start"] + 2)
    joint.to_csv(
        f"{p}.joint_site_stats.tsv.gz",
        sep="\t",
        index=False,
        float_format="%.5g",
        na_rep="NA",
        compression="gzip",
    )

    edges = np.linspace(0, 0.5, GRID_BINS + 1)
    grid_rows, delta_rows = [], []
    contexts = ["all"] + [c for c in CONTEXTS if c in set(joint.context.astype(str))]
    for metric, col in (("raw", "sd_beta"), ("matched", "sd_beta_matched")):
        xa, xb = joint[f"{col}_{labels[0]}"], joint[f"{col}_{labels[1]}"]
        ok = xa.notna() & xb.notna()
        h2, _, _ = np.histogram2d(
            np.minimum(xa[ok], 0.5 - 1e-9),
            np.minimum(xb[ok], 0.5 - 1e-9),
            bins=[edges, edges],
        )
        for i in range(GRID_BINS):
            for j in range(GRID_BINS):
                if h2[i, j]:
                    grid_rows.append(
                        (
                            metric,
                            edges[i],
                            edges[i + 1],
                            edges[j],
                            edges[j + 1],
                            int(h2[i, j]),
                        )
                    )
        for ctx in contexts:
            sel = ok if ctx == "all" else ok & (joint.context.astype(str) == ctx)
            d = (xa - xb)[sel].to_numpy()
            delta_rows.append(
                (
                    metric,
                    ctx,
                    int(sel.sum()),
                    float(np.median(xa[sel])) if sel.any() else np.nan,
                    float(np.median(xb[sel])) if sel.any() else np.nan,
                    *(
                        np.quantile(d, [0.05, 0.25, 0.5, 0.75, 0.95])
                        if len(d)
                        else [np.nan] * 5
                    ),
                    float((d > 0).mean()) if len(d) else np.nan,
                )
            )
    joint_grid = pd.DataFrame(
        grid_rows,
        columns=[
            "metric",
            f"sd_{labels[0]}_lo",
            f"sd_{labels[0]}_hi",
            f"sd_{labels[1]}_lo",
            f"sd_{labels[1]}_hi",
            "count",
        ],
    )
    joint_grid.to_csv(
        f"{p}.joint_sd_grid.tsv", sep="\t", index=False, float_format="%.6g"
    )
    delta = pd.DataFrame(
        delta_rows,
        columns=[
            "metric",
            "context",
            "n_sites",
            f"median_sd_{labels[0]}",
            f"median_sd_{labels[1]}",
            "delta_q05",
            "delta_q25",
            "delta_median",
            "delta_q75",
            "delta_q95",
            f"frac_sites_{labels[0]}_more_variable",
        ],
    )
    delta.to_csv(
        f"{p}.joint_delta_sd_summary.tsv", sep="\t", index=False, float_format="%.6g"
    )
    del joint

    # Joint PCA on the shared PCA site subset.
    pa = np.load(args.pca_matrix[0], allow_pickle=False)
    pb = np.load(args.pca_matrix[1], allow_pickle=False)
    key_a = pd.Index(
        pd.Series(pa["chrom"]).astype(str) + ":" + pd.Series(pa["pos"]).astype(str)
    )
    key_b = pd.Index(
        pd.Series(pb["chrom"]).astype(str) + ":" + pd.Series(pb["pos"]).astype(str)
    )
    shared = key_a.intersection(key_b)
    ia, ib = key_a.get_indexer(shared), key_b.get_indexer(shared)
    auto = np.array([is_autosome(k.split(":")[0]) for k in shared])
    mat = np.vstack([pa["beta"][:, ia[auto]], pb["beta"][:, ib[auto]]])
    ids = list(pa["sample_ids"]) + list(pb["sample_ids"])
    plats = [labels[0]] * len(pa["sample_ids"]) + [labels[1]] * len(pb["sample_ids"])
    coords, var_frac, n_used = run_pca(mat)
    write_pca(f"{p}.joint", ids, plats, coords, var_frac, n_used)

    qc = merged["sample_qc"]
    plot_compare_sd(
        plt, merged["site_sd_hist"], labels, f"{p}.compare_site_sd", written
    )
    plot_compare_sd_by_mean(
        plt, merged["site_sd_by_mean"], labels, f"{p}.compare_site_sd_by_mean", written
    )
    plot_joint_grid(plt, joint_grid, labels, f"{p}.compare_joint_site_sd", written)
    plot_sample_qc(plt, qc, f"{p}.compare_sample_qc_metrics", written)
    plot_beta_distributions(
        plt, merged["sample_hists"], f"{p}.compare_beta_distributions", written
    )
    plot_coverage(plt, merged["sample_hists"], qc, f"{p}.compare_coverage", written)
    plot_neighbor(
        plt, merged["sample_hists"], qc, f"{p}.compare_neighbor_discordance", written
    )
    plot_intermediate_vs_cov(plt, qc, f"{p}.compare_intermediate_vs_coverage", written)
    plot_compare_excess(
        plt,
        merged["site_excess_var_hist"],
        labels,
        f"{p}.compare_excess_variance",
        written,
    )
    plot_pca(
        plt,
        read_tsv(f"{p}.joint.pca_coords.tsv"),
        read_tsv(f"{p}.joint.pca_variance.tsv"),
        qc,
        f"{p}.compare_joint_pca",
        written,
    )
    log("\n".join(written))


def plot_compare_sd(plt, sd_hist, labels, out, written):
    """The headline figure: per-CpG SD across samples, by platform."""
    contexts = [c for c in ["all"] + CONTEXTS if c in set(sd_hist.context)]
    fig, axes = plt.subplots(
        2,
        len(contexts),
        figsize=(max(3.6 * len(contexts), 7), 7),
        squeeze=False,
        sharey="row",
    )
    for col_i, ctx in enumerate(contexts):
        for row_i, metric in enumerate(("raw", "matched")):
            ax = axes[row_i, col_i]
            for li, lab in enumerate(labels):
                h = sd_hist[sd_hist.platform == lab]
                x, _, cdf = sd_density(h, metric, ctx)
                color = PLATFORM_COLORS.get(lab, SERIES[li])
                ax.plot(x, cdf, color=color, label=lab)
                if len(cdf) and cdf[-1] > 0:
                    med = x[np.searchsorted(cdf, 0.5)]
                    ax.axvline(med, color=color, linewidth=1, linestyle=":")
            ax.set_xlim(0, 0.35)
            ax.set_title(f"{ctx} - {metric}", fontsize=9)
            if row_i == 1:
                ax.set_xlabel("Per-site SD of beta across samples")
            if col_i == 0:
                ax.set_ylabel("Cumulative fraction of sites")
    axes[0, 0].legend(loc="lower right")
    fig.suptitle(
        "Per-CpG methylation variation by platform (ECDF; dotted = median)\n"
        "Top: raw. Bottom: every sample downsampled to the same depth"
    )
    savefig(fig, out, written)
    plt.close(fig)


def plot_compare_sd_by_mean(plt, by_mean, labels, out, written):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, metric in zip(axes, ("raw", "matched")):
        for li, lab in enumerate(labels):
            bm = by_mean[
                (by_mean.platform == lab)
                & (by_mean.metric == metric)
                & (by_mean.context == "all")
                & (by_mean.mean_bin != "all")
            ]
            mid = (bm.mean_lo + bm.mean_hi) / 2
            color = PLATFORM_COLORS.get(lab, SERIES[li])
            ax.fill_between(
                mid, bm.sd_q25, bm.sd_q75, color=color, alpha=0.2, linewidth=0
            )
            ax.plot(
                mid,
                bm.sd_median,
                color=color,
                marker="o",
                markersize=4,
                label=f"{lab} median, IQR band",
            )
            if metric == "raw":
                ax.plot(
                    mid,
                    bm.binom_sd_mean,
                    color=color,
                    linestyle="--",
                    linewidth=1.2,
                    label=f"{lab} expected binomial SD",
                )
        ax.set_title(f"{metric} beta")
        ax.set_xlabel("Per-site mean beta")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("Per-site SD of beta")
    fig.suptitle("Per-CpG SD by methylation level (autosomes)")
    savefig(fig, out, written)
    plt.close(fig)


def plot_joint_grid(plt, grid, labels, out, written):
    from matplotlib.colors import LogNorm

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    for ax, metric in zip(axes, ("raw", "matched")):
        g = grid[grid.metric == metric]
        if g.empty:
            ax.set_visible(False)
            continue
        mat = np.zeros((GRID_BINS, GRID_BINS))
        i = np.clip(
            np.rint(g.iloc[:, 1].to_numpy() * 2 * GRID_BINS).astype(int),
            0,
            GRID_BINS - 1,
        )
        j = np.clip(
            np.rint(g.iloc[:, 3].to_numpy() * 2 * GRID_BINS).astype(int),
            0,
            GRID_BINS - 1,
        )
        np.add.at(mat, (i, j), g["count"].to_numpy())
        im = ax.imshow(
            np.ma.masked_equal(mat, 0).T,
            origin="lower",
            extent=[0, 0.5, 0, 0.5],
            cmap=sequential_cmap(),
            norm=LogNorm(),
        )
        ax.plot([0, 0.5], [0, 0.5], color=TEXT_SECONDARY, linewidth=1, linestyle="--")
        fig.colorbar(im, ax=ax, label="CpG sites")
        ax.set_xlabel(f"Per-site SD, {labels[0]}")
        ax.set_ylabel(f"Per-site SD, {labels[1]}")
        ax.set_title(f"{metric} beta (below diagonal: {labels[0]} more variable)")
    fig.suptitle("Same CpG, both platforms (autosomal sites covered in both cohorts)")
    savefig(fig, out, written)
    plt.close(fig)


def plot_compare_excess(plt, excess, labels, out, written):
    fig, ax = plt.subplots(figsize=(7, 4))
    for li, lab in enumerate(labels):
        h = excess[(excess.platform == lab) & (excess.context == "all")]
        tot = h["count"].sum()
        if tot:
            ax.plot(
                (h.bin_lo + h.bin_hi) / 2,
                h["count"] / tot,
                color=PLATFORM_COLORS.get(lab, SERIES[li]),
                label=lab,
            )
    ax.axvline(0, color=TEXT_SECONDARY, linewidth=1)
    ax.set_xlabel("Excess variance (observed - expected binomial)")
    ax.set_ylabel("Fraction of autosomal sites")
    ax.legend()
    ax.set_title("Variance beyond read-sampling noise")
    savefig(fig, out, written)
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("build-index")
    p.add_argument("--ref-fa", required=True)
    p.add_argument("--ref-fai", required=True)
    p.add_argument("--contigs", nargs="+", required=True)
    p.add_argument("--cpg-islands-bed")
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_build_index)

    p = sub.add_parser("sample")
    p.add_argument("--bed", required=True)
    p.add_argument("--sample-id", required=True)
    p.add_argument("--platform", required=True, choices=sorted(PLATFORM_FORMATS))
    p.add_argument("--cpg-index", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--min-coverage", type=int, default=10)
    p.add_argument("--matched-coverage", type=int, default=10)
    p.add_argument("--max-neighbor-distance", type=int, default=50)
    p.add_argument("--min-ref-cpg-match-frac", type=float, default=0.8)
    p.add_argument(
        "--pacbio-beta-source",
        choices=["model", "counts"],
        default="model",
        help="PacBio raw beta: pb-cpg-tools mod_score (model) or "
        "est_mod_count / cov (counts). Downsampling always uses counts.",
    )
    p.set_defaults(func=cmd_sample)

    p = sub.add_parser("accumulate")
    p.add_argument("--cpg-index", required=True)
    p.add_argument("--samples", nargs="+", required=True)
    p.add_argument("--min-coverage", type=int, default=10)
    p.add_argument("--pca-site-stride", type=int, default=100)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_accumulate)

    p = sub.add_parser("finalize")
    p.add_argument("--cpg-index", required=True)
    p.add_argument("--partials", nargs="+", required=True)
    p.add_argument("--sample-qc", required=True)
    p.add_argument("--platform", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--min-coverage", type=int, default=10)
    p.add_argument("--matched-coverage", type=int, default=10)
    p.add_argument("--min-samples-per-site", type=int, default=10)
    p.add_argument("--outlier-mad", type=float, default=4.0)
    p.set_defaults(func=cmd_finalize)

    p = sub.add_parser("plot")
    for name in (
        "sample-qc",
        "sample-hists",
        "site-sd-hist",
        "site-sd-by-mean",
        "site-mean-sd-grid",
        "site-excess-var-hist",
        "pca-coords",
        "pca-variance",
    ):
        p.add_argument(f"--{name}", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--outdir", default="figures")
    p.set_defaults(func=cmd_plot)

    p = sub.add_parser("compare")
    for name in (
        "sample-qc",
        "sample-hists",
        "site-sd-hist",
        "site-sd-by-mean",
        "site-mean-sd-grid",
        "site-excess-var-hist",
        "cohort-summary",
        "per-site-stats",
        "pca-matrix",
    ):
        p.add_argument(f"--{name}", nargs=2, required=True, metavar=("A", "B"))
    p.add_argument("--label-a", default="PacBio")
    p.add_argument("--label-b", default="ONT")
    p.add_argument("--min-samples-per-site", type=int, default=10)
    p.add_argument("--prefix", required=True)
    p.add_argument("--outdir", default="figures")
    p.set_defaults(func=cmd_compare)
    return ap


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
