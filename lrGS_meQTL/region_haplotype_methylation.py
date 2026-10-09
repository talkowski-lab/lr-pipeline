#!/usr/bin/env python3
"""Summarise per-haplotype CpG methylation within target regions.

Subcommands:
  contig     for one contig of a wide per-haplotype methylation table (#chrom start end <sample>_hap1 ... ; percent
             methylation, '.' = missing), take every target region on that contig with >= --min-cpg CpG sites and
             write, per haplotype column, the number of called CpGs, mean, median and SD (SD needs >= 2 called CpGs);
             also writes per-haplotype CpG-site statistics for sample QC; samples in --excluded-samples are removed
             first, so they are absent from every output
  sample-qc  combine the per-contig site statistics over autosomes and flag poor-quality samples
  outliers   per region, call samples with outlying methylation relative to the QC-passing cohort, following the
             MethBat background-comparison rule (see cmd_outliers)
  concat     concatenate gzipped TSVs that share a header (header written once)

contig outputs (<prefix> = --prefix):
  <prefix>.regions.tsv.gz    region_index (0-based row of the region in --regions, header excluded), the region's
                             original columns, region_n_cpg (CpG sites in the methylation table within
                             [start, end)), region_n_haplotypes_called (haplotypes with >= 1 called CpG in the region)
  <prefix>.n_called.tsv.gz, <prefix>.mean.tsv.gz, <prefix>.median.tsv.gz, <prefix>.sd.tsv.gz
                             region_index, chrom, start, end, then one column per haplotype (NA = not computable)
  <prefix>.site_qc_stats.tsv.gz   per haplotype column: contig, n_sites, n_called, sum, n_intermediate (20-80%),
                             correlation sums against the per-CpG cohort median, n_both / sum_hap_diff (sample level)
  <prefix>.site_hist.tsv.gz  per haplotype column: counts of called values in 0.1%-wide bins (for medians)
"""
import argparse
import gzip
import re
import warnings

import numpy as np
import pandas as pd


def load_matrix(path, chunksize=200000):
    """Return (chrom, CpG start positions, haplotype column names, float32 matrix CpGs x haplotypes)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").lstrip("#").split("\t")
    names = header[3:]
    dtypes = {c: np.float32 for c in names}
    dtypes.update({header[0]: str, header[1]: np.int64, header[2]: np.int64})
    chroms, pos, blocks = [], [], []
    reader = pd.read_csv(path, sep="\t", na_values=["."], dtype=dtypes, chunksize=chunksize, engine="c")
    for chunk in reader:
        chroms.append(chunk.iloc[:, 0].to_numpy())
        pos.append(chunk.iloc[:, 1].to_numpy(np.int64))
        blocks.append(chunk.iloc[:, 3:].to_numpy(np.float32))
    if not blocks:
        return np.array([], dtype=str), np.array([], dtype=np.int64), names, np.zeros((0, len(names)), np.float32)
    pos = np.concatenate(pos)
    mat = np.concatenate(blocks)
    order = np.argsort(pos, kind="stable")
    if not np.all(order == np.arange(len(pos))):
        pos, mat = pos[order], mat[order]
    return np.concatenate(chroms), pos, names, mat


def load_regions(path):
    """Return the region table (all columns as text) with a region_index column, and its column names."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        first = fh.readline().rstrip("\n")
    has_header = first.startswith("#")
    regions = pd.read_csv(path, sep="\t", header=None, skiprows=1 if has_header else 0, dtype=str,
                          keep_default_na=False)
    cols = first.lstrip("#").split("\t") if has_header else \
        ["chrom", "start", "end"] + [f"col{i}" for i in range(4, regions.shape[1] + 1)]
    regions.columns = cols
    regions.insert(0, "region_index", np.arange(len(regions)))
    return regions, cols


def haplotype_pairs(names):
    """Map sample -> (hap1 column index or None, hap2 column index or None) from <sample>_hap1/_hap2 names."""
    pairs = {}
    for j, n in enumerate(names):
        m = re.match(r"^(.*)_hap([12])$", n)
        if m:
            pairs.setdefault(m.group(1), [None, None])[int(m.group(2)) - 1] = j
    return pairs


def read_excluded(paths):
    """Union of sample IDs from exclusion lists (first tab-separated field per line; blank lines ignored)."""
    out = set()
    for path in paths or []:
        with open(path) as fh:
            out |= {line.split("\t")[0].strip() for line in fh if line.strip()}
    return out


def drop_samples(names, mat, excluded):
    """Remove both haplotype columns of excluded samples."""
    keep = [j for j, n in enumerate(names) if re.sub(r"_hap[12]$", "", n) not in excluded]
    return [names[j] for j in keep], mat[:, keep]


def robust_z(values):
    v = np.asarray(values, dtype=float)
    med = np.nanmedian(v)
    mad = 1.4826 * np.nanmedian(np.abs(v - med))
    if not np.isfinite(mad) or mad == 0:
        mad = np.nanstd(v) if np.nanstd(v) > 0 else 1.0
    return (v - med) / mad


def site_qc_stats(contig, names, mat):
    """Per haplotype column CpG-site statistics for one contig (summed across contigs in sample-qc).

    Works one column at a time so memory stays close to the size of the matrix itself.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ref = np.nanmedian(mat, axis=1) if len(mat) else np.zeros(0, dtype=np.float32)
    ref_ok = np.isfinite(ref)
    ref64 = ref.astype(np.float64)
    rows, hist = [], np.zeros((len(names), 1001), dtype=np.int64)
    for j in range(len(names)):
        v = mat[:, j]
        called = ~np.isnan(v)
        x = v[called].astype(np.float64)
        pair = called & ref_ok
        xp, y = v[pair].astype(np.float64), ref64[pair]
        rows.append({"n_called": int(called.sum()), "sum": x.sum(), "n_intermediate": int(((x >= 20) & (x <= 80)).sum()),
                     "n_pair": int(pair.sum()), "sx": xp.sum(), "sy": y.sum(), "sxx": (xp * xp).sum(),
                     "syy": (y * y).sum(), "sxy": (xp * y).sum()})
        hist[j] = np.bincount(np.clip(np.round(x * 10).astype(np.int64), 0, 1000), minlength=1001)
    stats = pd.DataFrame(rows)
    stats.insert(0, "n_sites", len(mat))
    stats.insert(0, "contig", contig)
    stats.insert(0, "column", names)
    stats["n_both"], stats["sum_hap_diff"] = np.nan, np.nan
    for sample, (j1, j2) in haplotype_pairs(names).items():
        if j1 is None or j2 is None:
            continue
        both = ~np.isnan(mat[:, j1]) & ~np.isnan(mat[:, j2])
        stats.loc[[j1, j2], "n_both"] = both.sum()
        stats.loc[j1, "sum_hap_diff"] = (mat[both, j1].astype(np.float64) - mat[both, j2]).sum()
    hist_df = pd.DataFrame(hist, columns=[f"b{k}" for k in range(1001)])
    hist_df.insert(0, "contig", contig)
    hist_df.insert(0, "column", names)
    return stats, hist_df


def cmd_contig(args):
    regions, cols = load_regions(args.regions)
    chrom_col, start_col, end_col = cols[0], cols[1], cols[2]
    regions = regions[regions[chrom_col] == args.contig].copy()
    _, pos, names, mat = load_matrix(args.bed)
    excluded = read_excluded(args.excluded_samples)
    n_before = len(names)
    names, mat = drop_samples(names, mat, excluded)
    if excluded:
        print(f"{args.contig}: removed {n_before - len(names)} haplotype columns of {len(excluded)} excluded samples")
    qc_stats, qc_hist = site_qc_stats(args.contig, names, mat)
    qc_stats.to_csv(f"{args.prefix}.site_qc_stats.tsv.gz", sep="\t", index=False, na_rep="NA", compression="gzip")
    qc_hist.to_csv(f"{args.prefix}.site_hist.tsv.gz", sep="\t", index=False, compression="gzip")
    starts = regions[start_col].astype(np.int64).to_numpy()
    ends = regions[end_col].astype(np.int64).to_numpy()
    i0 = np.searchsorted(pos, starts, side="left")
    i1 = np.searchsorted(pos, ends, side="left")
    n_cpg = i1 - i0
    keep = n_cpg >= args.min_cpg
    regions, i0, i1, n_cpg = regions[keep], i0[keep], i1[keep], n_cpg[keep]

    n_hap = len(names)
    out = {k: np.full((len(regions), n_hap), np.nan, dtype=np.float32) for k in ("n_called", "mean", "median", "sd")}
    cache = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for r, (a, b) in enumerate(zip(i0, i1)):
            key = (a, b)
            if key not in cache:
                block = mat[a:b]
                n = (~np.isnan(block)).sum(axis=0).astype(np.float32)
                sd = np.nanstd(block, axis=0, ddof=1)
                sd[n < 2] = np.nan
                cache = {key: (n, np.nanmean(block, axis=0), np.nanmedian(block, axis=0), sd)}
            n, mean, median, sd = cache[key]
            out["n_called"][r], out["mean"][r], out["median"][r], out["sd"][r] = n, mean, median, sd

    regions["region_n_cpg"] = n_cpg
    regions["region_n_haplotypes_called"] = (out["n_called"] > 0).sum(axis=1)
    regions.rename(columns={chrom_col: "#" + chrom_col}).to_csv(
        f"{args.prefix}.regions.tsv.gz", sep="\t", index=False, compression="gzip")
    key_cols = regions[["region_index", chrom_col, start_col, end_col]].reset_index(drop=True)
    key_cols.columns = ["region_index", "#chrom", "start", "end"]
    for stat, arr in out.items():
        df = pd.concat([key_cols, pd.DataFrame(arr, columns=names)], axis=1)
        df.to_csv(f"{args.prefix}.{stat}.tsv.gz", sep="\t", index=False, na_rep="NA",
                  float_format="%.0f" if stat == "n_called" else "%.2f", compression="gzip")
    print(f"{args.contig}: {int(keep.sum()):,} of {len(keep):,} regions with >= {args.min_cpg} CpG sites; "
          f"{len(pos):,} CpG sites x {n_hap} haplotypes")


def cmd_sample_qc(args):
    """Autosomal per-haplotype metrics; a sample is excluded if either haplotype is an outlier or the haplotypes disagree.

    Haplotype-level flags (robust z across all haplotype columns, |z| >= --z unless noted): mean, median, fraction of
    intermediate CpGs (20-80%), correlation with the per-CpG cohort median (low side only: z <= -z), and call rate
    below --min-call-rate-frac x the cohort median call rate. Sample-level flag: autosomal hap1 - hap2 mean over CpGs
    called on both haplotypes with |robust z| >= --z.
    """
    auto_re = re.compile(args.autosome_regex)
    st = pd.concat([pd.read_csv(f, sep="\t") for f in args.stats], ignore_index=True)
    st = st[st["contig"].astype(str).apply(lambda c: bool(auto_re.match(c)))]
    hist = pd.concat([pd.read_csv(f, sep="\t") for f in args.hists], ignore_index=True)
    hist = hist[hist["contig"].astype(str).apply(lambda c: bool(auto_re.match(c)))]
    g = st.groupby("column", sort=False)
    sums = g[["n_sites", "n_called", "sum", "n_intermediate", "n_pair", "sx", "sy", "sxx", "syy", "sxy",
              "n_both", "sum_hap_diff"]].sum(min_count=1)
    qc = pd.DataFrame(index=sums.index)
    qc["call_rate"] = sums["n_called"] / sums["n_sites"]
    qc["mean"] = sums["sum"] / sums["n_called"]
    h = hist.drop(columns="contig").groupby("column", sort=False).sum().reindex(qc.index).to_numpy()
    cum = h.cumsum(axis=1)
    half = cum[:, -1:] / 2
    qc["median"] = np.where(cum[:, -1] > 0, np.argmax(cum >= half, axis=1) / 10, np.nan)
    qc["frac_intermediate"] = sums["n_intermediate"] / sums["n_called"]
    n = sums["n_pair"]
    cov = sums["sxy"] - sums["sx"] * sums["sy"] / n
    var = (sums["sxx"] - sums["sx"] ** 2 / n) * (sums["syy"] - sums["sy"] ** 2 / n)
    qc["corr_to_cpg_median"] = cov / np.sqrt(var)
    for m in ("mean", "median", "frac_intermediate", "corr_to_cpg_median"):
        qc[f"z_{m}"] = robust_z(qc[m])
    qc["call_rate_vs_cohort_median"] = qc["call_rate"] / qc["call_rate"].median()
    reasons = []
    for _, r in qc.iterrows():
        why = [m for m in ("mean", "median", "frac_intermediate") if abs(r[f"z_{m}"]) >= args.z]
        if r["z_corr_to_cpg_median"] <= -args.z:
            why.append("low_corr_to_cohort")
        if r["call_rate_vs_cohort_median"] < args.min_call_rate_frac:
            why.append("low_call_rate")
        reasons.append(",".join(why))
    qc["haplotype_flags"] = reasons
    qc.index.name = "column"
    qc = qc.reset_index()
    qc["sample"] = qc["column"].str.replace(r"_hap[12]$", "", regex=True)
    qc["haplotype"] = qc["column"].str.extract(r"_(hap[12])$")[0]

    diff = (sums["sum_hap_diff"].groupby(sums.index.str.replace(r"_hap[12]$", "", regex=True)).sum(min_count=1) /
            sums["n_both"].groupby(sums.index.str.replace(r"_hap[12]$", "", regex=True)).max())
    samp = pd.DataFrame({"hap1_minus_hap2_mean": diff})
    samp["z_hap_mean_diff"] = robust_z(samp["hap1_minus_hap2_mean"])
    flags = qc[qc["haplotype_flags"] != ""].groupby("sample")
    samp["exclude_reasons"] = ""
    for sample, rows in flags:
        samp.loc[sample, "exclude_reasons"] = ";".join(f"{r.haplotype}:{r.haplotype_flags}" for r in rows.itertuples())
    imb = samp["z_hap_mean_diff"].abs() >= args.z
    samp.loc[imb, "exclude_reasons"] = (samp.loc[imb, "exclude_reasons"] + ";haplotype_mean_imbalance").str.lstrip(";")
    samp["exclude"] = samp["exclude_reasons"] != ""
    samp.index.name = "sample"
    qc.merge(samp.reset_index(), on="sample", how="left").to_csv(
        f"{args.prefix}.sample_qc.tsv", sep="\t", index=False, na_rep="NA", float_format="%.5g")
    with open(f"{args.prefix}.excluded_samples.txt", "w") as fh:
        for sample, r in samp[samp["exclude"]].iterrows():
            fh.write(f"{sample}\t{r['exclude_reasons']}\n")
    print(f"sample QC: {int(samp['exclude'].sum())} of {len(samp)} samples excluded")


def cmd_outliers(args):
    """Per region, label samples whose methylation is unusual relative to the QC-passing cohort (MethBat rule).

    Haplotype call-rate filter: a haplotype passes in a region if it has a called value at >= --min-hap-call-frac of
    the region's CpG sites (n_called / region_n_cpg, from --n-called-matrix and --regions-table); a sample is examined
    in a region only if both of its haplotypes pass, otherwise it is left out of the background and of the calls.
    For each examined sample: combined = mean of the hap1 and hap2 region means; abs_hap_delta = |hap2 - hap1|. The
    background for a region is the mean and SD (ddof=1) of each quantity across examined QC-passing samples; labels
    need >= --min-samples examined samples.
      HyperMethylated / HypoMethylated: z(combined) >= --min-z / <= -min-z and combined - background mean
                                        >= --min-delta / <= -min-delta
      HyperASM / HypoASM:               the same on abs_hap_delta
    As in MethBat, the sample's own region state must match the direction: HyperMethylated needs combined >=
    --methylated-min (MethBat "Methylated", 80%), HypoMethylated needs combined <= --unmethylated-max ("Unmethylated",
    20%), HyperASM needs abs_hap_delta >= --asm-min-abs-delta ("AlleleSpecificMethylation", 50%; MethBat's extra
    Fisher's test on read counts is not possible from the haplotype tables). --no-state-requirement drops these checks.
    Hyper/HypoMethylated take priority over the ASM labels, and HypoASM has the lowest priority.
    """
    excluded = read_excluded(args.excluded_samples)
    mean = pd.read_csv(args.mean_matrix, sep="\t")
    key = mean.iloc[:, :4]
    names = list(mean.columns[4:])
    vals = mean.iloc[:, 4:].to_numpy(np.float64)
    pairs = {s: p for s, p in haplotype_pairs(names).items() if s not in excluded}
    samples = list(pairs)

    def col(j):
        return vals[:, j] if j is not None else np.full(len(vals), np.nan)

    h1 = np.column_stack([col(pairs[s][0]) for s in samples]) if samples else np.zeros((len(vals), 0))
    h2 = np.column_stack([col(pairs[s][1]) for s in samples]) if samples else np.zeros((len(vals), 0))
    ncalled = pd.read_csv(args.n_called_matrix, sep="\t")
    if not np.array_equal(ncalled["region_index"].to_numpy(), mean["region_index"].to_numpy()):
        raise SystemExit("n_called and mean matrices list different regions")
    regions = pd.read_csv(args.regions_table, sep="\t", usecols=["region_index", "region_n_cpg"])
    n_cpg = mean[["region_index"]].merge(regions, on="region_index", how="left")["region_n_cpg"].to_numpy(np.float64)
    nc = ncalled.iloc[:, 4:]

    def call_frac(j):
        return nc.iloc[:, j].to_numpy(np.float64) / n_cpg if j is not None else np.zeros(len(vals))

    passed = np.column_stack([(call_frac(pairs[s][0]) >= args.min_hap_call_frac) &
                              (call_frac(pairs[s][1]) >= args.min_hap_call_frac) for s in samples]) \
        if samples else np.zeros((len(vals), 0), dtype=bool)
    h1 = np.where(passed, h1, np.nan)
    h2 = np.where(passed, h2, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        combined = np.nanmean(np.stack([h1, h2]), axis=0)
        absd = np.abs(h2 - h1)
        out = key.copy()
        out["region_n_cpg"] = n_cpg.astype(int)
        out["n_samples_both_haplotypes_pass"] = passed.sum(1)
        calls = []
        for name, x in (("combined", combined), ("abs_hap_delta", absd)):
            n = (~np.isnan(x)).sum(1)
            mu = np.nanmean(x, axis=1)
            sd = np.nanstd(x, axis=1, ddof=1)
            out[f"n_samples_{name}"] = n
            out[f"pop_mean_{name}"] = mu
            out[f"pop_sd_{name}"] = sd
            z = (x - mu[:, None]) / sd[:, None]
            z[(n < args.min_samples) | ~(sd > 0)] = np.nan
            calls.append((x, mu, z))
    (cx, cmu, cz), (ax, amu, az) = calls
    cdelta, adelta = cx - cmu[:, None], ax - amu[:, None]
    label = np.full(cx.shape, "", dtype=object)
    with np.errstate(invalid="ignore"):
        if args.no_state_requirement:
            meth_ok = unmeth_ok = asm_ok = np.ones(cx.shape, dtype=bool)
        else:
            meth_ok, unmeth_ok, asm_ok = cx >= args.methylated_min, cx <= args.unmethylated_max, ax >= args.asm_min_abs_delta
        label[(az <= -args.min_z) & (adelta <= -args.min_delta)] = "HypoASM"
        label[(az >= args.min_z) & (adelta >= args.min_delta) & asm_ok] = "HyperASM"
        label[(cz <= -args.min_z) & (cdelta <= -args.min_delta) & unmeth_ok] = "HypoMethylated"
        label[(cz >= args.min_z) & (cdelta >= args.min_delta) & meth_ok] = "HyperMethylated"
    labels = ("HyperMethylated", "HypoMethylated", "HyperASM", "HypoASM")
    samples_arr = np.array(samples, dtype=object)
    for lab in labels:
        hit = label == lab
        out[f"n_{lab}"] = hit.sum(1)
        out[f"samples_{lab}"] = [";".join(samples_arr[row]) for row in hit]
    out["n_outlier_samples"] = (label != "").sum(1)
    out["has_outlier"] = out["n_outlier_samples"] > 0
    out.to_csv(f"{args.prefix}.region_outliers.tsv.gz", sep="\t", index=False, na_rep="NA", float_format="%.3f",
               compression="gzip")
    r, c = np.nonzero(label != "")
    long = key.iloc[r].reset_index(drop=True)
    long["sample"] = samples_arr[c]
    long["label"] = label[r, c]
    long["combined"], long["combined_delta"], long["combined_z"] = cx[r, c], cdelta[r, c], cz[r, c]
    long["abs_hap_delta"], long["abs_hap_delta_delta"], long["abs_hap_delta_z"] = ax[r, c], adelta[r, c], az[r, c]
    long.to_csv(f"{args.prefix}.outlier_calls.tsv.gz", sep="\t", index=False, na_rep="NA", float_format="%.3f",
                compression="gzip")
    print(f"{len(out):,} regions; {int(out['has_outlier'].sum()):,} with >= 1 outlier sample; {len(long):,} outlier calls; "
          f"{len(samples)} QC-passing samples ({len(excluded)} excluded); sample x region pairs with both haplotypes "
          f">= {args.min_hap_call_frac:.0%} called: {passed.mean():.1%}")


def cmd_concat(args):
    with gzip.open(args.out, "wt") as out:
        wrote_header = False
        for path in args.inputs:
            with gzip.open(path, "rt") as fh:
                header = fh.readline()
                if not wrote_header:
                    out.write(header)
                    wrote_header = True
                for line in fh:
                    out.write(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("contig")
    c.add_argument("--bed", required=True, help="per-contig wide per-haplotype methylation table (bgzip)")
    c.add_argument("--regions", required=True, help="target regions (BED-like, optional '#' header, optional gzip)")
    c.add_argument("--contig", required=True)
    c.add_argument("--min-cpg", type=int, default=1)
    c.add_argument("--excluded-samples", nargs="*", default=[],
                   help="sample lists to remove before all statistics (sample<TAB>anything)")
    c.add_argument("--prefix", required=True)
    c.set_defaults(func=cmd_contig)
    q = sub.add_parser("sample-qc")
    q.add_argument("--stats", nargs="+", required=True, help="per-contig *.site_qc_stats.tsv.gz")
    q.add_argument("--hists", nargs="+", required=True, help="per-contig *.site_hist.tsv.gz")
    q.add_argument("--z", type=float, default=5.0)
    q.add_argument("--min-call-rate-frac", type=float, default=0.5)
    q.add_argument("--autosome-regex", default=r"^chr[0-9]+$")
    q.add_argument("--prefix", required=True)
    q.set_defaults(func=cmd_sample_qc)
    o = sub.add_parser("outliers")
    o.add_argument("--mean-matrix", required=True, help="per-contig or genome-wide *.mean.tsv.gz from contig")
    o.add_argument("--n-called-matrix", required=True, help="matching *.n_called.tsv.gz from contig")
    o.add_argument("--regions-table", required=True, help="matching *.regions.tsv.gz from contig (region_n_cpg)")
    o.add_argument("--min-hap-call-frac", type=float, default=0.8,
                   help="a haplotype passes in a region if called at >= this fraction of the region's CpGs")
    o.add_argument("--excluded-samples", nargs="*", default=[],
                   help="sample lists to leave out of the background and the calls (e.g. sample-qc output)")
    o.add_argument("--min-z", type=float, default=3.0)
    o.add_argument("--min-delta", type=float, default=20.0)
    o.add_argument("--methylated-min", type=float, default=80.0)
    o.add_argument("--unmethylated-max", type=float, default=20.0)
    o.add_argument("--asm-min-abs-delta", type=float, default=50.0)
    o.add_argument("--no-state-requirement", action="store_true",
                   help="label outliers on z-score and delta alone, without MethBat's sample-state check")
    o.add_argument("--min-samples", type=int, default=10)
    o.add_argument("--prefix", required=True)
    o.set_defaults(func=cmd_outliers)
    m = sub.add_parser("concat")
    m.add_argument("--inputs", nargs="+", required=True)
    m.add_argument("--out", required=True)
    m.set_defaults(func=cmd_concat)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
