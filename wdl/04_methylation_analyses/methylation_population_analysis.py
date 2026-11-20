#!/usr/bin/env python3
"""Population-level analysis of long-read CpG methylation tables.

Input tables are gzipped, tab-separated, wide-format beds: `#chrom`, `start`,
`end`, then one column per sample (per-sample tables) or per haplotype
(`<sample>_hap1` / `<sample>_hap2`, haplotype tables). Values are percent
methylation (0-100); missing values are `.`.

Subcommands (each run by one WDL task in MethylationPopulationAnalysis.wdl):
  contig-stats      per-sample summary statistics for one contig
  sample-qc         contig x sample mean/median tables + outlier sample calls
  variable-regions  regions whose methylation varies across samples more than
                    expected, per-sample hyper/hypo calls, per-region plots
  haplotype-stats   per-sample hap1 vs hap2 summary statistics for one contig
  haplotype-qc      haplotype-level outlier sample calls
  allele-specific   per-sample allele-specific methylation (ASM) regions and
                    their population-level merge, with plots
  concat            concatenate per-contig tables, keeping one header
"""

import argparse
import gzip
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

HYPER_COLOR = "#1f4fbf"
HYPO_COLOR = "#c62828"
BACKGROUND_COLOR = "#9e9e9e"


def set_font(font_ttf):
    font_manager.fontManager.addfont(font_ttf)
    name = font_manager.FontProperties(fname=font_ttf).get_name()
    plt.rcParams["font.family"] = name
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42


def load_matrix(path, chunksize=200000):
    """Return (chrom, pos, column_names, float32 matrix rows=CpGs x cols)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    names = header[3:]
    dtypes = {c: np.float32 for c in names}
    dtypes[header[0]] = str
    dtypes[header[1]] = np.int64
    dtypes[header[2]] = np.int64
    chroms, positions, blocks = [], [], []
    reader = pd.read_csv(path, sep="\t", na_values=["."], dtype=dtypes, chunksize=chunksize, engine="c")
    for chunk in reader:
        chroms.append(chunk.iloc[:, 0].to_numpy())
        positions.append(chunk.iloc[:, 1].to_numpy())
        blocks.append(chunk.iloc[:, 3:].to_numpy(np.float32))
    if not blocks:
        return np.array([], dtype=str), np.array([], dtype=np.int64), names, np.zeros((0, len(names)), np.float32)
    pos = np.concatenate(positions)
    X = np.vstack(blocks)
    order = np.argsort(pos, kind="stable")
    if not np.all(order == np.arange(len(pos))):
        pos, X = pos[order], X[order]
    return np.concatenate(chroms), pos, names, X


def robust_z(values):
    v = np.asarray(values, dtype=float)
    med = np.nanmedian(v)
    mad = 1.4826 * np.nanmedian(np.abs(v - med))
    if not np.isfinite(mad) or mad == 0:
        mad = np.nanstd(v) if np.nanstd(v) > 0 else 1.0
    return (v - med) / mad


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def wilcoxon_p(d):
    d = d[np.isfinite(d)]
    if len(d) < 2 or np.all(d == 0):
        return 1.0
    return float(wilcoxon(d, zero_method="wilcox", alternative="two-sided").pvalue)


def contig_sort_key(contig):
    c = contig.replace("chr", "")
    return (0, int(c), "") if c.isdigit() else (1, 0, c)


def make_tiles(pos, tile_cpgs, max_span):
    """Non-overlapping tiles of tile_cpgs consecutive CpGs; drop tiles spanning > max_span bp."""
    n_tiles = len(pos) // tile_cpgs
    idx = np.arange(n_tiles * tile_cpgs).reshape(n_tiles, tile_cpgs)
    if n_tiles == 0:
        return idx
    span = pos[idx[:, -1]] - pos[idx[:, 0]] + 1
    return idx[span <= max_span]


def merge_flagged(starts, ends, gap, signs=None):
    """Merge sorted flagged intervals within gap bp (and of equal sign when signs given).

    Returns a list of (start, end, member_indices)."""
    merged = []
    for i in range(len(starts)):
        if merged and starts[i] - merged[-1][1] <= gap and (signs is None or signs[i] == signs[merged[-1][2][0]]):
            merged[-1][1] = max(merged[-1][1], ends[i])
            merged[-1][2].append(i)
        else:
            merged.append([starts[i], ends[i], [i]])
    return [(s, e, m) for s, e, m in merged]


def write_empty_pdf(path, message):
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.axis("off")
        ax.text(0.5, 0.5, message, ha="center", va="center")
        pdf.savefig(fig)
        plt.close(fig)


# contig-stats

def cmd_contig_stats(args):
    _, pos, samples, X = load_matrix(args.bed)
    n_rows = X.shape[0]
    called = np.isfinite(X)
    n_called = called.sum(0)
    stats = pd.DataFrame({"sample": samples, "contig": args.contig, "n_cpg_rows": n_rows, "n_called": n_called})
    stats["call_rate"] = n_called / n_rows if n_rows else np.nan
    with np.errstate(invalid="ignore", divide="ignore"):
        stats["mean"] = np.nanmean(X, 0) if n_rows else np.nan
        stats["median"] = np.nanmedian(X, 0) if n_rows else np.nan
        stats["frac_low"] = ((X < 20) & called).sum(0) / n_called
        stats["frac_intermediate"] = ((X >= 20) & (X <= 80) & called).sum(0) / n_called
        stats["frac_high"] = ((X > 80) & called).sum(0) / n_called
    corr = np.full(len(samples), np.nan)
    n_corr = np.zeros(len(samples), dtype=np.int64)
    dense = called.mean(1) >= args.min_call_rate_for_corr if n_rows else np.zeros(0, bool)
    if dense.sum() > 100:
        Xd = X[dense]
        ref = np.nanmedian(Xd, 1)
        for j in range(len(samples)):
            ok = np.isfinite(Xd[:, j])
            if ok.sum() > 100 and np.std(Xd[ok, j]) > 0:
                corr[j] = np.corrcoef(Xd[ok, j], ref[ok])[0, 1]
                n_corr[j] = ok.sum()
    stats["corr_to_cpg_median"] = corr
    stats["n_corr"] = n_corr
    stats.to_csv(f"{args.prefix}.sample_stats.tsv", sep="\t", index=False, na_rep="NA")
    hist = np.zeros((len(samples), 1001), dtype=np.int64)
    for j in range(len(samples)):
        x = X[called[:, j], j]
        if len(x):
            hist[j] = np.bincount(np.clip(np.rint(x * 10).astype(np.int64), 0, 1000), minlength=1001)
    h = pd.DataFrame(hist, columns=[f"b{k}" for k in range(1001)])
    h.insert(0, "contig", args.contig)
    h.insert(0, "sample", samples)
    h.to_csv(f"{args.prefix}.sample_hist.tsv.gz", sep="\t", index=False)


# sample-qc

def cmd_sample_qc(args):
    set_font(args.font_ttf)
    long = pd.concat([pd.read_csv(f, sep="\t", dtype={"contig": str}) for f in args.stats], ignore_index=True)
    hists = pd.concat([pd.read_csv(f, sep="\t", dtype={"contig": str}) for f in args.hists], ignore_index=True)
    contigs = sorted(long["contig"].unique(), key=contig_sort_key)
    samples = list(dict.fromkeys(long["sample"]))
    auto_re = re.compile(args.autosome_regex)
    autosomes = [c for c in contigs if auto_re.match(c)]
    long.to_csv(f"{args.prefix}.per_contig_sample_stats.tsv", sep="\t", index=False, na_rep="NA")
    for metric, name in [("mean", "mean"), ("median", "median"), ("call_rate", "call_rate")]:
        tab = long.pivot(index="contig", columns="sample", values=metric).reindex(index=contigs, columns=samples)
        tab.to_csv(f"{args.prefix}.{name}_methylation_by_contig.tsv" if metric != "call_rate"
                   else f"{args.prefix}.call_rate_by_contig.tsv", sep="\t", na_rep="NA", float_format="%.4f")

    auto = long[long["contig"].isin(autosomes)].copy()
    g = auto.assign(w_mean=auto["mean"] * auto["n_called"],
                    w_int=auto["frac_intermediate"] * auto["n_called"],
                    w_corr=auto["corr_to_cpg_median"] * auto["n_corr"]).groupby("sample", sort=False)
    qc = pd.DataFrame(index=pd.Index(samples, name="sample"))
    qc["n_called_autosomal"] = g["n_called"].sum()
    qc["call_rate_autosomal"] = g["n_called"].sum() / g["n_cpg_rows"].sum()
    qc["mean_autosomal"] = g["w_mean"].sum() / g["n_called"].sum()
    qc["frac_intermediate_autosomal"] = g["w_int"].sum() / g["n_called"].sum()
    qc["corr_to_cpg_median_autosomal"] = g["w_corr"].sum() / g["n_corr"].sum()
    ah = hists[hists["contig"].isin(autosomes)].groupby("sample", sort=False)[[f"b{k}" for k in range(1001)]].sum()
    ah = ah.reindex(samples).to_numpy()
    cum = np.cumsum(ah, 1)
    qc["median_autosomal"] = [np.searchsorted(c, c[-1] / 2.0) / 10.0 if c[-1] > 0 else np.nan for c in cum]

    for contig in ["chrX", "chrY"]:
        sub = long[long["contig"] == contig].set_index("sample").reindex(samples)
        qc[f"{contig}_call_rate"] = sub["call_rate"].to_numpy()
        qc[f"{contig}_mean"] = sub["mean"].to_numpy()
    if qc["chrY_call_rate"].notna().any():
        y = qc["chrY_call_rate"].fillna(0)
        qc["inferred_sex"] = np.where(y >= args.male_chrY_frac * np.quantile(y, 0.95), "male", "female")
    else:
        qc["inferred_sex"] = "unknown"

    qc["z_mean"] = robust_z(qc["mean_autosomal"])
    qc["z_median"] = robust_z(qc["median_autosomal"])
    qc["z_frac_intermediate"] = robust_z(qc["frac_intermediate_autosomal"])
    qc["z_corr"] = robust_z(qc["corr_to_cpg_median_autosomal"])
    qc["z_call_rate"] = robust_z(qc["call_rate_autosomal"])
    contig_z = pd.DataFrame({c: robust_z(long[long["contig"] == c].set_index("sample").reindex(samples)["mean"])
                             for c in autosomes}, index=samples)
    qc["n_outlier_autosomes"] = (contig_z.abs() >= args.z).sum(1)
    contig_z.T.to_csv(f"{args.prefix}.mean_methylation_robust_z_by_contig.tsv", sep="\t", float_format="%.3f")

    reasons = []
    for s, r in qc.iterrows():
        why = []
        if abs(r["z_mean"]) >= args.z:
            why.append("mean")
        if abs(r["z_median"]) >= args.z:
            why.append("median")
        if abs(r["z_frac_intermediate"]) >= args.z:
            why.append("frac_intermediate")
        if r["z_corr"] <= -args.z:
            why.append("low_corr_to_cohort")
        if r["call_rate_autosomal"] < args.min_call_rate:
            why.append("low_call_rate")
        reasons.append(",".join(why))
    qc["exclude_reasons"] = reasons
    qc["exclude"] = qc["exclude_reasons"] != ""
    warnings = []
    for s, r in qc.iterrows():
        w = []
        if r["z_call_rate"] <= -args.z and "low_call_rate" not in r["exclude_reasons"]:
            w.append("relatively_low_call_rate")
        if r["n_outlier_autosomes"] >= args.min_outlier_autosomes_warn:
            w.append(f"outlier_on_{int(r['n_outlier_autosomes'])}_autosomes")
        warnings.append(",".join(w))
    qc["warnings"] = warnings
    qc.reset_index().to_csv(f"{args.prefix}.sample_qc.tsv", sep="\t", index=False, na_rep="NA", float_format="%.5g")
    with open(f"{args.prefix}.excluded_samples.txt", "w") as fh:
        for s in qc.index[qc["exclude"]]:
            fh.write(f"{s}\t{qc.loc[s, 'exclude_reasons']}\n")

    ex = qc["exclude"].to_numpy()
    with PdfPages(f"{args.prefix}.sample_qc.pdf") as pdf:
        fig, ax = plt.subplots(figsize=(6.5, 5))
        ax.scatter(qc["mean_autosomal"][~ex], qc["corr_to_cpg_median_autosomal"][~ex], s=12, c=BACKGROUND_COLOR, label="pass")
        ax.scatter(qc["mean_autosomal"][ex], qc["corr_to_cpg_median_autosomal"][ex], s=16, c=HYPO_COLOR, label="excluded")
        for s in qc.index[ex]:
            ax.annotate(s, (qc.loc[s, "mean_autosomal"], qc.loc[s, "corr_to_cpg_median_autosomal"]), fontsize=6)
        ax.set_xlabel("Autosomal mean methylation (%)")
        ax.set_ylabel("Correlation with cohort per-CpG median")
        ax.set_title(f"Sample QC ({ex.sum()} of {len(ex)} excluded)")
        ax.legend(frameon=False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        metrics = [("mean_autosomal", "Mean (%)"), ("median_autosomal", "Median (%)"),
                   ("frac_intermediate_autosomal", "Fraction CpGs 20-80%"),
                   ("corr_to_cpg_median_autosomal", "Corr. to cohort median"),
                   ("call_rate_autosomal", "Call rate"), ("n_outlier_autosomes", "# outlier autosomes")]
        fig, axes = plt.subplots(2, 3, figsize=(11, 6.5))
        rng = np.random.default_rng(0)
        jitter = rng.uniform(-0.3, 0.3, len(qc))
        for ax, (m, lab) in zip(axes.flat, metrics):
            ax.scatter(jitter[~ex], qc[m][~ex], s=8, c=BACKGROUND_COLOR)
            ax.scatter(jitter[ex], qc[m][ex], s=10, c=HYPO_COLOR)
            ax.set_xticks([])
            ax.set_ylabel(lab)
        fig.suptitle("Per-sample autosomal QC metrics (red = excluded)")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        mean_tab = long.pivot(index="contig", columns="sample", values="mean").reindex(index=contigs, columns=samples)
        fig, ax = plt.subplots(figsize=(11, 5))
        data = [mean_tab.loc[c, ~ex].dropna().to_numpy() for c in contigs]
        ax.boxplot(data, tick_labels=contigs, showfliers=False)
        for i, c in enumerate(contigs):
            v = mean_tab.loc[c, ex].dropna().to_numpy()
            ax.scatter(np.full(len(v), i + 1), v, s=10, c=HYPO_COLOR, zorder=3)
        ax.set_ylabel("Mean methylation (%)")
        ax.set_title("Per-contig mean methylation across passing samples (red = excluded samples)")
        plt.setp(ax.get_xticklabels(), rotation=90)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(12, 5))
        zc = np.clip(contig_z.T.to_numpy(dtype=float), -10, 10)
        im = ax.imshow(zc, aspect="auto", cmap="RdBu", vmin=-10, vmax=10, interpolation="nearest")
        ax.set_yticks(range(len(autosomes)))
        ax.set_yticklabels(autosomes, fontsize=6)
        ax.set_xlabel("Sample")
        ax.set_title("Robust z of per-contig mean methylation (clipped at +/-10)")
        fig.colorbar(im, ax=ax, fraction=0.02)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


# variable-regions

def cmd_variable_regions(args):
    set_font(args.font_ttf)
    qc = pd.read_csv(args.sample_qc, sep="\t")
    keep_set = set(qc.loc[~qc["exclude"], "sample"])
    if args.contig == "chrY":
        keep_set &= set(qc.loc[qc["inferred_sex"] == "male", "sample"])
    sex = dict(zip(qc["sample"], qc["inferred_sex"]))
    _, pos, samples, X = load_matrix(args.bed)
    cols = [j for j, s in enumerate(samples) if s in keep_set]
    samples = [samples[j] for j in cols]
    X = X[:, cols]
    tested = np.isfinite(X).mean(1) >= args.min_call_rate if X.shape[0] else np.zeros(0, bool)
    Xt, pt = X[tested], pos[tested]
    del X
    tiles = make_tiles(pt, args.tile_cpgs, args.max_tile_span)
    bed_path = f"{args.prefix}.variable_regions.bed"
    calls_path = f"{args.prefix}.variable_region_sample_calls.tsv.gz"
    pdf_path = f"{args.prefix}.variable_regions.pdf"
    bed_cols = ["#chrom", "start", "end", "region_id", "n_tiles", "n_tested_cpgs", "max_tile_sd_z", "mean_beta",
                "sd_sample_means", "min_sample_mean", "max_sample_mean", "n_samples", "n_hyper", "n_hypo",
                "female_minus_male", "hyper_samples", "hypo_samples"]
    call_cols = ["region_id", "sample", "regional_mean", "delta_vs_cpg_median", "n_cpg", "p", "q", "call"]
    tile_cols = ["#chrom", "start", "end", "mean_beta", "sd_across_samples", "bin_median_sd", "bin_mad_sd", "sd_z"]
    if len(tiles) == 0:
        pd.DataFrame(columns=bed_cols).to_csv(bed_path, sep="\t", index=False)
        pd.DataFrame(columns=call_cols).to_csv(calls_path, sep="\t", index=False)
        pd.DataFrame(columns=tile_cols).to_csv(f"{args.prefix}.tile_stats.tsv.gz", sep="\t", index=False)
        write_empty_pdf(pdf_path, f"{args.contig}: no testable tiles ({len(samples)} samples kept)")
        return
    T = Xt[tiles]
    cnt = np.isfinite(T).sum(1)
    with np.errstate(invalid="ignore"):
        tm = np.nanmean(T, 1)
    del T
    tm[cnt < args.tile_min_called] = np.nan
    t_start = pt[tiles[:, 0]]
    t_end = pt[tiles[:, -1]] + 1
    m = np.nanmean(tm, 1)
    sd = np.nanstd(tm, 1, ddof=1)
    nbins = int(np.ceil(100 / args.mean_bin_width))
    b = np.clip((m // args.mean_bin_width).astype(int), 0, nbins - 1)
    bin_med = np.full(nbins, np.nan)
    bin_mad = np.full(nbins, np.nan)
    for k in range(nbins):
        v = sd[(b == k) & np.isfinite(sd)]
        if len(v) >= args.min_tiles_per_bin:
            bin_med[k] = np.median(v)
            bin_mad[k] = 1.4826 * np.median(np.abs(v - bin_med[k]))
    populated = np.where(np.isfinite(bin_med) & (bin_mad > 0))[0]
    for k in range(nbins):
        if not (np.isfinite(bin_med[k]) and bin_mad[k] > 0) and len(populated):
            nearest = populated[np.argmin(np.abs(populated - k))]
            bin_med[k], bin_mad[k] = bin_med[nearest], bin_mad[nearest]
    z = (sd - bin_med[b]) / bin_mad[b]
    pd.DataFrame({"#chrom": args.contig, "start": t_start, "end": t_end, "mean_beta": m, "sd_across_samples": sd,
                  "bin_median_sd": bin_med[b], "bin_mad_sd": bin_mad[b], "sd_z": z}).to_csv(
        f"{args.prefix}.tile_stats.tsv.gz", sep="\t", index=False, float_format="%.3f")
    flagged = np.where(np.isfinite(z) & (z >= args.z_var))[0]
    regions = merge_flagged(t_start[flagged], t_end[flagged], args.merge_gap)

    bed_rows, call_rows, plot_data = [], [], []
    sexes = np.array([sex.get(s, "unknown") for s in samples])
    for ri, (rs, re_, members) in enumerate(regions):
        rid = f"{args.contig}:{rs}-{re_}"
        sel = slice(np.searchsorted(pt, rs), np.searchsorted(pt, re_))
        R = Xt[sel]
        cpg_med = np.nanmedian(R, 1)
        D = R - cpg_med[:, None]
        n_cpg = np.isfinite(R).sum(0)
        with np.errstate(invalid="ignore"):
            reg_mean = np.nanmean(R, 0)
            delta = np.nanmean(D, 0)
        min_n = max(5, int(np.ceil(0.5 * R.shape[0])))
        pvals = np.array([wilcoxon_p(D[:, j]) if n_cpg[j] >= min_n else np.nan for j in range(len(samples))])
        qvals = bh(pvals)
        sig = np.isfinite(qvals) & (qvals < args.fdr)
        hyper = sig & (delta >= args.min_delta)
        hypo = sig & (delta <= -args.min_delta)
        calls = np.where(hyper, "hyper", np.where(hypo, "hypo", "none"))
        f_mean = np.nanmean(reg_mean[sexes == "female"]) if (sexes == "female").any() else np.nan
        m_mean = np.nanmean(reg_mean[sexes == "male"]) if (sexes == "male").any() else np.nan
        tile_z = z[flagged[members]]
        bed_rows.append([args.contig, rs, re_, rid, len(members), int(R.shape[0]), float(np.max(tile_z)),
                         float(np.nanmean(reg_mean)), float(np.nanstd(reg_mean, ddof=1)), float(np.nanmin(reg_mean)),
                         float(np.nanmax(reg_mean)), int(np.isfinite(reg_mean).sum()), int(hyper.sum()),
                         int(hypo.sum()), f_mean - m_mean,
                         ",".join(np.array(samples)[hyper]) or ".", ",".join(np.array(samples)[hypo]) or "."])
        for j, s in enumerate(samples):
            call_rows.append([rid, s, reg_mean[j], delta[j], int(n_cpg[j]), pvals[j], qvals[j], calls[j]])
        plot_data.append((float(np.max(tile_z)), rid, pt[sel], R, calls))
    bed = pd.DataFrame(bed_rows, columns=bed_cols)
    bed.to_csv(bed_path, sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(call_rows, columns=call_cols).to_csv(calls_path, sep="\t", index=False, float_format="%.4g")

    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        ok = np.isfinite(sd) & np.isfinite(m)
        ax.hexbin(m[ok], sd[ok], gridsize=60, bins="log", cmap="Greys", mincnt=1)
        centers = (np.arange(nbins) + 0.5) * args.mean_bin_width
        ax.plot(centers, bin_med, c="black", lw=1, label="bin median SD")
        ax.plot(centers, bin_med + args.z_var * bin_mad, c=HYPO_COLOR, lw=1, label=f"threshold (z = {args.z_var:g})")
        ax.set_xlabel("Tile mean methylation across samples (%)")
        ax.set_ylabel("SD of tile methylation across samples")
        ax.set_title(f"{args.contig}: {len(sd)} tiles, {len(flagged)} flagged, {len(regions)} regions; "
                     f"{len(samples)} samples", fontsize=9)
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        for _, rid, xs, R, calls in sorted(plot_data, key=lambda t: -t[0])[:args.max_plots]:
            fig, ax = plt.subplots(figsize=(8, 4.5))
            line_styles = [("none", BACKGROUND_COLOR, 0.6, 0.35, 1), ("hyper", HYPER_COLOR, 1.1, 0.9, 3),
                           ("hypo", HYPO_COLOR, 1.1, 0.9, 3)]
            for group, color, lw, alpha, zo in line_styles:
                for j in np.where(calls == group)[0]:
                    y = R[:, j]
                    ok = np.isfinite(y)
                    ax.plot(xs[ok], y[ok], c=color, lw=lw, alpha=alpha, zorder=zo, marker="o", ms=1.5)
            ax.set_ylim(-2, 102)
            ax.set_xlabel(f"CpG position on {args.contig} (bp)")
            ax.set_ylabel("Methylation level (%)")
            n_hyper, n_hypo = int((calls == "hyper").sum()), int((calls == "hypo").sum())
            ax.set_title(f"{rid}  ({len(xs)} CpGs)   hyper n={n_hyper} (blue), hypo n={n_hypo} (red), "
                         f"other n={len(calls) - n_hyper - n_hypo} (grey)", fontsize=9)
            ax.ticklabel_format(axis="x", style="plain", useOffset=False)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


# haplotype-stats / haplotype-qc / allele-specific

def haplotype_pairs(names):
    idx = {c: i for i, c in enumerate(names)}
    samples = []
    for c in names:
        if c.endswith("_hap1") and c[:-5] + "_hap2" in idx:
            samples.append(c[:-5])
    return [(s, idx[s + "_hap1"], idx[s + "_hap2"]) for s in samples]


def sample_tile_deltas(pos, h1, h2, tile_cpgs, max_span):
    ok = np.isfinite(h1) & np.isfinite(h2)
    p = pos[ok]
    d = (h1 - h2)[ok]
    tiles = make_tiles(p, tile_cpgs, max_span)
    td = d[tiles].mean(1) if len(tiles) else np.zeros(0)
    return ok, p, d, tiles, td


def cmd_haplotype_stats(args):
    _, pos, names, H = load_matrix(args.bed)
    rows = []
    for s, i1, i2 in haplotype_pairs(names):
        h1, h2 = H[:, i1], H[:, i2]
        ok, p, d, tiles, td = sample_tile_deltas(pos, h1, h2, args.tile_cpgs, args.max_tile_span)
        sigma = 1.4826 * np.median(np.abs(td - np.median(td))) if len(td) >= 50 else np.nan
        with np.errstate(invalid="ignore"):
            rows.append([s, args.contig, H.shape[0], int(np.isfinite(h1).sum()), int(np.isfinite(h2).sum()), int(ok.sum()),
                         np.nanmean(h1) if np.isfinite(h1).any() else np.nan,
                         np.nanmean(h2) if np.isfinite(h2).any() else np.nan,
                         np.median(np.abs(d)) if len(d) else np.nan, len(td), sigma])
    pd.DataFrame(rows, columns=["sample", "contig", "n_cpg_rows", "n_hap1", "n_hap2", "n_both", "mean_hap1", "mean_hap2",
                                "median_abs_delta_cpg", "n_tiles", "tile_delta_sigma"]).to_csv(
        f"{args.prefix}.haplotype_stats.tsv", sep="\t", index=False, na_rep="NA", float_format="%.5g")


def cmd_haplotype_qc(args):
    set_font(args.font_ttf)
    long = pd.concat([pd.read_csv(f, sep="\t", dtype={"contig": str}) for f in args.stats], ignore_index=True)
    sample_qc = pd.read_csv(args.sample_qc, sep="\t").set_index("sample")
    sample_qc["exclude_reasons"] = sample_qc["exclude_reasons"].fillna("")
    auto_re = re.compile(args.autosome_regex)
    contigs = sorted(long["contig"].unique(), key=contig_sort_key)
    samples = list(dict.fromkeys(long["sample"]))
    long.to_csv(f"{args.prefix}.per_contig_haplotype_stats.tsv", sep="\t", index=False, na_rep="NA")
    long.pivot(index="contig", columns="sample", values="n_both").reindex(index=contigs, columns=samples).to_csv(
        f"{args.prefix}.n_both_haplotypes_by_contig.tsv", sep="\t", na_rep="NA")
    auto = long[long["contig"].apply(lambda c: bool(auto_re.match(c)))].copy()
    auto["w_sigma"] = auto["tile_delta_sigma"] * auto["n_tiles"]
    auto["w_h1"] = auto["mean_hap1"] * auto["n_hap1"]
    auto["w_h2"] = auto["mean_hap2"] * auto["n_hap2"]
    g = auto.groupby("sample", sort=False)
    qc = pd.DataFrame(index=pd.Index(samples, name="sample"))
    qc["n_both_autosomal"] = g["n_both"].sum()
    qc["frac_both_vs_cohort_median"] = qc["n_both_autosomal"] / qc["n_both_autosomal"].median()
    qc["tile_delta_sigma_autosomal"] = g["w_sigma"].sum() / g["n_tiles"].sum()
    qc["mean_hap1_autosomal"] = g["w_h1"].sum() / g["n_hap1"].sum()
    qc["mean_hap2_autosomal"] = g["w_h2"].sum() / g["n_hap2"].sum()
    qc["hap1_minus_hap2_mean"] = qc["mean_hap1_autosomal"] - qc["mean_hap2_autosomal"]
    qc["z_sigma"] = robust_z(qc["tile_delta_sigma_autosomal"])
    qc["z_hap_mean_diff"] = robust_z(qc["hap1_minus_hap2_mean"])
    qc["sample_level_exclude_reasons"] = sample_qc.reindex(samples)["exclude_reasons"].fillna("not_in_sample_qc").to_numpy()
    reasons = []
    for s, r in qc.iterrows():
        why = [f"sample_qc:{r['sample_level_exclude_reasons']}"] if isinstance(r["sample_level_exclude_reasons"], str) \
            and r["sample_level_exclude_reasons"] else []
        if r["frac_both_vs_cohort_median"] < args.min_both_frac:
            why.append("few_cpgs_with_both_haplotypes")
        if r["z_sigma"] >= args.z:
            why.append("noisy_haplotype_delta")
        if abs(r["z_hap_mean_diff"]) >= args.z:
            why.append("haplotype_mean_imbalance")
        reasons.append(";".join(why))
    qc["exclude_reasons"] = reasons
    qc["exclude"] = qc["exclude_reasons"] != ""
    qc["inferred_sex"] = sample_qc.reindex(samples)["inferred_sex"].fillna("unknown").to_numpy()
    qc.reset_index().to_csv(f"{args.prefix}.haplotype_sample_qc.tsv", sep="\t", index=False, na_rep="NA",
                            float_format="%.5g")
    with open(f"{args.prefix}.haplotype_excluded_samples.txt", "w") as fh:
        for s in qc.index[qc["exclude"]]:
            fh.write(f"{s}\t{qc.loc[s, 'exclude_reasons']}\n")
    ex = qc["exclude"].to_numpy()
    with PdfPages(f"{args.prefix}.haplotype_sample_qc.pdf") as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        axes[0].scatter(qc["frac_both_vs_cohort_median"][~ex], qc["tile_delta_sigma_autosomal"][~ex], s=10, c=BACKGROUND_COLOR)
        axes[0].scatter(qc["frac_both_vs_cohort_median"][ex], qc["tile_delta_sigma_autosomal"][ex], s=12, c=HYPO_COLOR)
        axes[0].set_xlabel("CpGs with both haplotypes (fraction of cohort median)")
        axes[0].set_ylabel("Noise of hap1-hap2 tile delta (robust SD)")
        axes[1].scatter(qc["mean_hap1_autosomal"][~ex], qc["mean_hap2_autosomal"][~ex], s=10, c=BACKGROUND_COLOR)
        axes[1].scatter(qc["mean_hap1_autosomal"][ex], qc["mean_hap2_autosomal"][ex], s=12, c=HYPO_COLOR)
        axes[1].set_xlabel("Autosomal mean, hap1 (%)")
        axes[1].set_ylabel("Autosomal mean, hap2 (%)")
        fig.suptitle(f"Haplotype QC ({ex.sum()} of {len(ex)} excluded, red)")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def cmd_allele_specific(args):
    set_font(args.font_ttf)
    qc = pd.read_csv(args.haplotype_qc, sep="\t")
    keep_set = set(qc.loc[~qc["exclude"], "sample"])
    _, pos, names, H = load_matrix(args.bed)
    pairs = [p for p in haplotype_pairs(names) if p[0] in keep_set]
    reg_cols = ["#chrom", "start", "end", "sample", "n_tiles", "n_cpg", "mean_hap1", "mean_hap2", "delta_hap1_minus_hap2",
                "higher_haplotype", "tile_delta_sigma", "p", "q"]
    rows = []
    informative = {}
    for s, i1, i2 in pairs:
        h1, h2 = H[:, i1], H[:, i2]
        ok, p, d, tiles, td = sample_tile_deltas(pos, h1, h2, args.tile_cpgs, args.max_tile_span)
        informative[s] = p
        if len(td) < 50:
            continue
        ok_rows = np.where(ok)[0]
        sigma = 1.4826 * np.median(np.abs(td - np.median(td)))
        thr = max(args.min_delta, args.sigma_k * sigma)
        fl = np.where(np.abs(td) >= thr)[0]
        if len(fl) == 0:
            continue
        ts = p[tiles[fl, 0]]
        te = p[tiles[fl, -1]] + 1
        for rs, re_, members in merge_flagged(ts, te, args.merge_gap, signs=np.sign(td[fl])):
            if len(members) < args.min_tiles:
                continue
            a, b = np.searchsorted(p, rs), np.searchsorted(p, re_)
            dd = d[a:b]
            rows_in = ok_rows[a:b]
            m1, m2 = float(np.mean(h1[rows_in])), float(np.mean(h2[rows_in]))
            rows.append([args.contig, rs, re_, s, len(members), b - a, m1, m2, m1 - m2, "hap1" if m1 > m2 else "hap2",
                         sigma, wilcoxon_p(dd), np.nan])
    asm = pd.DataFrame(rows, columns=reg_cols)
    asm["q"] = bh(asm["p"].to_numpy()) if len(asm) else []
    asm = asm[asm["q"] < args.fdr].sort_values(["start", "sample"]).reset_index(drop=True)
    asm.to_csv(f"{args.prefix}.asm_regions_per_sample.tsv.gz", sep="\t", index=False, float_format="%.4g")

    pop_cols = ["#chrom", "start", "end", "region_id", "n_asm_samples", "n_informative_samples", "frac_asm_samples",
                "median_abs_delta", "asm_samples"]
    pop_rows, plot_data = [], []
    if len(asm):
        order = np.argsort(asm["start"].to_numpy(), kind="stable")
        st, en = asm["start"].to_numpy()[order], asm["end"].to_numpy()[order]
        groups, cur_e = [], -1
        for k in range(len(st)):
            if groups and st[k] < cur_e:
                groups[-1].append(order[k])
                cur_e = max(cur_e, en[k])
            else:
                groups.append([order[k]])
                cur_e = en[k]
        for gi in groups:
            sub = asm.iloc[gi]
            rs, re_ = int(sub["start"].min()), int(sub["end"].max())
            ss = sorted(set(sub["sample"]))
            n_inf = sum(1 for s, p in informative.items() if np.searchsorted(p, re_) - np.searchsorted(p, rs) >= args.tile_cpgs)
            rid = f"{args.contig}:{rs}-{re_}"
            pop_rows.append([args.contig, rs, re_, rid, len(ss), n_inf, len(ss) / n_inf if n_inf else np.nan,
                             float(np.median(np.abs(sub["delta_hap1_minus_hap2"]))), ",".join(ss)])
            plot_data.append((len(ss), rid, rs, re_, sub))
    pop = pd.DataFrame(pop_rows, columns=pop_cols)
    pop.to_csv(f"{args.prefix}.asm_regions_population.bed", sep="\t", index=False, float_format="%.4g")

    pdf_path = f"{args.prefix}.asm_regions.pdf"
    if not plot_data:
        write_empty_pdf(pdf_path, f"{args.contig}: no ASM regions ({len(pairs)} samples kept)")
        return
    col = {s: (i1, i2) for s, i1, i2 in pairs}
    with PdfPages(pdf_path) as pdf:
        for n_s, rid, rs, re_, sub in sorted(plot_data, key=lambda t: -t[0])[:args.max_plots]:
            sel = (pos >= rs) & (pos < re_)
            xs = pos[sel]
            fig, ax = plt.subplots(figsize=(8, 4.5))
            for s in sorted(set(sub["sample"])):
                i1, i2 = col[s]
                a, b = H[sel, i1], H[sel, i2]
                hi, lo = (a, b) if np.nanmean(a) > np.nanmean(b) else (b, a)
                for y, c in [(hi, HYPER_COLOR), (lo, HYPO_COLOR)]:
                    okp = np.isfinite(y)
                    ax.plot(xs[okp], y[okp], c=c, lw=0.8, alpha=0.6, marker="o", ms=1.2)
            ax.set_ylim(-2, 102)
            ax.set_xlabel(f"CpG position on {args.contig} (bp)")
            ax.set_ylabel("Methylation level (%)")
            n_inf = int(pop.loc[pop["region_id"] == rid, "n_informative_samples"].iloc[0])
            ax.set_title(f"{rid}  ASM in {n_s}/{n_inf} informative samples; higher allele blue, lower allele red",
                         fontsize=9)
            ax.ticklabel_format(axis="x", style="plain", useOffset=False)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


# concat

def cmd_concat(args):
    opener = gzip.open if args.out.endswith(".gz") else open
    header_written = False
    with opener(args.out, "wt") as out:
        for f in args.inputs:
            fo = gzip.open if f.endswith(".gz") else open
            with fo(f, "rt") as fh:
                header = fh.readline()
                if not header_written:
                    out.write(header)
                    header_written = True
                for line in fh:
                    out.write(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("contig-stats")
    p.add_argument("--bed", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--min-call-rate-for-corr", type=float, default=0.5)
    p.set_defaults(func=cmd_contig_stats)

    p = sub.add_parser("sample-qc")
    p.add_argument("--stats", nargs="+", required=True)
    p.add_argument("--hists", nargs="+", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--font-ttf", required=True)
    p.add_argument("--z", type=float, default=5.0)
    p.add_argument("--min-call-rate", type=float, default=0.5)
    p.add_argument("--min-outlier-autosomes-warn", type=int, default=3)
    p.add_argument("--male-chrY-frac", type=float, default=0.25)
    p.add_argument("--autosome-regex", default=r"^chr[0-9]+$")
    p.set_defaults(func=cmd_sample_qc)

    p = sub.add_parser("variable-regions")
    p.add_argument("--bed", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--sample-qc", required=True)
    p.add_argument("--font-ttf", required=True)
    p.add_argument("--min-call-rate", type=float, default=0.8)
    p.add_argument("--tile-cpgs", type=int, default=10)
    p.add_argument("--tile-min-called", type=int, default=7)
    p.add_argument("--max-tile-span", type=int, default=2000)
    p.add_argument("--mean-bin-width", type=float, default=5.0)
    p.add_argument("--min-tiles-per-bin", type=int, default=50)
    p.add_argument("--z-var", type=float, default=4.0)
    p.add_argument("--merge-gap", type=int, default=1000)
    p.add_argument("--min-delta", type=float, default=20.0)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--max-plots", type=int, default=1000)
    p.set_defaults(func=cmd_variable_regions)

    p = sub.add_parser("haplotype-stats")
    p.add_argument("--bed", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--tile-cpgs", type=int, default=10)
    p.add_argument("--max-tile-span", type=int, default=2000)
    p.set_defaults(func=cmd_haplotype_stats)

    p = sub.add_parser("haplotype-qc")
    p.add_argument("--stats", nargs="+", required=True)
    p.add_argument("--sample-qc", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--font-ttf", required=True)
    p.add_argument("--z", type=float, default=5.0)
    p.add_argument("--min-both-frac", type=float, default=0.25)
    p.add_argument("--autosome-regex", default=r"^chr[0-9]+$")
    p.set_defaults(func=cmd_haplotype_qc)

    p = sub.add_parser("allele-specific")
    p.add_argument("--bed", required=True)
    p.add_argument("--contig", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--haplotype-qc", required=True)
    p.add_argument("--font-ttf", required=True)
    p.add_argument("--tile-cpgs", type=int, default=10)
    p.add_argument("--max-tile-span", type=int, default=2000)
    p.add_argument("--min-delta", type=float, default=20.0)
    p.add_argument("--sigma-k", type=float, default=4.0)
    p.add_argument("--merge-gap", type=int, default=1000)
    p.add_argument("--min-tiles", type=int, default=1)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--max-plots", type=int, default=200)
    p.set_defaults(func=cmd_allele_specific)

    p = sub.add_parser("concat")
    p.add_argument("--inputs", nargs="+", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_concat)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
