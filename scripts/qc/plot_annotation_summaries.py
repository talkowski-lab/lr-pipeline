#!/usr/bin/env python3

import argparse
import os
import re

import matplotlib
import matplotlib.colors as mc
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from matplotlib.gridspec import GridSpec
from scipy.ndimage import uniform_filter1d
from scipy.stats import chi2 as chi2_dist

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Fall back to the metric-compatible Liberation Sans where Arial is not installed
matplotlib.rcParams["font.family"] = ["Arial", "Liberation Sans", "DejaVu Sans"]

# Font size constants
FTITLE = 14  # main figure title
FSUB = 10  # n=X subtitle
FSUBTITLE = 12  # subplot title
FLABEL = 10  # axis labels
FTICK = 8  # axis tick labels
FHEAT = 10  # heatmap cell text

# Font size constants for the hwe_ plots specifically (larger, since these are dense grids of small panels)
HWE_FSUBTITLE = FSUBTITLE + 3
HWE_FSUB = FSUB + 2
HWE_FLABEL = FLABEL + 3
HWE_FTICK = FTICK + 3
HWE_FLEGEND = 12

GROUPED_TYPE_ORDER = ["SNV", "DEL_SHORT", "INS_SHORT", "DUP_SHORT", "DEL_SV", "INS_SV", "DUP_SV", "TRV", "OTH"]

# HWE type panel order — TR broken into subtypes. Grouped into columns (SHORT/TR_INS row above
# SV/TR_DEL row) to match QcAnnotations' gt_type_distribution layout, e.g. INS_SV below INS_SHORT.
HWE_TYPE_ORDER = ["SNV", "INS_SHORT", "DEL_SHORT", "DUP_SHORT", "TR_INS", "TR_SNV", "INS_SV", "DEL_SV", "DUP_SV", "TR_DEL", "OTH"]

GROUPED_COLORS = {
    "SNV": "#7459B2",
    "DEL_SHORT": "#D43925",
    "DEL_SV": "#8B1A11",
    "INS_SHORT": "#D474E0",
    "INS_SV": "#9B2DB4",
    "DUP_SHORT": "#2376B2",
    "DUP_SV": "#1A5B8A",
    "TRV": "#FA931E",
    "TR_SNV": "#FA931E",
    "TR_INS": "#FFBC5E",
    "TR_DEL": "#B56A00",
    "OTH": "#397246",
}

# Base type → SHORT-equivalent color for stacked plots
STACK_COLORS = {
    "SNV": GROUPED_COLORS["SNV"],
    "DEL": GROUPED_COLORS["DEL_SHORT"],
    "INS": GROUPED_COLORS["INS_SHORT"],
    "DUP": GROUPED_COLORS["DUP_SHORT"],
    "TRV": GROUPED_COLORS["TRV"],
    "Other": GROUPED_COLORS["OTH"],
}

META_COLS = {"category", "sub_category", "tr_status", "region"}
REGION_ORDER = ["US", "RM", "SD", "SR"]


def format_count(n):
    if abs(n) >= 1e6:
        return f"{n / 1e6:.1f}M"
    if abs(n) >= 1e3:
        return f"{n / 1e3:.1f}K"
    return str(int(round(n)))


def lighten_color(color, factor=0.45):
    c = mc.to_rgb(color)
    return tuple(ci + factor * (1 - ci) for ci in c)


def get_variant_buckets(sites_df):
    return [c for c in sites_df.columns if c not in META_COLS and "[" not in c and not c.endswith(" (Matched)")]


def _get_af_labels(df):
    seen, result = set(), []
    for c in df.columns:
        m = re.search(r"\[(.+)\]", c)
        if m:
            lab = m.group(1)
            if lab not in seen:
                result.append(lab)
                seen.add(lab)
    return result


def format_af_label(label):
    """Convert raw AF label like '0.0-0.01' or '>0.5' to percentage string."""
    if label.startswith(">"):
        val = float(label[1:]) * 100
        return f">{val:.3g}%"
    lo_str, hi_str = label.split("-", 1)
    lo, hi = float(lo_str), float(hi_str)
    if lo == 0.0:
        return f"≤{hi * 100:.3g}%"
    return f"{lo * 100:.3g}-{hi * 100:.3g}%"


def _af_bin_mask(af_values, label):
    """Boolean mask for a raw AF bucket label (matches the af_bins_plotting convention:
    first bin inclusive-inclusive, middle bins exclusive-inclusive, last bin open)."""
    if label.startswith(">"):
        return af_values > float(label[1:])
    lo, hi = (float(v) for v in label.split("-", 1))
    if lo == 0.0:
        return (af_values >= lo) & (af_values <= hi)
    return (af_values > lo) & (af_values <= hi)


def _title_with_subtitle(ax, title, subtitle=None, subtitle_offset=0, title_fontsize=FSUBTITLE, subtitle_fontsize=None):
    """Subplot title with an optional smaller 'n=X'-style subtitle below it, spaced to avoid overlap.
    subtitle_offset raises the subtitle (and title) above a panel's twiny MEI marker labels."""
    ax.set_title(title, fontsize=title_fontsize, fontweight="bold", pad=18 + subtitle_offset)
    if subtitle:
        ax.annotate(
            subtitle,
            xy=(0.5, 1.0),
            xycoords="axes fraction",
            xytext=(0, subtitle_offset),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=subtitle_fontsize if subtitle_fontsize is not None else FSUB - 1,
            color="#555555",
            clip_on=False,
        )


def get_size_bins(buckets):
    seen, result = set(), []
    for b in buckets:
        parts = b.split(" ", 1)
        if len(parts) == 2 and parts[0] in ("DEL", "INS", "DUP") and parts[1] not in seen:
            result.append(parts[1])
            seen.add(parts[1])
    return result


def get_total_row(sites_df):
    mask = (sites_df["category"] == "All") & (sites_df["sub_category"] == "All") & (sites_df["tr_status"] == "All")
    if "region" in sites_df.columns:
        mask &= sites_df["region"] == "All"
    return sites_df[mask].iloc[0]


def parse_lower_bound(size_label):
    if size_label.endswith("+"):
        return int(size_label[:-1])
    if "-" in size_label:
        return int(size_label.split("-")[0])
    return 0


def _bp_parts(n):
    if n < 1000:
        return str(n), "bp"
    if n < 1_000_000:
        v = n / 1000
        return (f"{int(v)}" if v == int(v) else f"{v:.1f}"), "kbp"
    v = n / 1_000_000
    return (f"{int(v)}" if v == int(v) else f"{v:.1f}"), "mbp"


def format_size_label(size_label):
    if size_label.endswith("+"):
        val, unit = _bp_parts(int(size_label[:-1]))
        return f"≥{val}{unit}"
    lo_str, hi_str = size_label.split("-")
    lo, hi = int(lo_str), int(hi_str)
    if lo == hi:
        val, unit = _bp_parts(lo)
        return f"{val}{unit}"
    if str(hi).endswith("9"):
        hi += 1
    lo_val, lo_unit = _bp_parts(lo)
    hi_val, hi_unit = _bp_parts(hi)
    if lo_unit == hi_unit:
        return f"{lo_val}-{hi_val}{hi_unit}"
    return f"{lo_val}{lo_unit}-{hi_val}{hi_unit}"


def bucket_to_grouped_type(bucket, manuscript_mode=False):
    b = bucket.strip()
    if b == "SNV":
        return "SNV"
    if b == "TRV":
        return "TRV"
    if b in ("Other", "OTH"):
        return "OTH"
    parts = b.split(" ", 1)
    if len(parts) != 2 or parts[0] not in ("DEL", "INS", "DUP"):
        return "OTH"
    typ, size_label = parts
    if manuscript_mode and typ == "DUP":
        typ = "INS"
    return f"{typ}_SV" if parse_lower_bound(size_label) >= 50 else f"{typ}_SHORT"


def _group_counts(active_buckets, counts_series, manuscript_mode=False):
    """Aggregate per-bucket counts into grouped variant types."""
    grouped = {}
    for b in active_buckets:
        gt = bucket_to_grouped_type(b, manuscript_mode)
        grouped[gt] = grouped.get(gt, 0.0) + float(counts_series.get(b, 0))
    return grouped


# ==============================
# Plot 1: counts_overall_distribution.png
# ==============================
def plot_overall_distribution(sites_df, outpath, manuscript_mode=False):
    buckets = get_variant_buckets(sites_df)
    total_row = get_total_row(sites_df)
    counts = total_row[buckets].astype(float)
    active_buckets = [b for b in buckets if counts[b] > 0]
    active_counts = counts[active_buckets]
    size_bins = [sz for sz in get_size_bins(active_buckets) if any(counts.get(f"{t} {sz}", 0) > 0 for t in ("DEL", "INS", "DUP"))]

    matched_counts = pd.Series({b: float(total_row.get(f"{b} (Matched)", 0)) for b in active_buckets})

    has_region = "region" in sites_df.columns and sites_df["region"].isin(REGION_ORDER).any()
    af_labels = _get_af_labels(sites_df)
    has_af = len(af_labels) > 0

    n_sites = format_count(int(counts[active_buckets].sum()))
    ncols = 2 + int(has_af) + int(has_region)
    fig, axes = plt.subplots(2, ncols, figsize=(7 * ncols, 9))
    fig.patch.set_facecolor("white")

    _plot_variant_count(
        axes[0, 0],
        active_buckets,
        active_counts,
        matched_counts,
        title="Site Counts by Type",
        subtitle=f"N = {n_sites} sites",
        manuscript_mode=manuscript_mode,
    )
    size_slot = axes[0, 1].get_subplotspec()
    axes[0, 1].remove()
    _plot_stacked_by_size_broken(fig, size_slot, counts, size_bins, manuscript_mode=manuscript_mode)
    _plot_variant_proportion(axes[1, 0], active_buckets, active_counts, manuscript_mode=manuscript_mode)
    _plot_stacked_by_size(axes[1, 1], counts, size_bins, scaled=True, show_legend=False, manuscript_mode=manuscript_mode)
    axes[1, 1].set_title("Site Counts by Size (Scaled)", fontsize=FSUBTITLE, fontweight="bold")

    col = 2
    if has_af:
        _plot_stacked_by_af(axes[0, col], sites_df, af_labels, scaled=False, show_legend=False, manuscript_mode=manuscript_mode)
        axes[0, col].set_title("Site Counts by AF", fontsize=FSUBTITLE, fontweight="bold")
        _plot_stacked_by_af(axes[1, col], sites_df, af_labels, scaled=True, show_legend=False, manuscript_mode=manuscript_mode)
        axes[1, col].set_title("Site Counts by AF (Scaled)", fontsize=FSUBTITLE, fontweight="bold")
        col += 1

    if has_region:
        _plot_stacked_by_region(axes[0, col], sites_df, size_bins, scaled=False, manuscript_mode=manuscript_mode)
        axes[0, col].set_title("Site Counts by Region", fontsize=FSUBTITLE, fontweight="bold")
        _plot_stacked_by_region(axes[1, col], sites_df, size_bins, scaled=True, manuscript_mode=manuscript_mode)
        axes[1, col].set_title("Site Counts by Region (Scaled)", fontsize=FSUBTITLE, fontweight="bold")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def _plot_variant_count(
    ax, active_buckets, active_counts, matched_counts, title="Site Count", subtitle=None, manuscript_mode=False
):
    grouped = _group_counts(active_buckets, active_counts, manuscript_mode)
    grouped_matched = _group_counts(active_buckets, matched_counts, manuscript_mode)

    gtypes = [gt for gt in GROUPED_TYPE_ORDER if gt in grouped]
    total_vals = np.array([grouped[gt] for gt in gtypes])
    matched_vals = np.array([grouped_matched.get(gt, 0.0) for gt in gtypes])
    unmatched_vals = total_vals - matched_vals
    colors = [GROUPED_COLORS.get(gt, "#888") for gt in gtypes]

    x = np.arange(len(gtypes))
    lighter_colors = [lighten_color(c) for c in colors]
    ax.bar(x, matched_vals, color=colors, width=0.7, edgecolor="white", linewidth=0.5)
    ax.bar(x, unmatched_vals, bottom=matched_vals, color=lighter_colors, width=0.7, edgecolor="white", linewidth=0.5)

    ax.set_yscale("log")
    ax.set_ylim(top=max(total_vals) * 30)
    ax.set_xticks(x)
    ax.set_xticklabels(gtypes, rotation=45, ha="right", fontsize=FTICK)
    ax.set_xlabel("Type", fontsize=FLABEL)
    ax.set_ylabel("Count", fontsize=FLABEL)
    _title_with_subtitle(ax, title, subtitle)

    for i, (gt, tot, matched) in enumerate(zip(gtypes, total_vals, matched_vals)):
        label = format_count(int(tot))
        if matched > 0:
            label += f"\n({format_count(int(matched))})"
        ax.text(i, tot * 1.4, label, ha="center", va="bottom", fontsize=6, color=GROUPED_COLORS.get(gt, "#888"), linespacing=1.1)

    ax.spines[["top", "right"]].set_visible(False)


def _plot_variant_proportion(ax, active_buckets, active_counts, manuscript_mode=False):
    grouped = _group_counts(active_buckets, active_counts, manuscript_mode)

    gtypes = [gt for gt in GROUPED_TYPE_ORDER if gt in grouped]
    total = sum(grouped.values())
    props = [grouped[gt] / total * 100 for gt in gtypes]
    colors = [GROUPED_COLORS.get(gt, "#888") for gt in gtypes]

    ax.bar(range(len(gtypes)), props, color=colors, width=0.7, edgecolor="white", linewidth=0.5)
    ax.set_xticks(range(len(gtypes)))
    ax.set_xticklabels(gtypes, rotation=45, ha="right", fontsize=FTICK)
    ax.set_xlabel("Type", fontsize=FLABEL)
    ax.set_ylabel("Proportion (%)", fontsize=FLABEL)
    ax.set_title("Site Counts by Type (Scaled)", fontsize=FSUBTITLE, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)


def _annotate_bar_totals(ax, x, totals):
    """Black count labels above each bar in an unscaled stacked bar plot, with headroom to avoid the title."""
    if len(totals) == 0 or max(totals) <= 0:
        return
    top = max(totals)
    ax.set_ylim(top=top * 1.2)
    for xi, ti in zip(x, totals):
        if ti > 0:
            ax.text(xi, ti + top * 0.02, format_count(ti), ha="center", va="bottom", fontsize=6, color="black")


def _plot_stacked_by_size(ax, counts, size_bins, scaled=False, show_legend=True, manuscript_mode=False):
    type_order = ["DEL", "INS", "DUP"]

    mat = {typ: [] for typ in type_order}
    for sz_lab in size_bins:
        for typ in type_order:
            mat[typ].append(float(counts.get(f"{typ} {sz_lab}", 0)))

    if manuscript_mode:
        mat["INS"] = [i + d for i, d in zip(mat["INS"], mat["DUP"])]
        type_order = ["DEL", "INS"]

    data = np.array([mat[t] for t in type_order], dtype=float)
    if scaled:
        col_sums = np.where(data.sum(axis=0) > 0, data.sum(axis=0), 1)
        data = data / col_sums * 100

    x = np.arange(len(size_bins))
    bottom = np.zeros(len(size_bins))
    handles = []
    for i, typ in enumerate(type_order):
        b = ax.bar(x, data[i], bottom=bottom, color=STACK_COLORS[typ], label=typ, width=0.7, edgecolor="white", linewidth=0.3)
        handles.append(b[0])
        bottom += data[i]

    ax.set_ylabel("Proportion (%)" if scaled else "Count", fontsize=FLABEL)

    x_labels = [format_size_label(s) for s in size_bins]
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=FTICK)
    ax.set_xlabel("Size", fontsize=FLABEL)
    if show_legend:
        ax.legend(handles=handles, labels=type_order, fontsize=7, loc="upper right", framealpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def _plot_stacked_by_size_broken(fig, ax_slot, counts, size_bins, manuscript_mode=False):
    """Count by Size with a broken y-axis: top axis scaled to the first size bucket,
    bottom axis scaled to the remaining buckets, so a dominant first bucket doesn't
    flatten the rest of the bars."""
    type_order = ["DEL", "INS", "DUP"]

    mat = {typ: [] for typ in type_order}
    for sz_lab in size_bins:
        for typ in type_order:
            mat[typ].append(float(counts.get(f"{typ} {sz_lab}", 0)))
    if manuscript_mode:
        mat["INS"] = [i + d for i, d in zip(mat["INS"], mat["DUP"])]
        type_order = ["DEL", "INS"]
    data = np.array([mat[t] for t in type_order], dtype=float)
    totals = data.sum(axis=0)

    x = np.arange(len(size_bins))
    x_labels = [format_size_label(s) for s in size_bins]

    if len(size_bins) < 2 or totals[0] <= max(totals[1:]):
        ax = fig.add_subplot(ax_slot)
        bottom = np.zeros(len(size_bins))
        for i, typ in enumerate(type_order):
            ax.bar(x, data[i], bottom=bottom, color=STACK_COLORS[typ], label=typ, width=0.7, edgecolor="white", linewidth=0.3)
            bottom += data[i]
        _annotate_bar_totals(ax, x, totals)
        ax.set_ylabel("Count", fontsize=FLABEL)
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=FTICK)
        ax.set_xlabel("Size", fontsize=FLABEL)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title("Site Counts by Size", fontsize=FSUBTITLE, fontweight="bold")
        return

    rest_max = max(totals[1:])
    gs = ax_slot.subgridspec(2, 1, height_ratios=[1, 3], hspace=0.08)
    ax_hi = fig.add_subplot(gs[0])
    ax_lo = fig.add_subplot(gs[1], sharex=ax_hi)

    for ax in (ax_hi, ax_lo):
        bottom = np.zeros(len(size_bins))
        for i, typ in enumerate(type_order):
            ax.bar(x, data[i], bottom=bottom, color=STACK_COLORS[typ], width=0.7, edgecolor="white", linewidth=0.3)
            bottom += data[i]

    ax_lo.set_ylim(0, rest_max * 1.25)
    ax_hi.set_ylim(rest_max * 1.15, totals[0] * 1.1)

    for xi, ti in zip(x[1:], totals[1:]):
        if ti > 0:
            ax_lo.text(xi, ti + rest_max * 0.03, format_count(ti), ha="center", va="bottom", fontsize=6, color="black")
    if totals[0] > 0:
        ax_hi.text(
            x[0], totals[0] + totals[0] * 0.01, format_count(totals[0]), ha="center", va="bottom", fontsize=6, color="black"
        )

    ax_hi.spines["bottom"].set_visible(False)
    ax_hi.spines[["top", "right"]].set_visible(False)
    ax_hi.tick_params(bottom=False, labelbottom=False)
    ax_lo.spines["top"].set_visible(False)
    ax_lo.spines["right"].set_visible(False)

    d = 0.5
    break_kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10, linestyle="none", color="k", mec="k", mew=1, clip_on=False)
    ax_hi.plot([0], [0], transform=ax_hi.transAxes, **break_kwargs)
    ax_lo.plot([0], [1], transform=ax_lo.transAxes, **break_kwargs)

    ax_lo.set_xticks(x)
    ax_lo.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=FTICK)
    ax_lo.set_xlabel("Size", fontsize=FLABEL)
    ax_lo.set_ylabel("Count", fontsize=FLABEL)
    ax_hi.set_title("Site Counts by Size", fontsize=FSUBTITLE, fontweight="bold", pad=10)


def _plot_stacked_by_region(ax, sites_df, size_bins, scaled=False, manuscript_mode=False):
    type_order = ["SNV", "DEL", "INS", "DUP", "TRV", "Other"]

    def type_count_from_series(series, typ):
        if typ in ("SNV", "TRV", "Other"):
            return float(series.get(typ, 0))
        return sum(float(series.get(f"{typ} {s}", 0)) for s in size_bins)

    region_rows = {}
    for region in REGION_ORDER:
        mask = (
            (sites_df["category"] == "All")
            & (sites_df["sub_category"] == "All")
            & (sites_df["tr_status"] == "All")
            & (sites_df["region"] == region)
        )
        rows = sites_df[mask]
        if not rows.empty:
            region_rows[region] = rows.iloc[0]

    active_regions = [r for r in REGION_ORDER if r in region_rows]

    mat = {typ: [] for typ in type_order}
    for region in active_regions:
        for typ in type_order:
            mat[typ].append(type_count_from_series(region_rows[region], typ))

    if manuscript_mode:
        mat["INS"] = [i + d for i, d in zip(mat["INS"], mat["DUP"])]
        type_order = ["SNV", "DEL", "INS", "TRV", "Other"]

    data = np.array([mat[t] for t in type_order], dtype=float)
    if scaled:
        col_sums = np.where(data.sum(axis=0) > 0, data.sum(axis=0), 1)
        data = data / col_sums * 100

    x = np.arange(len(active_regions))
    bottom = np.zeros(len(active_regions))
    for i, typ in enumerate(type_order):
        ax.bar(x, data[i], bottom=bottom, color=STACK_COLORS[typ], label=typ, width=0.7, edgecolor="white", linewidth=0.3)
        bottom += data[i]

    if not scaled:
        _annotate_bar_totals(ax, x, data.sum(axis=0))

    ax.set_ylabel("Proportion (%)" if scaled else "Count", fontsize=FLABEL)

    ax.set_xticks(x)
    ax.set_xticklabels(active_regions, rotation=45, ha="right", fontsize=FTICK)
    ax.set_xlabel("Region", fontsize=FLABEL)
    ax.spines[["top", "right"]].set_visible(False)


def _plot_stacked_by_af(ax, sites_df, af_labels, scaled=False, show_legend=True, manuscript_mode=False):
    type_order = ["SNV", "DEL", "INS", "DUP", "TRV", "Other"]
    total_row = get_total_row(sites_df)

    mat = {typ: [] for typ in type_order}
    for label in af_labels:
        for typ in type_order:
            mat[typ].append(float(total_row.get(f"{typ} [{label}]", 0)))

    if manuscript_mode:
        mat["INS"] = [i + d for i, d in zip(mat["INS"], mat["DUP"])]
        type_order = ["SNV", "DEL", "INS", "TRV", "Other"]

    data = np.array([mat[t] for t in type_order], dtype=float)
    if scaled:
        col_sums = np.where(data.sum(axis=0) > 0, data.sum(axis=0), 1)
        data = data / col_sums * 100

    x = np.arange(len(af_labels))
    bottom = np.zeros(len(af_labels))
    handles = []
    for i, typ in enumerate(type_order):
        b = ax.bar(x, data[i], bottom=bottom, color=STACK_COLORS[typ], label=typ, width=0.7, edgecolor="white", linewidth=0.3)
        handles.append(b[0])
        bottom += data[i]

    if not scaled:
        _annotate_bar_totals(ax, x, data.sum(axis=0))

    ax.set_ylabel("Proportion (%)" if scaled else "Count", fontsize=FLABEL)

    ax.set_xticks(x)
    ax.set_xticklabels([format_af_label(label) for label in af_labels], rotation=45, ha="right", fontsize=FTICK)
    ax.set_xlabel("AF", fontsize=FLABEL)
    if show_legend:
        ax.legend(handles=handles, labels=type_order, fontsize=7, loc="upper right", framealpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)


# ==============================
# Plot 1b: srwgs_overlap.png
# ==============================
OVERLAP_SIZE_CLASSES = [("1-49", 1, 50), ("50-499", 50, 500), ("500+", 500, None)]


def plot_srwgs_overlap(sites_df, outpath, manuscript_mode=False):
    """Site counts for SNVs and DEL/INS by size class (log scale), with the srWGS-matched portion
    of each bar hatched and labeled with its percentage of the bar."""
    total_row = get_total_row(sites_df)

    class_buckets = {}
    for b in get_variant_buckets(sites_df):
        parts = b.split(" ", 1)
        if len(parts) != 2 or parts[0] not in ("DEL", "INS", "DUP"):
            continue
        typ = bucket_to_grouped_type(b, manuscript_mode).split("_")[0]
        lb = parse_lower_bound(parts[1])
        for size_class, lo, hi in OVERLAP_SIZE_CLASSES:
            if lb >= lo and (hi is None or lb < hi):
                class_buckets.setdefault((size_class, typ), []).append(b)

    def _counts(cols):
        return (sum(float(total_row.get(c, 0)) for c in cols), sum(float(total_row.get(f"{c} (Matched)", 0)) for c in cols))

    # (x, label, color, total, matched); size classes are separated by a gap
    bars = [(0.0, "SNV", GROUPED_COLORS["SNV"], *_counts(["SNV"]))]
    groups = []
    pos = 0.0
    for size_class, lo, _ in OVERLAP_SIZE_CLASSES:
        pos += 0.5
        start = pos + 1
        for typ in ("DEL", "INS", "DUP"):
            tot, matched = _counts(class_buckets.get((size_class, typ), []))
            if tot > 0:
                pos += 1
                color = GROUPED_COLORS[f"{typ}_SHORT" if lo < 50 else f"{typ}_SV"]
                bars.append((pos, typ, color, tot, matched))
        if pos >= start:
            groups.append((format_size_label(size_class), start, pos))

    x = np.array([b[0] for b in bars])
    colors = [b[2] for b in bars]
    total_vals = np.array([b[3] for b in bars])
    matched_vals = np.array([b[4] for b in bars])

    fig, ax = plt.subplots(figsize=(1 + 0.9 * len(bars), 5))
    fig.patch.set_facecolor("white")
    ax.bar(x, matched_vals, color=colors, width=0.7, edgecolor="black", linewidth=0.8, hatch="//")
    ax.bar(x, total_vals - matched_vals, bottom=matched_vals, color=[lighten_color(c) for c in colors], width=0.7, linewidth=0)

    for xi, tot, matched in zip(x, total_vals, matched_vals):
        ax.text(
            xi,
            matched if matched > 0 else tot,
            f"{matched / tot * 100:.0f}%",
            ha="center",
            va="bottom",
            fontsize=FLABEL,
            fontweight="bold",
        )

    ax.set_yscale("log")
    ax.set_ylim(1, max(total_vals) * 5)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: format_count(v)))
    ax.tick_params(axis="y", labelsize=FTICK)
    ax.set_ylabel("Count", fontsize=FLABEL)
    ax.set_xticks(x)
    ax.set_xticklabels([b[1] for b in bars], fontsize=FLABEL)
    for label, x0, x1 in groups:
        ax.plot(
            [x0 - 0.35, x1 + 0.35],
            [-0.08, -0.08],
            color="black",
            linewidth=0.8,
            transform=ax.get_xaxis_transform(),
            clip_on=False,
        )
        ax.text((x0 + x1) / 2, -0.1, label, ha="center", va="top", fontsize=FLABEL, transform=ax.get_xaxis_transform())

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#888888", edgecolor="black", linewidth=0.8, hatch="//"),
        plt.Rectangle((0, 0), 1, 1, facecolor=lighten_color("#888888"), linewidth=0),
    ]
    ax.legend(legend_handles, ["Matched in srWGS", "Not matched"], fontsize=FTICK, loc="upper right", frameon=False)
    _title_with_subtitle(ax, "srWGS Overlap by Type and Size", f"N = {format_count(int(total_vals.sum()))} sites")
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


# ==============================
# Plot 2: counts_per_genome_distribution.png
# ==============================
def plot_per_genome_distribution(sample_df, outpath, allele_df=None, random_seed=42, manuscript_mode=False):
    active_buckets = _get_active_buckets(sample_df)
    has_allele = allele_df is not None
    nrows = 2 if has_allele else 1
    af_labels = _get_af_labels(sample_df)
    has_af = len(af_labels) > 0
    ncols = 4 if has_af else 3
    n_samples = len(sample_df)

    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 7 * nrows))
    if nrows == 1:
        axes = axes[np.newaxis, :]
    fig.patch.set_facecolor("white")

    _plot_genome_row(
        axes[0],
        sample_df,
        active_buckets,
        label="Sites",
        n_samples=n_samples,
        has_af=has_af,
        random_seed=random_seed,
        manuscript_mode=manuscript_mode,
    )
    if has_allele:
        allele_buckets = _get_active_buckets(allele_df)
        _plot_genome_row(
            axes[1],
            allele_df,
            allele_buckets,
            label="Alleles",
            n_samples=n_samples,
            has_af=has_af,
            random_seed=random_seed,
            manuscript_mode=manuscript_mode,
        )

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def _get_active_buckets(df):
    tr_region_suffixes = [" - TR", " - Not TR"] + [f" - {r}" for r in REGION_ORDER]
    return [
        c.replace("All - ", "")
        for c in df.columns
        if c.startswith("All - ") and not any(c.endswith(suf) for suf in tr_region_suffixes) and "[" not in c and df[c].sum() > 0
    ]


def _plot_genome_row(
    axes_row, df, active_buckets, label="Sites", n_samples=None, has_af=False, random_seed=42, manuscript_mode=False
):
    ax_jit, ax_size = axes_row[0], axes_row[1]
    ax_af, ax_reg = (axes_row[2], axes_row[3]) if has_af else (None, axes_row[2])
    rng = np.random.default_rng(random_seed)

    bucket_gtype = {b: bucket_to_grouped_type(b, manuscript_mode) for b in active_buckets}
    active_gtypes = [gt for gt in GROUPED_TYPE_ORDER if any(v == gt for v in bucket_gtype.values())]

    n_str = format_count(n_samples) if n_samples else "?"

    # Panel 1: Jitter plot
    plot_labels = ["All"] + active_gtypes
    col_spacing = 1.35
    for i, gt in enumerate(plot_labels):
        pos = i * col_spacing
        if gt == "All":
            vals = np.zeros(len(df))
            for b in active_buckets:
                col_name = f"All - {b}"
                if col_name in df.columns:
                    vals += df[col_name].values.astype(float)
            color = "#555555"
        else:
            gt_bkts = [b for b, g in bucket_gtype.items() if g == gt]
            vals = np.zeros(len(df))
            for b in gt_bkts:
                col_name = f"All - {b}"
                if col_name in df.columns:
                    vals += df[col_name].values.astype(float)
            color = GROUPED_COLORS.get(gt, "#888")
        jitter = rng.uniform(-0.25, 0.25, len(vals))
        ax_jit.scatter(pos + jitter, vals, s=8, color=color, alpha=0.4, linewidths=0, zorder=3)
        median_val = np.median(vals)
        ax_jit.hlines(median_val, pos - 0.35, pos + 0.35, colors="black", linewidth=2, zorder=4)
        ax_jit.text(pos + 0.4, median_val, format_count(median_val), va="center", fontsize=5.5, color=color)

    ax_jit.set_xlim(-0.7, (len(plot_labels) - 1) * col_spacing + 0.7)
    ax_jit.set_xticks([i * col_spacing for i in range(len(plot_labels))])
    ax_jit.set_xticklabels(plot_labels, rotation=45, ha="right", fontsize=FTICK)
    ax_jit.set_xlabel("Type", fontsize=FLABEL)
    ax_jit.set_ylabel(f"{label} Per Genome", fontsize=FLABEL)
    _title_with_subtitle(ax_jit, f"{label} per Genome by Type", f"N = {n_str} samples")
    ax_jit.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: format_count(v)))
    ax_jit.tick_params(axis="y", labelsize=FTICK)
    ax_jit.spines[["top", "right"]].set_visible(False)

    # Panel 2: Heatmap by size
    base_types = [t for t in ("DEL", "INS", "DUP") if any(b.startswith(t + " ") for b in active_buckets)]
    if manuscript_mode:
        base_types = [t for t in base_types if t != "DUP"]
    size_bins = get_size_bins(active_buckets)
    has_snv_size = "SNV" in active_buckets
    if has_snv_size:
        size_bins = ["0-0"] + size_bins
    display_types = base_types + (["SNV"] if has_snv_size else [])

    if display_types and size_bins:
        size_data = np.full((len(size_bins), len(display_types)), np.nan)
        for j, bt in enumerate(display_types):
            for i, sb in enumerate(size_bins):
                if bt == "SNV":
                    col_name = "All - SNV" if sb == "0-0" else None
                elif sb == "0-0":
                    col_name = None
                else:
                    col_name = f"All - {bt} {sb}"
                if col_name and col_name in df.columns:
                    vals = df[col_name].values.astype(float)
                    if manuscript_mode and bt == "INS":
                        dup_col = f"All - DUP {sb}"
                        if dup_col in df.columns:
                            vals = vals + df[dup_col].values.astype(float)
                    size_data[i, j] = np.median(vals)

        valid = size_data[~np.isnan(size_data) & (size_data > 0)]
        vmin = max(float(valid.min()), 0.5) if len(valid) > 0 else 0.5
        vmax = max(float(valid.max()), 1.0) if len(valid) > 0 else 1.0
        norm = mc.LogNorm(vmin=vmin, vmax=vmax)
        im_s = ax_size.imshow(size_data, aspect="auto", cmap="YlOrRd", norm=norm)
        ax_size.set_xticks(range(len(display_types)))
        ax_size.set_xticklabels(display_types, fontsize=FTICK)
        ax_size.set_xlabel("Type", fontsize=FLABEL)
        ax_size.set_yticks(range(len(size_bins)))
        ax_size.set_yticklabels([format_size_label(s) for s in size_bins], fontsize=FTICK)
        ax_size.set_ylabel("Size", fontsize=FLABEL)
        ax_size.set_title(f"{label} per Genome by Size", fontsize=FSUBTITLE, fontweight="bold")
        plt.colorbar(im_s, ax=ax_size, shrink=0.8)
        threshold = np.nanmedian(size_data) * 4 if not np.all(np.isnan(size_data)) else vmax
        for i in range(len(size_bins)):
            for j in range(len(display_types)):
                v = size_data[i, j]
                if not np.isnan(v) and v > 0:
                    text_col = "white" if v > threshold else "black"
                    ax_size.text(j, i, format_count(v), ha="center", va="center", fontsize=FHEAT, color=text_col)
    else:
        ax_size.set_visible(False)

    # Panel 3: Heatmap by AF bin
    if has_af:
        af_labels = _get_af_labels(df)
        af_base_types = ["SNV", "DEL", "INS", "DUP", "TRV", "Other"]
        af_base_active = [t for t in af_base_types if any(f"All - {t} [{afl}]" in df.columns for afl in af_labels)]
        if manuscript_mode:
            af_base_active = [t for t in af_base_active if t != "DUP"]
        if af_labels and af_base_active:
            af_data = np.full((len(af_labels), len(af_base_active)), np.nan)
            for i, afl in enumerate(af_labels):
                for j, bt in enumerate(af_base_active):
                    col_name = f"All - {bt} [{afl}]"
                    if col_name in df.columns:
                        vals = df[col_name].values.astype(float)
                        if manuscript_mode and bt == "INS":
                            dup_col = f"All - DUP [{afl}]"
                            if dup_col in df.columns:
                                vals = vals + df[dup_col].values.astype(float)
                        af_data[i, j] = np.median(vals)

            valid = af_data[~np.isnan(af_data) & (af_data > 0)]
            vmin = max(float(valid.min()), 0.5) if len(valid) > 0 else 0.5
            vmax = max(float(valid.max()), 1.0) if len(valid) > 0 else 1.0
            norm = mc.LogNorm(vmin=vmin, vmax=vmax)
            im_af = ax_af.imshow(af_data, aspect="auto", cmap="YlOrRd", norm=norm)
            ax_af.set_xticks(range(len(af_base_active)))
            ax_af.set_xticklabels(af_base_active, rotation=45, ha="right", fontsize=FTICK)
            ax_af.set_xlabel("Type", fontsize=FLABEL)
            ax_af.set_yticks(range(len(af_labels)))
            ax_af.set_yticklabels([format_af_label(label) for label in af_labels], fontsize=FTICK)
            ax_af.set_ylabel("AF", fontsize=FLABEL)
            ax_af.set_title(f"{label} per Genome by AF", fontsize=FSUBTITLE, fontweight="bold")
            plt.colorbar(im_af, ax=ax_af, shrink=0.8)
            threshold = np.nanmedian(af_data) * 4 if not np.all(np.isnan(af_data)) else vmax
            for i in range(len(af_labels)):
                for j in range(len(af_base_active)):
                    v = af_data[i, j]
                    if not np.isnan(v) and v > 0:
                        text_col = "white" if v > threshold else "black"
                        ax_af.text(j, i, format_count(v), ha="center", va="center", fontsize=FHEAT, color=text_col)
        else:
            ax_af.set_visible(False)

    # Panel 4: Heatmap by region
    active_regions = [r for r in REGION_ORDER if any(f"All - {b} - {r}" in df.columns for b in active_buckets)]

    if active_regions:
        reg_data = np.full((len(active_regions), len(active_gtypes)), np.nan)
        for i, region in enumerate(active_regions):
            for j, gt in enumerate(active_gtypes):
                gt_bkts = [b for b, g in bucket_gtype.items() if g == gt]
                vals = np.zeros(len(df))
                found = False
                for b in gt_bkts:
                    col_name = f"All - {b} - {region}"
                    if col_name in df.columns:
                        vals += df[col_name].values.astype(float)
                        found = True
                if found:
                    reg_data[i, j] = np.median(vals)

        valid = reg_data[~np.isnan(reg_data) & (reg_data > 0)]
        vmin = max(float(valid.min()), 0.5) if len(valid) > 0 else 0.5
        vmax = max(float(valid.max()), 1.0) if len(valid) > 0 else 1.0
        norm = mc.LogNorm(vmin=vmin, vmax=vmax)
        im_r = ax_reg.imshow(reg_data, aspect="auto", cmap="YlOrRd", norm=norm)
        ax_reg.set_xticks(range(len(active_gtypes)))
        ax_reg.set_xticklabels(active_gtypes, rotation=45, ha="right", fontsize=FTICK)
        ax_reg.set_xlabel("Type", fontsize=FLABEL)
        ax_reg.set_yticks(range(len(active_regions)))
        ax_reg.set_yticklabels(active_regions, fontsize=FTICK)
        ax_reg.set_ylabel("Region", fontsize=FLABEL)
        ax_reg.set_title(f"{label} per Genome by Region", fontsize=FSUBTITLE, fontweight="bold")
        plt.colorbar(im_r, ax=ax_reg, shrink=0.8)
        threshold = np.nanmedian(reg_data) * 4 if not np.all(np.isnan(reg_data)) else vmax
        for i in range(len(active_regions)):
            for j in range(len(active_gtypes)):
                v = reg_data[i, j]
                if not np.isnan(v) and v > 0:
                    text_col = "white" if v > threshold else "black"
                    ax_reg.text(j, i, format_count(v), ha="center", va="center", fontsize=FHEAT, color=text_col)
    else:
        ax_reg.set_visible(False)


# ==============================
# Plot 3: trio_inheritance_distribution.png
# ==============================
def plot_trio_inheritance(denovo_df, outpath, random_seed=42, manuscript_mode=False):
    trio_suffixes = [" - TR (Proband)", " - Not TR (Proband)"] + [f" - {r} (Proband)" for r in REGION_ORDER]
    active_buckets = [
        c.replace("All - ", "").replace(" (Proband)", "")
        for c in denovo_df.columns
        if c.startswith("All - ")
        and c.endswith(" (Proband)")
        and not any(c.endswith(suf) for suf in trio_suffixes)
        and denovo_df[c].sum() > 0
    ]

    size_bins = get_size_bins(active_buckets)
    if "SNV" in active_buckets:
        size_bins = ["0-0"] + size_bins

    type_groups = {}
    for b in active_buckets:
        type_groups.setdefault(bucket_to_grouped_type(b, manuscript_mode), []).append(b)
    type_order = [gt for gt in GROUPED_TYPE_ORDER if gt in type_groups]

    n_trios = format_count(len(denovo_df))

    trio_allele_suffixes = [" - TR (Proband Alleles)", " - Not TR (Proband Alleles)"] + [
        f" - {r} (Proband Alleles)" for r in REGION_ORDER
    ]
    has_allele = any(
        c.startswith("All - ")
        and c.endswith(" (Proband Alleles)")
        and not any(c.endswith(suf) for suf in trio_allele_suffixes)
        and denovo_df[c].sum() > 0
        for c in denovo_df.columns
    )

    has_region = any(f"All - {b} - {r} (Proband)" in denovo_df.columns for b in active_buckets for r in REGION_ORDER)
    ncols = 4 if has_region else 3

    nrows = 2 if has_allele else 1
    fig = plt.figure(figsize=(22 if has_region else 17, 7 * nrows))
    fig.patch.set_facecolor("white")

    gs = GridSpec(
        nrows, ncols, figure=fig, wspace=0.45, hspace=0.2, top=0.95, bottom=0.07 if nrows > 1 else 0.14, left=0.10, right=0.97
    )
    ax_inh = fig.add_subplot(gs[0, 0])
    ax_type = fig.add_subplot(gs[0, 1])
    ax_size = fig.add_subplot(gs[0, 2])

    _plot_inheritance_panel(
        ax_inh,
        denovo_df,
        active_buckets,
        count_suffix="",
        title="Site Inheritance",
        subtitle=f"N = {n_trios} trios",
        random_seed=random_seed,
    )
    _plot_dnr_beeswarm(
        ax_type,
        denovo_df,
        active_buckets,
        type_order,
        strat_type="type",
        type_groups=type_groups,
        label="Site",
        random_seed=random_seed,
    )
    _plot_dnr_beeswarm(ax_size, denovo_df, active_buckets, size_bins, strat_type="size", label="Site", random_seed=random_seed)
    if has_region:
        ax_reg = fig.add_subplot(gs[0, 3])
        _plot_dnr_beeswarm(
            ax_reg, denovo_df, active_buckets, REGION_ORDER, strat_type="region", label="Site", random_seed=random_seed
        )

    if has_allele:
        ax_inh2 = fig.add_subplot(gs[1, 0])
        ax_type2 = fig.add_subplot(gs[1, 1])
        ax_size2 = fig.add_subplot(gs[1, 2])
        _plot_inheritance_panel(
            ax_inh2,
            denovo_df,
            active_buckets,
            count_suffix=" Alleles",
            title="Allele Inheritance",
            subtitle=f"N = {n_trios} trios",
            random_seed=random_seed,
        )
        _plot_dnr_beeswarm(
            ax_type2,
            denovo_df,
            active_buckets,
            type_order,
            strat_type="type",
            type_groups=type_groups,
            count_suffix=" Alleles",
            label="Allele",
            random_seed=random_seed,
        )
        _plot_dnr_beeswarm(
            ax_size2,
            denovo_df,
            active_buckets,
            size_bins,
            strat_type="size",
            count_suffix=" Alleles",
            label="Allele",
            random_seed=random_seed,
        )
        if has_region:
            ax_reg2 = fig.add_subplot(gs[1, 3])
            _plot_dnr_beeswarm(
                ax_reg2,
                denovo_df,
                active_buckets,
                REGION_ORDER,
                strat_type="region",
                count_suffix=" Alleles",
                label="Allele",
                random_seed=random_seed,
            )

    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def _compute_dnr_region(denovo_df, buckets, region, count_suffix=""):
    zero = pd.Series(0, index=denovo_df.index)
    proband = sum(denovo_df.get(f"All - {b} - {region} (Proband{count_suffix})", zero).astype(float) for b in buckets)
    mendelian = sum(denovo_df.get(f"All - {b} - {region} (Mendelian{count_suffix})", zero).astype(float) for b in buckets)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(proband > 0, (proband - mendelian) / proband, np.nan)


def _compute_dnr(denovo_df, buckets, count_suffix=""):
    proband = sum(
        denovo_df.get(f"All - {b} (Proband{count_suffix})", pd.Series(0, index=denovo_df.index)).astype(float) for b in buckets
    )
    mendelian = sum(
        denovo_df.get(f"All - {b} (Mendelian{count_suffix})", pd.Series(0, index=denovo_df.index)).astype(float) for b in buckets
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(proband > 0, (proband - mendelian) / proband, np.nan)


def _compute_inh_rate(denovo_df, buckets, count_suffix=""):
    proband = sum(
        denovo_df.get(f"All - {b} (Proband{count_suffix})", pd.Series(0, index=denovo_df.index)).astype(float) for b in buckets
    )
    mendelian = sum(
        denovo_df.get(f"All - {b} (Mendelian{count_suffix})", pd.Series(0, index=denovo_df.index)).astype(float) for b in buckets
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(proband > 0, mendelian / proband, np.nan)


def _compute_rate(denovo_df, buckets, numerator_type, denominator_type):
    num = sum(denovo_df.get(f"All - {b} ({numerator_type})", pd.Series(0, index=denovo_df.index)).astype(float) for b in buckets)
    denom = sum(
        denovo_df.get(f"All - {b} ({denominator_type})", pd.Series(0, index=denovo_df.index)).astype(float) for b in buckets
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denom > 0, num / denom, np.nan)


def _compute_parental_trans_rate(denovo_df, buckets, count_suffix=""):
    zero = pd.Series(0, index=denovo_df.index)
    paternal = sum(denovo_df.get(f"All - {b} (Paternal{count_suffix})", zero).astype(float) for b in buckets)
    maternal = sum(denovo_df.get(f"All - {b} (Maternal{count_suffix})", zero).astype(float) for b in buckets)
    pat_total = sum(denovo_df.get(f"All - {b} (Paternal Total{count_suffix})", zero).astype(float) for b in buckets)
    mat_total = sum(denovo_df.get(f"All - {b} (Maternal Total{count_suffix})", zero).astype(float) for b in buckets)
    denom = pat_total + mat_total
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denom > 0, (paternal + maternal) / denom, np.nan)


MAIN_INHERITANCE_METRICS = {"Proband\nInheritance Rate", "Proband\nDe Novo Rate", "Parental\nTransmission Rate"}


def _plot_inheritance_panel(
    ax, denovo_df, active_buckets, count_suffix="", title="Site Inheritance", subtitle=None, random_seed=42
):
    inh_rates = _compute_inh_rate(denovo_df, active_buckets, count_suffix)
    dnr_rates = _compute_dnr(denovo_df, active_buckets, count_suffix)
    rng = np.random.default_rng(random_seed)

    has_paternal = any(c.endswith(f" (Paternal{count_suffix})") for c in denovo_df.columns)
    has_paternal_total = any(c.endswith(f" (Paternal Total{count_suffix})") for c in denovo_df.columns)

    metrics = [
        ("Proband\nInheritance Rate", inh_rates, "#555555"),
        ("Proband\nDe Novo Rate", dnr_rates, "#F13D15"),
    ]
    if has_paternal:
        metrics += [
            (
                "Inheritance\nRate (Paternal)",
                _compute_rate(denovo_df, active_buckets, f"Paternal{count_suffix}", f"Mendelian{count_suffix}"),
                "#2376B2",
            ),
            (
                "Inheritance\nRate (Maternal)",
                _compute_rate(denovo_df, active_buckets, f"Maternal{count_suffix}", f"Mendelian{count_suffix}"),
                "#D43925",
            ),
            ("Parental\nTransmission Rate", _compute_parental_trans_rate(denovo_df, active_buckets, count_suffix), "#888888"),
        ]
    if has_paternal and has_paternal_total:
        metrics += [
            (
                "Transmission\nRate (Paternal)",
                _compute_rate(denovo_df, active_buckets, f"Paternal{count_suffix}", f"Paternal Total{count_suffix}"),
                "#2376B2",
            ),
            (
                "Transmission\nRate (Maternal)",
                _compute_rate(denovo_df, active_buckets, f"Maternal{count_suffix}", f"Maternal Total{count_suffix}"),
                "#D43925",
            ),
        ]

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-len(metrics) - 0.5, 0.5)
    ax.set_xlabel("Rate", fontsize=FLABEL)
    _title_with_subtitle(ax, title, subtitle)
    ax.axvline(x=0, color="gray", linewidth=0.8)
    ax.axvline(x=1, color="gray", linewidth=0.8)
    for x_grid in np.arange(0.2, 1.0, 0.2):
        ax.axvline(x=x_grid, color="lightgray", linewidth=0.5, linestyle="--")

    for i, (label, vals, col) in enumerate(metrics):
        y = -(i + 1)
        is_main = label in MAIN_INHERITANCE_METRICS
        valid = vals[~np.isnan(vals)]
        ax.axhspan(y - 0.25, y + 0.25, color=col, alpha=0.16 if is_main else 0.08)

        jitter = rng.uniform(-0.15, 0.15, len(valid))
        ax.scatter(valid, y + jitter, s=20, color=col, alpha=0.7, linewidths=0, zorder=3)

        mean_val = np.nanmean(vals)
        ax.vlines(mean_val, y - 0.3, y + 0.3, color=col, linewidth=3, zorder=4)
        ax.text(
            mean_val,
            y + 0.33,
            f"{mean_val * 100:.1f}%",
            ha="center",
            va="bottom",
            fontsize=7.5 if is_main else 6.5,
            fontweight="bold" if is_main else "normal",
            color=col,
        )

        ax.text(
            -0.04,
            y,
            label,
            ha="right",
            va="center",
            transform=ax.get_yaxis_transform(),
            clip_on=False,
            fontsize=8 if is_main else 7,
            fontweight="bold" if is_main else "normal",
            fontstyle="normal" if is_main else "italic",
        )

    ax.set_yticks([])
    ax.set_ylabel("Metric", fontsize=FLABEL)
    ax.yaxis.set_label_coords(-0.30, 0.5)
    ax.tick_params(axis="x", labelsize=FTICK)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v * 100:.0f}%"))


def _plot_dnr_beeswarm(
    ax,
    denovo_df,
    active_buckets,
    strat_labels,
    strat_type="size",
    type_groups=None,
    count_suffix="",
    label="Site",
    random_seed=42,
):
    rng = np.random.default_rng(random_seed)
    if strat_type == "size":
        size_bins = strat_labels
        strat_buckets = {"All": active_buckets}
        strat_buckets.update({s: (["SNV"] if s == "0-0" else [b for b in active_buckets if b.endswith(s)]) for s in size_bins})
        strat_order = ["All"] + size_bins
        title = f"{label} De Novo Rate by Size"
    elif strat_type == "region":
        active_regions = [
            r
            for r in strat_labels
            if any(f"All - {b} - {r} (Proband{count_suffix})" in denovo_df.columns for b in active_buckets)
        ]
        strat_order = ["All"] + active_regions
        strat_buckets = None
        title = f"{label} De Novo Rate by Region"
    else:
        strat_buckets = {"All": active_buckets}
        strat_buckets.update({t: type_groups[t] for t in strat_labels if t in type_groups})
        strat_order = ["All"] + [t for t in strat_labels if t in type_groups]
        title = f"{label} De Novo Rate by Type"

    all_dnr_data = {}
    for s in strat_order:
        if strat_type == "region":
            if s == "All":
                all_dnr_data[s] = _compute_dnr(denovo_df, active_buckets, count_suffix)
            else:
                all_dnr_data[s] = _compute_dnr_region(denovo_df, active_buckets, s, count_suffix)
        else:
            bkts = strat_buckets.get(s, [])
            if bkts:
                all_dnr_data[s] = _compute_dnr(denovo_df, bkts, count_suffix)
    all_valid = np.concatenate([v[~np.isnan(v)] for v in all_dnr_data.values()]) if all_dnr_data else np.array([0.1])
    max_dnr = float(np.nanmax(all_valid)) if len(all_valid) > 0 else 0.1
    x_max = max(max_dnr * 1.45, 0.05)
    grid_step = 0.02 if x_max < 0.12 else (0.05 if x_max < 0.3 else 0.1)

    n_strats = len(strat_order)
    ax.set_xlim(-0.005, x_max)
    ax.set_ylim(-n_strats - 0.5, 0.5)
    ax.set_xlabel("De Novo Rate", fontsize=FLABEL, labelpad=8)
    ax.set_title(title, fontsize=FSUBTITLE, fontweight="bold")
    for x_grid in np.arange(grid_step, x_max, grid_step):
        ax.axvline(x=x_grid, color="lightgray", linewidth=0.5, linestyle="--")
    ax.axvline(x=0, color="gray", linewidth=0.8)

    for i, s in enumerate(strat_order):
        y = -(i + 1)
        s_display = format_size_label(s) if strat_type == "size" and s != "All" else s
        if s not in all_dnr_data:
            ax.text(-0.025, y, s_display, ha="right", va="center", fontsize=7, transform=ax.get_yaxis_transform(), clip_on=False)
            continue

        dnr = all_dnr_data[s]
        valid = dnr[~np.isnan(dnr)]
        col = "#00B0CF" if s != "All" else "#555555"
        ax.axhspan(y - 0.25, y + 0.25, color=col, alpha=0.07)

        jitter = rng.uniform(-0.15, 0.15, len(valid))
        ax.scatter(valid, y + jitter, s=18, color=col, alpha=0.75, linewidths=0, zorder=3)

        mean_val = np.nanmean(dnr)
        ax.vlines(mean_val, y - 0.3, y + 0.3, color=col, linewidth=3, zorder=4)
        ax.text(mean_val + grid_step * 0.15, y + 0.33, f"{mean_val * 100:.1f}%", ha="left", va="bottom", fontsize=6.5, color=col)
        ax.text(-0.025, y, s_display, ha="right", va="center", fontsize=7, transform=ax.get_yaxis_transform(), clip_on=False)

    ax.set_yticks([])
    ylabel = "Size" if strat_type == "size" else ("Region" if strat_type == "region" else "Type")
    ax.set_ylabel(ylabel, fontsize=FLABEL)
    ax.yaxis.set_label_coords(-0.18, 0.5)
    ax.tick_params(axis="x", labelsize=FTICK)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v * 100:.0f}%"))


# ==============================
# Plots 4-9: variant list (AF/size/GT distributions)
# ==============================

VL_TYPE_ORDER = ["SNV", "DEL", "INS+DUP", "TR_SNV", "TR_INS", "TR_DEL"]
VL_TYPE_COLORS = {
    "SNV": "#7459B2",
    "DEL": "#D43925",
    "INS+DUP": "#D474E0",
    "TR_SNV": "#FA931E",
    "TR_INS": "#FA931E",
    "TR_DEL": "#B56A00",
}
VL_ALL_COLOR = "#262626"


def _vl_display_label(typ, manuscript_mode):
    return "INS" if manuscript_mode and typ == "INS+DUP" else typ


HWE_COLORS = {"in_hwe": "#4DAC26", "nominal": "#81F850", "bonferroni": "#AC26A1"}


def _freq_hist(vals_log, n_breaks=150, k=4, xlims=None):
    """Fraction-of-variants histogram + rolling mean, on log scale."""
    clean = vals_log[np.isfinite(vals_log)]
    if len(clean) == 0:
        return None, None, None
    if xlims is None:
        xlims = (float(clean.min()), float(clean.max()))
    breaks = np.linspace(xlims[0], xlims[1], n_breaks + 1)
    mids = (breaks[:-1] + breaks[1:]) / 2
    h, _ = np.histogram(clean, bins=breaks)
    frac = h / max(len(clean), 1)
    # Rolling mean over raw fracs (zeros included, matching R's rollapply behavior)
    roll = uniform_filter1d(frac.astype(float), size=k, mode="nearest")
    return mids, frac, roll


def _log_freq_yaxis(ax, ylims, fontsize=FTICK):
    """Log-scale y-axis with 'X%' labels."""
    ax.set_yscale("log")
    lo_e = max(int(np.floor(np.log10(ylims[0]))), -5)
    hi_e = max(lo_e + 1, int(np.floor(np.log10(ylims[1]))))  # floor not ceil → stay within data
    ticks = [t for t in [10.0**i for i in range(lo_e, hi_e + 1)] if ylims[0] * 0.5 <= t <= ylims[1] * 2]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t * 100:.3g}%" for t in ticks], fontsize=fontsize)
    ax.yaxis.set_minor_formatter(plt.NullFormatter())
    ax.set_ylim(ylims)  # enforce AFTER yticks to prevent axis extension
    ax.set_ylabel("Fraction of Variants", fontsize=FLABEL)


def _log_af_xaxis(ax, xlims_log, fontsize=FTICK):
    """Log-scale x-axis for AF with percentage labels."""
    candidates = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
    ticks_log, labels = [], []
    for v in candidates:
        lv = np.log10(v)
        if xlims_log[0] - 0.05 <= lv <= xlims_log[1] + 0.05:
            ticks_log.append(lv)
            pct = v * 100
            labels.append(f"{pct:.3g}%")
    if ticks_log:
        ax.set_xticks(ticks_log)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=fontsize)
    ax.set_xlabel("AF", fontsize=FLABEL)


def _log_size_xaxis(ax, xlims_log, fontsize=FTICK):
    """Log-scale x-axis for sizes with bp/kb labels."""
    pairs = [
        (1, "1bp"),
        (5, "5bp"),
        (10, "10bp"),
        (50, "50bp"),
        (100, "100bp"),
        (500, "500bp"),
        (1000, "1kb"),
        (5000, "5kb"),
        (10000, "10kb"),
        (50000, "50kb"),
        (100000, "100kb"),
        (500000, "500kb"),
        (1000000, "1Mb"),
    ]
    ticks_log, labels = [], []
    for v, lb in pairs:
        lv = np.log10(v)
        if xlims_log[0] - 0.05 <= lv <= xlims_log[1] + 0.05:
            ticks_log.append(lv)
            labels.append(lb)
    if ticks_log:
        ax.set_xticks(ticks_log)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=fontsize)
    ax.set_xlabel("Size", fontsize=FLABEL)


def _add_log_gridlines(ax, xlims_log, axis="x"):
    """Add gray gridlines at log-scale intervals."""
    lo = int(np.floor(xlims_log[0])) - 1
    hi = int(np.ceil(xlims_log[1])) + 1
    for exp in range(lo, hi + 1):
        for mantissa in range(1, 10):
            lv = np.log10(mantissa * 10**exp)
            if xlims_log[0] <= lv <= xlims_log[1]:
                col = "#d9d9d9" if mantissa == 1 else ("#ebebeb" if mantissa == 5 else "#f5f5f5")
                lw = 0.8 if mantissa == 1 else (0.6 if mantissa == 5 else 0.4)
                ax.axvline(lv, color=col, linewidth=lw, zorder=0)


def classify_hwe_type(allele_class, svlen, manuscript_mode=False):
    """Vectorized over full columns (not row-wise apply) since allele_class/svlen
    can span 100M+ rows."""
    cls = allele_class.str.lower()
    is_sv = svlen.abs() >= 50
    is_dup = cls.str.contains("dup") | (cls == "numt")
    is_del = cls.str.contains("del")
    is_ins = cls.str.contains("ins")
    dup_label = np.where(is_sv, "INS_SV", "INS_SHORT") if manuscript_mode else np.where(is_sv, "DUP_SV", "DUP_SHORT")
    return pd.Categorical(
        np.select(
            [cls == "snv", cls == "tr_snv", cls == "tr_ins", cls == "tr_del", is_dup, is_del, is_ins],
            [
                "SNV",
                "TR_SNV",
                "TR_INS",
                "TR_DEL",
                dup_label,
                np.where(is_sv, "DEL_SV", "DEL_SHORT"),
                np.where(is_sv, "INS_SV", "INS_SHORT"),
            ],
            default="OTH",
        )
    )


def _compute_hwe_pvalues_vec(n_ref, n_het, n_hom):
    N = n_ref + n_het + n_hom
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.where(N > 0, (n_het + 2 * n_hom) / (2 * N), 0.0)
        p = 1 - q
        e_ref = p**2 * N
        e_het = 2 * p * q * N
        e_hom = q**2 * N
        safe = (e_ref > 0) & (e_het > 0) & (e_hom > 0)
        chi2_stat = np.where(
            safe,
            (n_ref - e_ref) ** 2 / np.where(e_ref > 0, e_ref, 1)
            + (n_het - e_het) ** 2 / np.where(e_het > 0, e_het, 1)
            + (n_hom - e_hom) ** 2 / np.where(e_hom > 0, e_hom, 1),
            np.inf,
        )
    pvals = np.where(N > 0, chi2_dist.sf(chi2_stat, 1), np.nan)
    return pvals


def _plot_hwe_triangle(ax, df_sub, title, bonf_threshold, show_labels=True, sample_cap=100000, random_seed=42):
    n_full = len(df_sub)

    # Triangle boundary + HWE parabola
    ax.plot([0, 1, 0.5, 0], [0, 0, 1, 0], "k-", linewidth=0.8)
    q_vals = np.linspace(0, 1, 200)
    ax.plot(q_vals, 2 * q_vals * (1 - q_vals), "-", color=HWE_COLORS["in_hwe"], linewidth=1.2, alpha=0.9)

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("", fontsize=HWE_FLABEL)
    ax.set_ylabel("Fraction of Genotypes", fontsize=HWE_FLABEL)

    y_pct = np.arange(0, 1.01, 0.2)
    ax.set_yticks(y_pct)
    ax.set_yticklabels([f"{int(v * 100)}%" for v in y_pct], fontsize=HWE_FTICK)
    ax.set_xticks([])
    ax.spines[["top", "right", "bottom"]].set_visible(False)

    ax.text(0, -0.04, "Ref.", ha="center", va="top", fontsize=HWE_FTICK, fontweight="bold")
    ax.text(1, -0.04, "Hom.", ha="center", va="top", fontsize=HWE_FTICK, fontweight="bold")
    ax.text(0.5, 1.04, "Het.", ha="center", va="bottom", fontsize=HWE_FTICK, fontweight="bold")

    _title_with_subtitle(ax, title, f"N = {format_count(n_full)} sites", title_fontsize=HWE_FSUBTITLE, subtitle_fontsize=HWE_FSUB)

    if n_full == 0:
        return

    n_total = (df_sub["n_ref"] + df_sub["n_het"] + df_sub["n_hom"]).values
    valid = n_total > 0
    n_ref = df_sub["n_ref"].values[valid].astype(float)
    n_het = df_sub["n_het"].values[valid].astype(float)
    n_hom = df_sub["n_hom"].values[valid].astype(float)
    nt = n_total[valid]
    x_all = n_hom / nt + n_het / (2 * nt)
    y_all = n_het / nt

    # Nominal (alpha=0.05) and Bonferroni-corrected acceptance-region boundaries around the HWE
    # curve: for allele frequency q, the chi-square statistic simplifies to N*d^2/(4*p^2*q^2), where
    # d is the deviation of the heterozygote fraction from its expected value 2pq — so the boundary
    # at a given critical chi-square value is just the HWE curve scaled by (1 +/- sqrt(chi2_crit/N)).
    n_typical = np.median(nt)
    if n_typical > 0:
        q_fine = np.linspace(0.001, 0.999, 200)
        main_curve = 2 * q_fine * (1 - q_fine)
        feasible_max = 2 * np.minimum(q_fine, 1 - q_fine)
        chi2_nom = chi2_dist.ppf(0.95, 1)
        chi2_bonf = chi2_dist.ppf(1 - bonf_threshold, 1)
        for chi2_crit, color in [(chi2_bonf, HWE_COLORS["bonferroni"]), (chi2_nom, HWE_COLORS["nominal"])]:
            d = main_curve * np.sqrt(chi2_crit / n_typical)
            for bound in (main_curve + d, main_curve - d):
                ax.plot(q_fine, np.clip(bound, 0, feasible_max), "-", color=color, linewidth=0.8, alpha=0.7, zorder=1)

    pvals = _compute_hwe_pvalues_vec(n_ref, n_het, n_hom)
    cats = np.where(pvals >= 0.05, "in_hwe", np.where(pvals >= bonf_threshold, "nominal", "bonferroni"))

    n_valid = valid.sum()
    n_pass = (cats == "in_hwe").sum()
    n_nom = (cats == "nominal").sum()
    n_bonf = (cats == "bonferroni").sum()

    if show_labels:
        legend_labels = {
            "in_hwe": f"In HWE ({100 * n_pass / max(n_valid, 1):.1f}%)",
            "nominal": f"Nominal ({100 * n_nom / max(n_valid, 1):.1f}%)",
            "bonferroni": f"Not in HWE ({100 * n_bonf / max(n_valid, 1):.1f}%)",
        }
    else:
        legend_labels = {
            "in_hwe": f"{100 * n_pass / max(n_valid, 1):.1f}%",
            "nominal": f"{100 * n_nom / max(n_valid, 1):.1f}%",
            "bonferroni": f"{100 * n_bonf / max(n_valid, 1):.1f}%",
        }

    # Subsample total to sample_cap, proportional across categories
    rng = np.random.default_rng(random_seed)
    idx_all = np.arange(len(x_all))
    if len(idx_all) > sample_cap:
        idx_all = rng.choice(idx_all, sample_cap, replace=False)

    # Draw in_hwe first (bottom layer), bonferroni last (top layer)
    for cat in ["in_hwe", "nominal", "bonferroni"]:
        mask = cats[idx_all] == cat
        idx = idx_all[mask]
        if len(idx) == 0:
            continue
        ax.scatter(
            x_all[idx], y_all[idx], s=4, color=HWE_COLORS[cat], alpha=0.4, linewidths=0, label=legend_labels[cat], rasterized=True
        )

    ax.legend(fontsize=HWE_FLEGEND, loc="upper right", frameon=False, markerscale=2, labelspacing=0.8)


def _read_variant_list_parquet(path, columns=None):
    """Read the (potentially 100M+ row) variant list parquet, dictionary-encoding the
    low-cardinality string columns so they land as pandas categoricals instead of
    blowing up memory as per-row Python string objects. Only `columns` are read when given."""
    if columns is not None:
        available = pq.ParquetFile(path).schema_arrow.names
        columns = [c for c in columns if c in available]
    table = pq.read_table(path, columns=columns)
    for col in ("allele_class", "region"):
        if col in table.column_names:
            idx = table.column_names.index(col)
            table = table.set_column(idx, col, table.column(col).dictionary_encode())
    df = table.to_pandas(split_blocks=True, self_destruct=True)
    del table
    return df


VL_TYPE_GROUPS = ["SNV", "DEL", "INS+DUP", "TR_SNV", "TR_INS", "TR_DEL"]


def _prep_vl_df(df):
    for col in ["n_ref", "n_het", "n_hom", "n_missing"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype("int32")
    if "af" in df.columns:
        df["af"] = pd.to_numeric(df["af"], errors="coerce")
    if "svlen" in df.columns:
        df["svlen"] = pd.to_numeric(df["svlen"], errors="coerce")
    if "allele_class" in df.columns:
        cls = df["allele_class"].str.lower()
        ins_dup = ("ins", "alu_ins", "sva_ins", "line_ins", "dup", "dup_interspersed", "complex_dup", "inv_dup", "numt")
        del_types = ("del", "alu_del", "sva_del", "line_del")
        # Store the group as int8 category codes rather than one Python string per row
        codes = np.select(
            [cls == "snv", cls.isin(del_types), cls.isin(ins_dup), cls == "tr_snv", cls == "tr_ins", cls == "tr_del"],
            np.arange(len(VL_TYPE_GROUPS), dtype="int8"),
            default=-1,
        ).astype("int8")
        df["type_group"] = pd.Categorical.from_codes(codes, categories=VL_TYPE_GROUPS)
    return df


def _plot_freq_distrib(
    ax,
    groups,
    colors,
    title,
    x_axis_fn,
    n_breaks=150,
    k=4,
    legend=True,
    all_log_vals=None,
    subtitle_offset=0,
    tick_fontsize=FTICK,
    title_fontsize=FSUBTITLE,
):
    """Shared log-log frequency histogram plot (matches R plotSizeDistrib/plotFreqDistrib style)."""
    all_fracs = []
    series = {}
    for label, log_vals, color in groups:
        if len(log_vals) < 2:
            continue
        result = _freq_hist(log_vals, n_breaks=n_breaks, k=k)
        if result[0] is None:
            continue
        mids, frac, roll = result
        series[label] = (mids, frac, roll, color, len(log_vals))
        all_fracs.extend(frac[frac > 0].tolist())
    if not series:
        ax.text(0.5, 0.5, "No Data", ha="center", va="center", transform=ax.transAxes)
        ax.spines[["top", "right"]].set_visible(False)
        return

    # Compute ALL series from combined log_vals if provided
    if all_log_vals is not None and len(all_log_vals) >= 2:
        r = _freq_hist(all_log_vals, n_breaks=n_breaks, k=k)
        if r[0] is not None:
            series["ALL"] = (r[0], r[1], r[2], VL_ALL_COLOR, len(all_log_vals))
            all_fracs.extend(r[1][r[1] > 0].tolist())

    pos_fracs = [v for v in all_fracs if v > 0]
    if not pos_fracs:
        return
    y_min = min(pos_fracs) * 0.5
    y_max = np.quantile(pos_fracs, 0.99)
    xlims = None
    for mids, *_ in series.values():
        if xlims is None:
            xlims = [mids[0], mids[-1]]
        else:
            xlims[0] = min(xlims[0], mids[0])
            xlims[1] = max(xlims[1], mids[-1])

    _add_log_gridlines(ax, xlims)
    # Scatter dots first, then rolling mean lines on top
    for mids, frac, _, color, _ in series.values():
        pos = frac > 0
        ax.scatter(mids[pos], frac[pos], s=1.5, color=color, alpha=0.5, linewidths=0, zorder=2)
    for label, (mids, frac, roll, color, n) in series.items():
        lwd = 2.5 if label == "ALL" else 1.5
        legend_label = f"{label} (N = {format_count(n)})" if legend else label
        # Plot line only where bins have data (frac > 0), using rolling mean values
        pos = frac > 0
        ax.plot(mids[pos], roll[pos], color=color, linewidth=lwd, label=legend_label, zorder=3)

    _log_freq_yaxis(ax, (y_min, y_max), fontsize=tick_fontsize)
    x_axis_fn(ax, xlims, fontsize=tick_fontsize)
    ax.spines[["top", "right"]].set_visible(False)
    if legend:
        ax.legend(fontsize=7, loc="upper right", frameon=True, framealpha=0.85, edgecolor="none")

    main_title, *sub_parts = title.split("\n", 1)
    _title_with_subtitle(
        ax, main_title, sub_parts[0] if sub_parts else None, subtitle_offset=subtitle_offset, title_fontsize=title_fontsize
    )


def plot_af_distribution(df, outpath, manuscript_mode=False):
    SIZE_BINS = [
        (None, None, "0bp"),
        (0, 50, format_size_label("1-49")),
        (50, 100, format_size_label("50-99")),
        (100, 500, format_size_label("100-499")),
        (500, 5000, format_size_label("500-4999")),
        (5000, 50000, format_size_label("5000-49999")),
        (50000, None, format_size_label("50000+")),
    ]
    SIZE_COLORS = plt.cm.viridis(np.linspace(0.05, 0.95, len(SIZE_BINS)))
    tr_types = ("TR_SNV", "TR_INS", "TR_DEL")
    type_order = [t for t in VL_TYPE_ORDER if not (manuscript_mode and t in tr_types)]
    sv_types = [t for t in type_order if t not in ("SNV", "TR_SNV")]

    df_all = df.dropna(subset=["af", "type_group"])
    df_all = df_all[(df_all["af"] > 0) & (df_all["af"] <= 1)]
    n_all = len(df_all)
    df_snv = df_all[df_all["type_group"] == "SNV"]
    df_sv = df_all[df_all["type_group"] != "SNV"].dropna(subset=["svlen"])
    df_sv = df_sv[df_sv["svlen"] > 0]

    def _size_bin_subset(lo, hi):
        if lo is None:
            return df_snv
        sub = df_sv if hi is None else df_sv[df_sv["svlen"] < hi]
        return sub[sub["svlen"] >= lo]

    n_bottom = len(SIZE_BINS)
    half = (n_bottom + 1) // 2
    fig = plt.figure(figsize=(16, 8))
    fig.patch.set_facecolor("white")
    gs = GridSpec(2, n_bottom, figure=fig, height_ratios=[3, 1], hspace=0.5, wspace=0.5)
    ax_type = fig.add_subplot(gs[0, :half])
    ax_size = fig.add_subplot(gs[0, half:])
    bottom_axes = [fig.add_subplot(gs[1, i]) for i in range(n_bottom)]
    for ax in bottom_axes:
        ax.set_box_aspect(1)

    # Top left: by variant type
    type_groups = []
    all_log_vals_type = []
    for typ in type_order:
        sub_af = df_all[df_all["type_group"] == typ]["af"].values
        if len(sub_af) < 5:
            continue
        lv = np.log10(sub_af[sub_af > 0])
        type_groups.append((_vl_display_label(typ, manuscript_mode), lv, VL_TYPE_COLORS[typ]))
        all_log_vals_type.extend(lv.tolist())
    _plot_freq_distrib(
        ax_type,
        type_groups,
        None,
        f"AF Distribution\nN = {format_count(n_all)} sites",
        _log_af_xaxis,
        all_log_vals=np.array(all_log_vals_type),
    )

    # Top right: by size bin
    size_groups = []
    for i, (lo, hi, label) in enumerate(SIZE_BINS):
        sub = _size_bin_subset(lo, hi)["af"].values
        if len(sub) < 5:
            continue
        lv = np.log10(sub[sub > 0])
        size_groups.append((label, lv, SIZE_COLORS[i]))
    _plot_freq_distrib(ax_size, size_groups, None, "AF by Size Bin", _log_af_xaxis, all_log_vals=None)

    # Bottom row: per-size AF distributions colored by type
    for i, (lo, hi, sz_label) in enumerate(SIZE_BINS):
        ax = bottom_axes[i]
        sub_sz = _size_bin_subset(lo, hi)
        types_here = ["SNV"] if lo is None else sv_types
        type_grps, all_log = [], []
        for typ in types_here:
            sub = sub_sz[sub_sz["type_group"] == typ]["af"].values
            if len(sub) < 5:
                continue
            lv = np.log10(sub[sub > 0])
            type_grps.append((_vl_display_label(typ, manuscript_mode), lv, VL_TYPE_COLORS[typ]))
            all_log.extend(lv.tolist())
        _plot_freq_distrib(
            ax,
            type_grps,
            None,
            sz_label,
            _log_af_xaxis,
            all_log_vals=np.array(all_log) if all_log else None,
            legend=False,
            tick_fontsize=FTICK - 2,
            title_fontsize=FSUBTITLE - 2,
        )

    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def _add_mei_lines(ax, xlims_log):
    """Add Alu/SVA/L1 reference lines on top axis like QcAnnotations."""
    mei_marks = [
        (np.log10(300), "Alu"),
        ((np.log10(1000) + np.log10(2000)) / 2, "SVA"),
        (np.log10(6000), "L1"),
    ]
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    tick_locs, tick_labs = [], []
    for lv, label in mei_marks:
        if xlims_log[0] <= lv <= xlims_log[1]:
            ax.axvline(lv, color="gray", linewidth=0.6, linestyle=":", alpha=0.6, zorder=1)
            tick_locs.append(lv)
            tick_labs.append(label)
    ax2.set_xticks(tick_locs)
    ax2.set_xticklabels(tick_labs, fontsize=7, fontstyle="italic")
    ax2.tick_params(length=3, top=True)
    ax2.spines[["top", "right"]].set_visible(False)


def plot_size_distribution(df, af_labels, outpath, manuscript_mode=False):
    AF_COLORS = plt.cm.viridis(np.linspace(0.95, 0.05, len(af_labels)))

    df_sv = df.dropna(subset=["svlen", "type_group"])
    df_sv = df_sv[(df_sv["svlen"] > 0) & (df_sv["type_group"] != "SNV")]
    df_af = df_sv.dropna(subset=["af"])
    tr_types = ("TR_SNV", "TR_INS", "TR_DEL")
    sv_types = [t for t in VL_TYPE_ORDER if t not in ("SNV", "TR_SNV") and not (manuscript_mode and t in tr_types)]

    def _xaxis_fn(ax, xlims_log, fontsize=FTICK):
        _log_size_xaxis(ax, xlims_log, fontsize=fontsize)
        _add_mei_lines(ax, xlims_log)

    n_bottom = len(af_labels)
    half = (n_bottom + 1) // 2
    fig = plt.figure(figsize=(16, 8))
    fig.patch.set_facecolor("white")
    gs = GridSpec(2, n_bottom, figure=fig, height_ratios=[3, 1], hspace=0.5, wspace=0.5)
    ax_type = fig.add_subplot(gs[0, :half])
    ax_af_top = fig.add_subplot(gs[0, half:])
    bottom_axes = [fig.add_subplot(gs[1, i]) for i in range(n_bottom)]
    for ax in bottom_axes:
        ax.set_box_aspect(1)

    # Top left: by variant type
    type_groups, all_log_sv = [], []
    for typ in sv_types:
        sub = df_sv[df_sv["type_group"] == typ]["svlen"].values
        if len(sub) < 5:
            continue
        lv = np.log10(sub)
        type_groups.append((_vl_display_label(typ, manuscript_mode), lv, VL_TYPE_COLORS[typ]))
        all_log_sv.extend(lv.tolist())
    _plot_freq_distrib(
        ax_type,
        type_groups,
        None,
        f"Size Distribution\nN = {format_count(len(all_log_sv))} sites (SNV/TR_SNV dropped)",
        _xaxis_fn,
        all_log_vals=np.array(all_log_sv),
        subtitle_offset=14,
    )

    # Top right: by AF bin
    af_groups = []
    for i, label in enumerate(af_labels):
        sub = df_af[_af_bin_mask(df_af["af"].values, label)]
        if len(sub) < 5:
            continue
        lv = np.log10(sub["svlen"].values)
        af_groups.append((format_af_label(label), lv, AF_COLORS[i]))
    _plot_freq_distrib(ax_af_top, af_groups, None, "Size by AF Bin", _xaxis_fn, all_log_vals=None)

    # Bottom row: per-AF size distributions colored by type
    for i, label in enumerate(af_labels):
        ax = bottom_axes[i]
        sub_af = df_af[_af_bin_mask(df_af["af"].values, label)]
        type_grps, all_log = [], []
        for typ in sv_types:
            sub = sub_af[sub_af["type_group"] == typ]["svlen"].values
            if len(sub) < 5:
                continue
            lv = np.log10(sub)
            type_grps.append((_vl_display_label(typ, manuscript_mode), lv, VL_TYPE_COLORS[typ]))
            all_log.extend(lv.tolist())
        _plot_freq_distrib(
            ax,
            type_grps,
            None,
            format_af_label(label),
            _log_size_xaxis,
            all_log_vals=np.array(all_log) if all_log else None,
            legend=False,
            tick_fontsize=FTICK - 2,
            title_fontsize=FSUBTITLE - 2,
        )

    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_gt_overall_distribution(df, outpath, hwe_sample_cap=100000, random_seed=42):
    df = df.dropna(subset=["n_ref", "n_het", "n_hom"])
    n_vars = len(df)
    bonf = 0.05 / n_vars if n_vars > 0 else 0.05

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor("white")

    _plot_hwe_triangle(
        axes[0], df, "Genotype Distribution", bonf, show_labels=True, sample_cap=hwe_sample_cap, random_seed=random_seed
    )
    axes[0].text(
        0.5,
        -0.08,
        "Autosomal variants only",
        ha="center",
        va="top",
        fontsize=HWE_FSUB,
        color="#555555",
        transform=axes[0].transAxes,
    )

    ax = axes[1]
    df_af = df.dropna(subset=["af"])
    n_total = df_af["n_ref"] + df_af["n_het"] + df_af["n_hom"]
    valid = n_total > 0
    df_v = df_af[valid].copy()
    carrier_freq = (df_v["n_het"] + df_v["n_hom"]) / n_total[valid]
    df_v = df_v.assign(carrier_freq=carrier_freq.values).sort_values("carrier_freq")

    # Grid lines at 20% intervals (matching QcAnnotations)
    for pct in np.arange(0.2, 1.01, 0.2):
        ax.axhline(pct, color="#d9d9d9", linewidth=0.6, zorder=0)
        ax.axvline(pct, color="#d9d9d9", linewidth=0.6, zorder=0)
    # Reference diagonal lines
    cf_line = np.linspace(0, 1, 200)
    ax.plot(cf_line, cf_line, color="#d9d9d9", linewidth=0.8, zorder=0)
    ax.plot(cf_line, cf_line / 2, color="#d9d9d9", linewidth=0.8, zorder=0)

    # Subsample scatter points (max 50k)
    rng = np.random.default_rng(random_seed)
    if len(df_v) > 50000:
        idx = rng.choice(len(df_v), 50000, replace=False)
        ax.scatter(
            df_v["carrier_freq"].values[idx],
            df_v["af"].values[idx],
            s=2,
            alpha=0.3,
            color="black",
            linewidths=0,
            rasterized=True,
            label="Site",
        )
    else:
        ax.scatter(
            df_v["carrier_freq"].values,
            df_v["af"].values,
            s=2,
            alpha=0.3,
            color="black",
            linewidths=0,
            rasterized=True,
            label="Site",
        )

    # Rolling mean on sorted carrier_freq (k=100 window)
    k = min(100, max(1, len(df_v) // 200))
    rm_cf = uniform_filter1d(df_v["carrier_freq"].values.astype(float), size=k)
    rm_af = uniform_filter1d(df_v["af"].values.astype(float), size=k)
    ax.plot(rm_cf, rm_af, "r-", linewidth=1.5, label="Rolling Mean", zorder=3)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    pct_ticks = np.arange(0, 1.01, 0.2)
    ax.set_xticks(pct_ticks)
    ax.set_xticklabels([f"{int(v * 100)}%" for v in pct_ticks], fontsize=HWE_FTICK)
    ax.set_yticks(pct_ticks)
    ax.set_yticklabels([f"{int(v * 100)}%" for v in pct_ticks], fontsize=HWE_FTICK)
    ax.set_xlabel("Carrier Frequency", fontsize=HWE_FLABEL)
    ax.set_ylabel("AF", fontsize=HWE_FLABEL)
    ax.text(1.02, 0.5, "All\nHet.", ha="left", va="center", fontsize=HWE_FSUB, color="gray", transform=ax.get_yaxis_transform())
    ax.text(1.02, 1.0, "All\nHom.", ha="left", va="center", fontsize=HWE_FSUB, color="gray", transform=ax.get_yaxis_transform())
    ax.legend(fontsize=HWE_FLEGEND, loc="upper left", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_gt_type_distribution(df, outpath, hwe_sample_cap=100000, random_seed=42, manuscript_mode=False):
    df = df.dropna(subset=["n_ref", "n_het", "n_hom"])
    df["hwe_type"] = classify_hwe_type(df["allele_class"], df["svlen"], manuscript_mode)
    active_types = [t for t in HWE_TYPE_ORDER if (df["hwe_type"] == t).sum() > 0]

    ncols = 5
    nrows = max(1, (len(active_types) + ncols - 1) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    fig.patch.set_facecolor("white")

    for i, typ in enumerate(active_types):
        ax = axes[i // ncols, i % ncols]
        sub = df[df["hwe_type"] == typ]
        bonf = 0.05 / len(sub) if len(sub) > 0 else 0.05
        _plot_hwe_triangle(ax, sub, typ, bonf, show_labels=False, sample_cap=hwe_sample_cap, random_seed=random_seed)
    for j in range(len(active_types), nrows * ncols):
        axes[j // ncols, j % ncols].set_visible(False)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_gt_size_distribution(df, outpath, hwe_sample_cap=100000, random_seed=42):
    SIZE_BINS = [
        ("0-0", 0, 0),
        ("1-49", 1, 50),
        ("50-99", 50, 100),
        ("100-499", 100, 500),
        ("500-4999", 500, 5000),
        ("5000-49999", 5000, 50000),
        ("50000+", 50000, None),
    ]
    df = df.dropna(subset=["n_ref", "n_het", "n_hom", "svlen"])
    df_sv = df[df["svlen"] >= 0]

    ncols, nrows = 3, 2
    col_width = 6
    col0_width = 5.5 * nrows
    fig = plt.figure(figsize=(col0_width + col_width * ncols, 5.5 * nrows))
    fig.patch.set_facecolor("white")
    gs = GridSpec(nrows, ncols + 1, figure=fig, width_ratios=[col0_width] + [col_width] * ncols)

    axes = [fig.add_subplot(gs[:, 0])]
    axes += [fig.add_subplot(gs[i, j + 1]) for i in range(nrows) for j in range(ncols)]

    for ax, (raw_label, lo, hi) in zip(axes, SIZE_BINS):
        if lo == hi:
            sub = df_sv[df_sv["svlen"] == lo]
        elif hi is None:
            sub = df_sv[df_sv["svlen"] >= lo]
        else:
            sub = df_sv[(df_sv["svlen"] >= lo) & (df_sv["svlen"] < hi)]
        bonf = 0.05 / len(sub) if len(sub) > 0 else 0.05
        _plot_hwe_triangle(
            ax, sub, format_size_label(raw_label), bonf, show_labels=False, sample_cap=hwe_sample_cap, random_seed=random_seed
        )

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def plot_gt_region_distribution(df, outpath, hwe_sample_cap=100000, random_seed=42):
    df = df.dropna(subset=["n_ref", "n_het", "n_hom"])
    active_regions = [r for r in REGION_ORDER if (df["region"] == r).sum() > 0]
    if not active_regions:
        return

    ncols = 2
    nrows = max(1, (len(active_regions) + ncols - 1) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows))
    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = axes[np.newaxis, :]
    fig.patch.set_facecolor("white")

    for i, region in enumerate(active_regions):
        ax = axes[i // ncols, i % ncols]
        sub = df[df["region"] == region]
        bonf = 0.05 / len(sub) if len(sub) > 0 else 0.05
        _plot_hwe_triangle(ax, sub, region, bonf, show_labels=False, sample_cap=hwe_sample_cap, random_seed=random_seed)
    for j in range(len(active_regions), nrows * ncols):
        axes[j // ncols, j % ncols].set_visible(False)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


# ==============================
# Main
# ==============================
# Variant list columns each plot reads, so a single-plot run loads only what it needs
VL_COLUMNS = {
    "distribution_af": ["allele_class", "svlen", "af"],
    "distribution_size": ["allele_class", "svlen", "af"],
    "hwe_type": ["allele_class", "svlen", "n_ref", "n_het", "n_hom"],
    "hwe_overall": ["n_ref", "n_het", "n_hom", "af"],
    "hwe_size": ["svlen", "n_ref", "n_het", "n_hom"],
    "hwe_region": ["region", "n_ref", "n_het", "n_hom"],
}

# Plot name (written as <name>.png) → input argument it needs
PLOT_INPUTS = {
    "site_counts": "sites_parquet",
    "srwgs_overlap": "sites_parquet",
    "sites_per_genome": "samples_parquet",
    "inheritance": "denovo_parquet",
    "distribution_af": "variant_list_parquet",
    "distribution_size": "variant_list_parquet",
    "hwe_type": "variant_list_parquet",
    "hwe_overall": "variant_list_parquet",
    "hwe_size": "variant_list_parquet",
    "hwe_region": "variant_list_parquet",
}


def main():
    parser = argparse.ArgumentParser(description="Plot annotation count distributions from CountAnnotations outputs.")
    parser.add_argument("--sites_parquet", help="plotting_sites_parquet from CountAnnotations")
    parser.add_argument("--samples_parquet", help="plotting_samples_parquet from CountAnnotations")
    parser.add_argument("--alleles_parquet", help="plotting_alleles_parquet from CountAnnotations")
    parser.add_argument("--denovo_parquet", help="plotting_denovo_parquet from CountAnnotations")
    parser.add_argument("--variant_list_parquet", help="plotting_variant_list_parquet from CountAnnotations")
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
        dest="random_seed",
        help="Seed for all random subsampling in this script (default: 42)",
    )
    parser.add_argument(
        "--hwe-subsample-count",
        type=int,
        default=500000,
        dest="hwe_subsample_count",
        help="Max variants to subsample when plotting HWE triangle plots",
    )
    parser.add_argument("--output_dir", default=".", help="Output directory for plots")
    parser.add_argument(
        "--plots",
        nargs="+",
        choices=list(PLOT_INPUTS),
        help="Generate only these plots, overwriting existing files (default: every plot whose "
        "input is given and whose <name>.png is not already in --output_dir)",
    )
    args = parser.parse_args()

    out_dir = args.output_dir
    if args.plots:
        missing = [p for p in args.plots if not getattr(args, PLOT_INPUTS[p])]
        if missing:
            parser.error(
                f"missing inputs for --plots {' '.join(missing)}: {', '.join(sorted({'--' + PLOT_INPUTS[p] for p in missing}))}"
            )
        plots = set(args.plots)
    else:
        plots = {
            p for p, inp in PLOT_INPUTS.items() if getattr(args, inp) and not os.path.exists(os.path.join(out_dir, f"{p}.png"))
        }
    print(f"Generating: {', '.join(p for p in PLOT_INPUTS if p in plots) or 'nothing'}")

    def outpath(name):
        return os.path.join(out_dir, f"{name}.png")

    sites_df, sample_df, allele_df, denovo_df, vl_df, af_labels = None, None, None, None, None, []

    if args.sites_parquet:
        sites_df = pd.read_parquet(args.sites_parquet)

    if args.samples_parquet:
        sample_df = pd.read_parquet(args.samples_parquet).set_index("sample_id")
        if args.alleles_parquet:
            allele_df = pd.read_parquet(args.alleles_parquet).set_index("sample_id")

    if args.denovo_parquet:
        denovo_df = pd.read_parquet(args.denovo_parquet).set_index("sample_id")

    # The variant list can span 100M+ rows, so only load it when one of its plots is pending
    vl_plots = [p for p in plots if PLOT_INPUTS[p] == "variant_list_parquet"]
    if vl_plots:
        vl_columns = sorted({c for p in vl_plots for c in VL_COLUMNS[p]})
        vl_df = _prep_vl_df(_read_variant_list_parquet(args.variant_list_parquet, columns=vl_columns))
        af_labels = _get_af_labels(sites_df if sites_df is not None else sample_df)

    os.makedirs(out_dir, exist_ok=True)

    if "site_counts" in plots:
        plot_overall_distribution(sites_df, outpath("site_counts"), manuscript_mode=True)
    if "srwgs_overlap" in plots:
        plot_srwgs_overlap(sites_df, outpath("srwgs_overlap"), manuscript_mode=True)

    if "sites_per_genome" in plots:
        plot_per_genome_distribution(
            sample_df, outpath("sites_per_genome"), allele_df=allele_df, random_seed=args.random_seed, manuscript_mode=True
        )

    if "inheritance" in plots:
        plot_trio_inheritance(denovo_df, outpath("inheritance"), random_seed=args.random_seed, manuscript_mode=True)

    if "distribution_af" in plots:
        plot_af_distribution(vl_df, outpath("distribution_af"), manuscript_mode=True)
    if "distribution_size" in plots:
        plot_size_distribution(vl_df, af_labels, outpath("distribution_size"), manuscript_mode=True)
    if "hwe_type" in plots:
        plot_gt_type_distribution(
            vl_df,
            outpath("hwe_type"),
            hwe_sample_cap=args.hwe_subsample_count,
            random_seed=args.random_seed,
            manuscript_mode=True,
        )
    if "hwe_overall" in plots:
        plot_gt_overall_distribution(
            vl_df, outpath("hwe_overall"), hwe_sample_cap=args.hwe_subsample_count, random_seed=args.random_seed
        )
    if "hwe_size" in plots:
        plot_gt_size_distribution(
            vl_df, outpath("hwe_size"), hwe_sample_cap=args.hwe_subsample_count, random_seed=args.random_seed
        )
    if "hwe_region" in plots and "region" in vl_df.columns and vl_df["region"].notna().any() and (vl_df["region"] != "").any():
        plot_gt_region_distribution(
            vl_df, outpath("hwe_region"), hwe_sample_cap=args.hwe_subsample_count, random_seed=args.random_seed
        )


if __name__ == "__main__":
    main()
