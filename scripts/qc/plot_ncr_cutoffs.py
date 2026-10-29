#!/usr/bin/env python3
"""Plot precision and recall across NCR cutoffs from an EvaluateNCRCutoffs TSV and suggest a cutoff per stratum.

Both metrics count non-ref calls. Precision is prec_numerator / prec_denominator (inherited over inherited plus
de novo trio proband calls, skipped calls excluded) and recall is rec_numerator / rec_denominator (calls of
every sample passing the cutoff over all their calls at sites with NCR). The ROC curve treats inherited proband
calls as positives and de novo ones as negatives. Rows whose cutoff is "missing" count calls at sites without
NCR; they sit outside every curve.

Allele types, size bins, FILTER sets and regions are read from the TSV; allele type and size bin `all` marginals
are skipped. Each allele type gets three figures with one row per size bin holding proband calls, each row showing
precision against recall, the ROC curve, and precision and recall against the cutoff. Only by_filter includes
non-PASS sites:
  <TYPE>.combined.png  - PASS sites across every region, with the suggested cutoffs marked.
  <TYPE>.by_filter.png - one line per FILTER set (PASS first, the rarest sets summed into Other).
  <TYPE>.by_region.png - PASS sites, one line per region.
  <TYPE>.distribution_sites.png, <TYPE>.distribution_calls.png - NCR histograms of sites and of non-ref calls, one
      per FILTER set, in one panel per size bin and region.
Precision axes whose values all sit high are broken to zoom on them while keeping 0 in view. Three suggestions are
made from PASS sites for every allele type, size bin and region and written with site counts to cutoffs.tsv:
  knee   - the precision-recall point furthest from the chord joining the curve's two ends.
  max_f1 - the smallest cutoff reaching the maximum F1.
  youden - the smallest cutoff reaching the maximum of TPR - FPR on the ROC curve.
With --recommended_tsv, one method's cutoff is also written per PASS allele type, size bin and region in the format
FilterNCRCutoffs reads. A stratum too sparse for the method takes the method's cutoff across every region, then
--fallback_cutoff when given, and is otherwise left out; the source column records which applied.
"""

import argparse
import os
from collections import namedtuple

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from brokenaxes import BrokenAxes, brokenaxes  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import NullLocator  # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES_COLOR = "#2a78d6"
# Categorical slots in fixed order, so a FILTER set or region keeps its colour across figures
CATEGORY_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]
SUGGESTION_STYLES = {"knee": ("Knee", "o", ":"), "max_f1": ("Max F1", "s", "--"), "youden": ("Youden J", "^", "-.")}

KEYS = ["allele_type", "size_bin", "filter", "region"]
COUNT_COLUMNS = [
    "prec_numerator",
    "prec_denominator",
    "prec_skipped",
    "rec_numerator",
    "rec_denominator",
    "sites_passing",
    "sites_total",
]
REGION_ORDER = ["US", "RM", "SD", "SR"]
MAX_NAMED_FILTERS = 6
OTHER = "Other"
# A breakdown group with fewer classified proband calls shows only its recall line
MIN_PRECISION_CALLS = 30
# Fewer de novo proband calls leave no knee, Youden J or breakdown ROC line
MIN_DE_NOVO_CALLS = 10
# Title, x column, y column, x label, y label
VIEWS = [
    ("Precision vs Recall", "recall", "precision", "Recall", "Precision"),
    ("ROC", "fpr", "tpr", "FPR", "TPR"),
    ("Precision vs NCR Cutoff", "cutoff", "precision", "NCR Cutoff", "Precision"),
    ("Recall vs NCR Cutoff", "cutoff", "recall", "NCR Cutoff", "Recall"),
]
# Combined figures mark suggestions on cutoff views with vertical lines; breakdowns mark points on the line itself
Series = namedtuple("Series", "label color linestyle rows suggestions cutoff_lines", defaults=(True,))

PANEL_SIZE = (2.4, 2.3)
INFO_WIDTH = 0.75
# Top and bottom figure margins in inches: title and column labels, then x labels and one legend row
MARGINS = (0.85, 1.0)
LEGEND_ROW_HEIGHT = 0.25
DEFINITION_LINE_HEIGHT = 0.15
DEFINITIONS = (
    "•  A cutoff keeps sites with NCR (no-call rate) ≤ the cutoff and removes the rest. Every metric counts non-ref calls "
    "at kept sites.\n"
    "•  Precision = inherited trio proband calls / (inherited + de novo trio proband calls) at kept sites; inherited means "
    "a parent also has a non-ref call at the site.\n"
    "    1 − precision is the de novo rate. Calls with no carrier parent but a missing parent genotype are skipped.\n"
    "•  Recall = share of all samples' calls kept.\n"
    "•  ROC = TPR against FPR, both over trio proband calls:\n"
    "      –  TPR = recall measured on inherited trio proband calls only.\n"
    "      –  FPR = share of de novo trio proband calls kept. These calls are rare, so FPR can be high while precision "
    "stays near 1."
)
# Cumulative column, total column, unit and the noun the caption counts
DISTRIBUTION_MEASURES = {
    "sites": ("sites_passing", "sites_total", "sites", "sites"),
    "calls": ("rec_numerator", "rec_denominator", "calls", "non-ref calls across all samples"),
}
# A FILTER set with fewer sites in a panel gets no histogram there, as a few sites only add noise
MIN_HISTOGRAM_SITES = 30
# Points below the panel's bottom edge that clear the tick labels and x label
NOTE_BELOW_XLABEL = 32
# Filled per measure with the noun it counts
DISTRIBUTION_DEFINITIONS = (
    "•  NCR (no-call rate) = share of a site's sample alleles without a call. Sites without NCR are left out.\n"
    "•  One histogram per FILTER set, as in the by_filter figures (PASS first, the rarest sets summed into Other), each "
    "the count of that set's {noun} per bin.\n"
    f"    A set with under {MIN_HISTOGRAM_SITES} sites in a panel is left out of it.\n"
    "•  The first bin, [0, b], includes NCR = 0; the rest are (a, b]. Bin width is set per panel from every FILTER set "
    "together by NumPy's auto rule\n"
    "    (the narrower of Freedman–Diaconis and Sturges), rounded up to a whole cutoff step.\n"
    "•  The y-axis is log scale, so small sets show beside PASS; empty bins drop to its floor, below 1."
)
# Left edges, in points right of the size bin label, of the by_region cutoff table's region, knee, max F1 and Youden J columns
TABLE_COLUMNS = (0, 24, 54, 92)
# Only the combined and by_region figures mark suggested cutoffs
SUGGESTION_DEFINITIONS = (
    "•  Suggested cutoffs:\n"
    "      –  Knee = where the precision-recall curve bends most (furthest from the straight line joining its ends).\n"
    "      –  Max F1 = smallest cutoff giving the best balance of precision and recall (F1, their harmonic mean).\n"
    "      –  Youden J = smallest cutoff giving the largest TPR − FPR.\n"
    f"    All three show n/a with under {MIN_PRECISION_CALLS} inherited + de novo proband calls, too few to trust precision.\n"
    f"    Knee and Youden J also show n/a with under {MIN_DE_NOVO_CALLS} de novo proband calls, too few to place them."
)
TICKS = [0, 0.5, 1]
PAD = 0.04
# Precision axes whose zoomed range would start at or above BREAK_ABOVE are broken, keeping a short 0 stub
BREAK_ABOVE = 0.5
ZOOM_STEP = 0.01
STUB_HEIGHT = 0.05
ZOOM_PAD = 0.08
BROKEN_XLABEL_PAD = 18
BROKEN_YLABEL_PAD = 38
# Upper length bound written for the open-ended size bin, as FilterNCRCutoffs takes inclusive integer bounds
MAX_LENGTH = 1000000000


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tsv", required=True, help="ncr_cutoffs_tsv output of EvaluateNCRCutoffs")
    parser.add_argument("--output_dir", required=True, help="Directory for plots and the cutoff table")
    parser.add_argument("--recommended_tsv", help="Path for the recommended cutoffs, in the format FilterNCRCutoffs reads")
    parser.add_argument(
        "--method",
        choices=list(SUGGESTION_STYLES),
        default="knee",
        help="Suggestion written to --recommended_tsv (default: knee)",
    )
    parser.add_argument(
        "--fallback_cutoff",
        type=float,
        help="Cutoff for strata too sparse for --method in their region and across regions; unset leaves them out",
    )
    return parser.parse_args()


def apply_style():
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": AXIS,
            "axes.linewidth": 0.8,
            "axes.labelcolor": INK_SECONDARY,
            "axes.titlecolor": INK,
            "axes.titlesize": 9.5,
            "axes.titleweight": "semibold",
            "axes.titlepad": 6,
            "axes.labelsize": 9,
            "axes.grid": True,
            "axes.axisbelow": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.frameon": False,
            "legend.fontsize": 9,
            "figure.titlesize": 11.5,
            "figure.titleweight": "semibold",
            "figure.dpi": 200,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.3,
        }
    )


def load(path):
    """The cutoff sweep, and separately the rows counting calls and sites without NCR."""
    table = pd.read_csv(path, sep="\t", dtype={column: str for column in KEYS + ["cutoff"]}, keep_default_na=False)
    missing_ncr = table[table["cutoff"] == "missing"]
    return table[table["cutoff"] != "missing"].astype({"cutoff": float}), missing_ncr


def derive(table):
    """Adds precision, recall, F1 and the ROC rates; a zero denominator leaves precision and recall undefined."""
    table = table.copy()
    table["precision"] = table["prec_numerator"] / table["prec_denominator"].replace(0, np.nan)
    table["recall"] = table["rec_numerator"] / table["rec_denominator"].replace(0, np.nan)
    precision, recall = table["precision"], table["recall"]
    table["f1"] = (2 * precision * recall / (precision + recall)).fillna(0)
    # Every site with NCR passes the loosest cutoff, so each stratum's maxima are its total inherited and called
    strata_groups = table.groupby(KEYS, sort=False)
    positives = strata_groups["prec_numerator"].transform("max")
    negatives = strata_groups["prec_denominator"].transform("max") - positives
    table["tpr"] = (table["prec_numerator"] / positives.replace(0, np.nan)).fillna(0)
    table["fpr"] = ((table["prec_denominator"] - table["prec_numerator"]) / negatives.replace(0, np.nan)).fillna(0)
    table["youden"] = table["tpr"] - table["fpr"]
    return table


def strata(table):
    """Size bins with proband calls for each allele type, in the TSV's own order, without the `all` marginals."""
    table = table.dropna(subset=["precision", "recall"])
    table = table[
        (table["allele_type"] != "all") & (table["size_bin"] != "all") & (table["filter"] == "all") & (table["region"] == "all")
    ]
    return {allele_type: list(pd.unique(rows["size_bin"])) for allele_type, rows in table.groupby("allele_type", sort=False)}


def rows_of(table, allele_type, size_bin, filter_value="all", region="all"):
    """Every cutoff row of one stratum, in cutoff order."""
    mask = (
        (table["allele_type"] == allele_type)
        & (table["size_bin"] == size_bin)
        & (table["filter"] == filter_value)
        & (table["region"] == region)
    )
    return table[mask].reset_index(drop=True)


def select(table, allele_type, size_bin, filter_value="all", region="all"):
    """Cutoff rows of one stratum that carry proband calls."""
    rows = rows_of(table, allele_type, size_bin, filter_value, region)
    return rows.dropna(subset=["precision", "recall"]).reset_index(drop=True)


def suggest(rows):
    """Row indices of the knee of the precision-recall curve, the max-F1 cutoff and the Youden J cutoff."""
    offsets = rows[["recall", "precision"]].to_numpy() - rows[["recall", "precision"]].iloc[0].to_numpy()
    chord = offsets[-1]
    # Perpendicular distance from the chord is proportional to the 2D cross product with it
    distances = np.abs(chord[0] * offsets[:, 1] - chord[1] * offsets[:, 0])
    # F1 and J plateau once their terms saturate, so take the smallest cutoff reaching the rounded maximum
    suggestions = {
        "knee": int(np.argmax(distances)),
        "max_f1": int(rows["f1"].round(3).idxmax()),
        "youden": int(rows["youden"].round(3).idxmax()),
    }
    # Knee and Youden J follow the de novo calls, so a handful of them makes both noise
    classified, inherited = classified_calls(rows)
    if classified < MIN_PRECISION_CALLS:
        return {}
    if classified - inherited < MIN_DE_NOVO_CALLS:
        return {"max_f1": suggestions["max_f1"]}
    return suggestions


def grouped(table, column, groups):
    """Counts summed within each group of FILTER sets or regions, so a group reads as one stratum."""
    labels = {value: label for label, values in groups.items() for value in values}
    subset = table[table[column].isin(labels)].copy()
    subset[column] = subset[column].map(labels)
    return derive(subset.groupby(KEYS + ["cutoff"], sort=False, as_index=False)[COUNT_COLUMNS].sum())


def filter_groups(table):
    """PASS first, then FILTER sets by call count; past MAX_NAMED_FILTERS the rest are summed into Other."""
    sets = table[
        (table["allele_type"] == "all") & (table["size_bin"] == "all") & (table["region"] == "all") & (table["filter"] != "all")
    ]
    calls = sets.groupby("filter")["rec_denominator"].max()
    ordered = sorted(calls.index, key=lambda value: (value != "PASS", -calls[value]))
    groups = {value: [value] for value in ordered[:MAX_NAMED_FILTERS]}
    if ordered[MAX_NAMED_FILTERS:]:
        groups[OTHER] = ordered[MAX_NAMED_FILTERS:]
    return groups


def region_groups(table):
    present = set(table["region"]) - {"all"}
    return {region: [region] for region in [r for r in REGION_ORDER if r in present] + sorted(present - set(REGION_ORDER))}


def category_colors(labels):
    return {label: MUTED if label == OTHER else color for label, color in zip(labels, CATEGORY_COLORS)}


def plot_grid(title, size_bins, series_for, info_for, legend_handles, definitions, path):
    """One row per size bin across VIEWS, with an info column naming the row and the legend and definitions below."""
    legend_columns = min(len(legend_handles), 3 if len(legend_handles) > 6 else len(legend_handles))
    legend_rows = -(-len(legend_handles) // legend_columns)
    definitions_height = DEFINITION_LINE_HEIGHT * (definitions.count("\n") + 2)
    bottom = MARGINS[1] + LEGEND_ROW_HEIGHT * (legend_rows - 1) + definitions_height
    height = PANEL_SIZE[1] * len(size_bins) + MARGINS[0] + bottom
    width_ratios = [1] * len(VIEWS) + [INFO_WIDTH]
    # brokenaxes places its own sub-axes, so the grid is laid out by hand rather than by constrained layout
    fig = plt.figure(figsize=(PANEL_SIZE[0] * sum(width_ratios), height), layout="none")
    grid = fig.add_gridspec(
        len(size_bins),
        len(width_ratios),
        width_ratios=width_ratios,
        top=1 - MARGINS[0] / height,
        bottom=bottom / height,
        left=0.06,
        right=0.98,
        hspace=0.55,
        wspace=0.5,
    )
    for row_index, size_bin in enumerate(size_bins):
        series = series_for(size_bin)
        lowest = min((value for item in series for value in [item.rows["precision"].min()] if np.isfinite(value)), default=np.nan)
        for col_index, (view_title, x, y, xlabel, ylabel) in enumerate(VIEWS):
            spec = grid[row_index, col_index]
            ax = precision_axes(fig, spec, lowest) if y == "precision" else fig.add_subplot(spec)
            ax.set_xlim(0, 1 + PAD)
            if y != "precision":
                ax.set_ylim(0, 1 + PAD)
            broken = isinstance(ax, BrokenAxes)
            # brokenaxes labels its whole frame, which carries no tick labels, so it needs a wider pad to line up
            ax.set_xlabel(xlabel, labelpad=BROKEN_XLABEL_PAD if broken else 4)
            ax.set_ylabel(ylabel, labelpad=BROKEN_YLABEL_PAD if broken else 4)
            if row_index == 0:
                ax.set_title(view_title, fontsize=10.5, fontweight="semibold", color=INK, pad=12)
            if x == "fpr":
                ax.plot([0, 1], [0, 1], color=AXIS, linestyle="--", linewidth=0.7, zorder=1)
            for item in series:
                draw_series(ax, item, x, y)
            set_ticks(ax, broken)
        info_for(fig.add_subplot(grid[row_index, -1]), size_bin)
    fig.suptitle(title, fontsize=14, y=1 - 0.2 / height)
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, definitions_height / height), ncol=legend_columns)
    fig.text(0.03, 0.1 / height, definitions, ha="left", va="bottom", fontsize=7.5, linespacing=1.5, color=INK_SECONDARY)
    fig.savefig(path)
    plt.close(fig)


def pass_only(table):
    """PASS rows standing in for the all-FILTER rows, so every other function reads them unchanged."""
    return table[table["filter"] == "PASS"].assign(filter="all")


def mask_sparse(rows):
    """Blanks lines too few proband calls back: precision and ROC, or just ROC, as FPR moves only with de novo calls."""
    classified, inherited = classified_calls(rows)
    if classified < MIN_PRECISION_CALLS:
        rows[["precision", "tpr", "fpr"]] = np.nan
    elif classified - inherited < MIN_DE_NOVO_CALLS:
        rows[["tpr", "fpr"]] = np.nan
    return rows


def plot_allele_type(table, allele_type, size_bins, path):
    """Every region together, with the suggested cutoffs and their counts."""

    def series_for(size_bin):
        rows = select(table, allele_type, size_bin)
        suggestions = suggest(rows)
        return [Series(allele_type, SERIES_COLOR, "-", mask_sparse(rows), suggestions)]

    def info_for(ax, size_bin):
        draw_info(ax, select(table, allele_type, size_bin), size_bin)

    plot_grid(
        f"{allele_type.upper()} (PASS only)",
        size_bins,
        series_for,
        info_for,
        suggestion_handles(),
        f"{DEFINITIONS}\n{SUGGESTION_DEFINITIONS}",
        path,
    )


def plot_breakdown(table, allele_type, size_bins, column, groups, path):
    """One line per FILTER set or region beside the dashed stratum across all of them; regions also mark their cutoffs."""
    by_group = grouped(table, column, groups)
    colors = category_colors(groups)

    def group_rows(size_bin, label):
        # Keep groups without proband calls, which still have a recall line
        rows = rows_of(by_group, allele_type, size_bin, **{"filter_value" if column == "filter" else column: label})
        rows = rows.dropna(subset=["recall"]).reset_index(drop=True)
        return rows if rows.empty else mask_sparse(rows)

    def group_suggestions(size_bin, label):
        """Suggested cutoffs of one region, from the same rows as its cutoffs.tsv row."""
        rows = select(by_group, allele_type, size_bin, region=label)
        return {name: rows.loc[index, "cutoff"] for name, index in suggest(rows).items()} if len(rows) else {}

    def series_for(size_bin):
        series = [Series("All", INK, "--", select(table, allele_type, size_bin), None)]
        for label in groups:
            rows = group_rows(size_bin, label)
            if rows.empty:
                continue
            # Only regions are marked; FILTER sets are compared, not tuned
            suggestions = None
            if column == "region":
                suggestions = {
                    name: int(np.flatnonzero(rows["cutoff"] == cutoff)[0])
                    for name, cutoff in group_suggestions(size_bin, label).items()
                }
            series.append(Series(label, colors[label], "-", rows, suggestions, cutoff_lines=False))
        return series

    def info_for(ax, size_bin):
        draw_label(ax, size_bin)
        if column != "region":
            return
        # A table of region rows by method columns; each column is one text block, so the shared line spacing aligns rows
        cutoffs = {label: group_suggestions(size_bin, label) for label in groups}
        columns = [("", list(groups))] + [
            (method, [f"{cutoffs[label][name]:.2f}" if name in cutoffs[label] else "n/a" for label in groups])
            for name, (method, _, _) in SUGGESTION_STYLES.items()
        ]
        for x, (header, values) in zip(TABLE_COLUMNS, columns):
            # Point offsets keep the columns apart however narrow the info axis is
            transform = offset_copy(ax.transAxes, fig=ax.figure, x=x, units="points")
            ax.text(
                0.05,
                0.56,
                header,
                transform=transform,
                ha="left",
                va="top",
                fontsize=8.5,
                linespacing=1.6,
                fontweight="semibold",
                color=INK,
            )
            ax.text(
                0.05,
                0.56,
                "\n" + "\n".join(values),
                transform=transform,
                ha="left",
                va="top",
                fontsize=8.5,
                linespacing=1.6,
                color=INK_SECONDARY,
            )

    def sites(label):
        """Sites with NCR across every size bin of the allele type."""
        rows = rows_of(by_group, allele_type, "all", **{"filter_value" if column == "filter" else column: label})
        return abbreviate(rows["sites_total"].iloc[-1])

    drawn = [label for label in groups if any(not group_rows(size_bin, label).empty for size_bin in size_bins)]
    handles = [Line2D([], [], color=colors[label], linewidth=1.6, label=f"{label} ({sites(label)} sites)") for label in drawn]
    total = abbreviate(rows_of(table, allele_type, "all")["sites_total"].iloc[-1])
    handles.append(Line2D([], [], color=INK, linestyle="--", linewidth=1.6, label=f"All ({total} sites)"))
    title = f"{allele_type.upper()} by FILTER" if column == "filter" else f"{allele_type.upper()} by Region (PASS only)"
    if column == "region":
        handles += suggestion_handles()[: len(SUGGESTION_STYLES)]
    definitions = f"{DEFINITIONS}\n{SUGGESTION_DEFINITIONS}" if column == "region" else DEFINITIONS
    plot_grid(title, size_bins, series_for, info_for, handles, definitions, path)


def plot_distribution(table, allele_type, size_bins, regions, filters, measure, path):
    """NCR histograms of one measure, one per FILTER set, in one panel per size bin and region."""
    passing, total_column, unit, noun = DISTRIBUTION_MEASURES[measure]
    by_group = grouped(table, "filter", filters)
    colors = category_colors(filters)
    definitions = DISTRIBUTION_DEFINITIONS.format(noun=noun)
    definitions_height = DEFINITION_LINE_HEIGHT * (definitions.count("\n") + 2)
    legend_rows = -(-len(filters) // 3)
    bottom = MARGINS[1] + LEGEND_ROW_HEIGHT * (legend_rows - 1) + definitions_height
    height = PANEL_SIZE[1] * len(size_bins) + MARGINS[0] + bottom
    fig = plt.figure(figsize=(PANEL_SIZE[0] * len(regions) + 0.8, height), layout="none")
    grid = fig.add_gridspec(
        len(size_bins),
        len(regions),
        top=1 - MARGINS[0] / height,
        bottom=bottom / height,
        left=0.12,
        right=0.98,
        hspace=0.75,
        wspace=0.35,
    )
    drawn = set()
    for row_index, size_bin in enumerate(size_bins):
        for col_index, region in enumerate(regions):
            ax = fig.add_subplot(grid[row_index, col_index])
            pooled = rows_of(table, allele_type, size_bin, region=region)
            total = int(pooled[total_column].iloc[-1]) if len(pooled) else 0
            ax.set(xlim=(0, 1 + PAD), ylim=(0, 1 + PAD), xticks=TICKS, yticks=TICKS)
            ax.set_title(region, pad=15)
            draw_note(ax, f"{abbreviate(total)} {unit}", 1, 4)
            ax.set_xlabel("NCR")
            if col_index == 0:
                ax.set_ylabel(unit.capitalize())
                ax.text(
                    -0.42,
                    0.5,
                    "SNV" if size_bin == "0-0" else f"{size_bin} bp",
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="semibold",
                    color=INK,
                )
            if total == 0:
                ax.text(
                    0.5, 0.5, f"No {unit}", transform=ax.transAxes, ha="center", va="center", fontsize=8.5, color=INK_SECONDARY
                )
                continue
            # Every FILTER set shares the panel's bins, set from all of them together, so their bars line up
            cutoffs = pooled["cutoff"].to_numpy()
            steps = ncr_bin_steps(cutoffs, pooled[passing].to_numpy(), total)
            edges = list(range(0, len(cutoffs), steps))
            if edges[-1] != len(cutoffs) - 1:
                edges.append(len(cutoffs) - 1)
            highest = 0
            for label in filters:
                rows = rows_of(by_group, allele_type, size_bin, filter_value=label, region=region)
                if rows.empty or rows["sites_total"].iloc[-1] < MIN_HISTOGRAM_SITES:
                    continue
                drawn.add(label)
                # Each bin holds what passes its upper edge but not its lower one; the first also holds NCR = 0
                upper = rows[passing].to_numpy()[edges[1:]]
                counts = upper - np.r_[0, upper[:-1]]
                ax.stairs(counts, cutoffs[edges], color=colors[label], linewidth=1.4, zorder=2)
                highest = max(highest, counts.max())
            draw_note(ax, f"bin = {round(steps * (cutoffs[1] - cutoffs[0]), 6):g}", 0, -NOTE_BELOW_XLABEL)
            # Empty bins fall below the log axis and are clipped to its floor, under a single site or call
            if highest > 0:
                ax.set_yscale("log", nonpositive="clip")
                ax.set_ylim(0.5, highest * 2)
                ax.yaxis.set_minor_locator(NullLocator())

    def type_total(label):
        """The FILTER set's sites or calls across every size bin and region of the allele type."""
        return abbreviate(rows_of(by_group, allele_type, "all", filter_value=label)[total_column].iloc[-1])

    handles = [
        Line2D([], [], color=colors[label], linewidth=1.6, label=f"{label} ({type_total(label)} {unit})")
        for label in filters
        if label in drawn
    ]
    fig.suptitle(
        f"{allele_type.upper()} NCR Distribution of {unit.capitalize()} by Region and FILTER", fontsize=14, y=1 - 0.2 / height
    )
    # Allele types too sparse for any histogram draw no lines to label
    if handles:
        fig.legend(
            handles=handles, loc="lower center", bbox_to_anchor=(0.5, definitions_height / height), ncol=min(len(handles), 3)
        )
    fig.text(0.03, 0.1 / height, definitions, ha="left", va="bottom", fontsize=7.5, linespacing=1.5, color=INK_SECONDARY)
    fig.savefig(path)
    plt.close(fig)


def draw_note(ax, text, y, offset):
    """Small italic text centred on the panel, offset in points from its top (y=1) or bottom (y=0) edge."""
    transform = offset_copy(ax.transAxes, fig=ax.figure, y=offset, units="points")
    ax.text(
        0.5,
        y,
        text,
        transform=transform,
        ha="center",
        va="bottom" if offset > 0 else "top",
        fontsize=7.5,
        fontstyle="italic",
        color=INK_SECONDARY,
    )


def ncr_bin_steps(cutoffs, passing, total):
    """Cutoff steps per bin by NumPy's auto rule, with quartiles read off the cumulative counts."""
    shares = passing / total
    first_quartile, third_quartile = (cutoffs[np.searchsorted(shares, q)] for q in (0.25, 0.75))
    width = 1 / (np.log2(total) + 1)
    # Freedman-Diaconis is undefined when most of the counts share one NCR, leaving Sturges
    if third_quartile > first_quartile:
        width = min(width, 2 * (third_quartile - first_quartile) / total ** (1 / 3))
    return max(1, int(np.ceil(width / (cutoffs[1] - cutoffs[0]) - 1e-9)))


def abbreviate(count):
    """A count in thousands or millions to 1 dp, such as 5.2M or 106.7K."""
    for divisor, suffix in ((1e6, "M"), (1e3, "K")):
        if count >= divisor:
            return f"{count / divisor:.1f}{suffix}"
    return f"{count:,}"


def suggestion_handles():
    """Markers on the precision-recall and ROC curves, then the matching vertical lines on the cutoff plots."""
    styles = SUGGESTION_STYLES.values()
    markers = [
        Line2D(
            [],
            [],
            color=INK,
            marker=marker,
            linestyle="none",
            markerfacecolor="none",
            markersize=7,
            markeredgewidth=1.2,
            label=f"{label} (point)",
        )
        for label, marker, _ in styles
    ]
    lines = [
        Line2D([], [], color=INK_SECONDARY, linestyle=linestyle, linewidth=0.9, label=f"{label} (cutoff)")
        for label, _, linestyle in styles
    ]
    return markers + lines


def precision_axes(fig, spec, lowest):
    """Precision axes; when every value sits high, the y-axis is broken to keep 0 in view but zoom on the top."""
    floor = np.floor((lowest - ZOOM_STEP / 4) / ZOOM_STEP) * ZOOM_STEP
    # A row whose precision lines are all masked as sparse has no values to zoom on
    if not np.isfinite(floor) or floor < BREAK_ABOVE:
        ax = fig.add_subplot(spec)
        ax.set_ylim(0, 1 + PAD)
        return ax
    return brokenaxes(
        ylims=((0, STUB_HEIGHT), (floor, 1 + (1 - floor) * ZOOM_PAD)),
        height_ratios=[4, 1],
        hspace=0.12,
        subplot_spec=spec,
        fig=fig,
        diag_color=AXIS,
    )


def set_ticks(ax, broken):
    """Fixed 0-0.5-1 ticks; a broken axis gets three ticks on its zoomed segment and 0 on its stub."""
    # Runs after drawing because every brokenaxes plotting call resets its tick locators
    if not broken:
        ax.set(xticks=TICKS, yticks=TICKS)
        return
    top, bottom = ax.axs
    for sub in ax.axs:
        sub.set_xticks(TICKS)
    top.set_yticks(np.linspace(top.get_ylim()[0], 1, 3))
    bottom.set_yticks([0])


def draw_series(ax, series, x, y):
    """A curve, plus a hollow marker per suggestion on curve views or a vertical line per suggestion on cutoff views."""
    rows = series.rows
    ax.plot(rows[x], rows[y], color=series.color, linestyle=series.linestyle, linewidth=1.6, zorder=2)
    for name, index in (series.suggestions or {}).items():
        _, marker, linestyle = SUGGESTION_STYLES[name]
        if x == "cutoff" and series.cutoff_lines:
            ax.axvline(rows.loc[index, "cutoff"], color=INK_SECONDARY, linewidth=0.9, linestyle=linestyle, zorder=1)
        else:
            # A breakdown's markers take its line colour, so each region's cutoffs read against its own line
            edge = INK if series.cutoff_lines else series.color
            ax.plot(
                rows.loc[index, x],
                rows.loc[index, y],
                marker=marker,
                markersize=7,
                markerfacecolor=SURFACE,
                markeredgecolor=edge,
                markeredgewidth=1.4,
                linestyle="none",
                zorder=4,
            )


def draw_label(ax, size_bin):
    """The row's size bin; SNVs are the only calls in the 0-0 bin."""
    ax.set_axis_off()
    ax.text(
        0.05,
        0.62,
        "SNV" if size_bin == "0-0" else f"{size_bin} bp",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        fontweight="semibold",
        color=INK,
    )


def draw_info(ax, rows, size_bin):
    """Size bin, suggested cutoffs and the precision, FPR and recall denominators beside the row's plots."""
    draw_label(ax, size_bin)
    suggestions = suggest(rows)
    lines = [
        f"{label} = {rows.loc[suggestions[name], 'cutoff']:g}" if name in suggestions else f"{label} = n/a"
        for name, (label, _, _) in SUGGESTION_STYLES.items()
    ]
    classified, inherited = classified_calls(rows)
    lines += [
        f"n (precision) = {classified:,}",
        f"n (de novo) = {classified - inherited:,}",
        f"n (recall) = {int(rows['rec_denominator'].iloc[-1]):,}",
    ]
    ax.text(
        0.05,
        0.56,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
        linespacing=1.6,
        color=INK_SECONDARY,
    )


def cutoff_table(table, regions):
    """Per allele type, size bin and region of PASS sites: site and call counts, then each method's cutoff and failures.

    The table holds PASS rows only. A stratum with too few classified calls gets no suggestions, and one without
    de novo calls gets no knee or Youden J, since its precision is flat and its FPR undefined.
    """
    records = []
    for allele_type, size_bins in strata(table).items():
        for size_bin in size_bins:
            for region in ["all", *regions]:
                rows = select(table, allele_type, size_bin, region=region)
                if rows.empty:
                    continue
                sites, calls = (int(rows[column].iloc[0]) for column in ("sites_total", "rec_denominator"))
                classified, inherited = classified_calls(rows)
                record = {
                    "allele_type": allele_type,
                    "size_bin": size_bin,
                    "region": region,
                    "input_pass_sites": sites,
                    "input_pass_calls": calls,
                    "classified_proband_calls": classified,
                    "de_novo_proband_calls": classified - inherited,
                    "skipped_proband_calls": int(rows["prec_skipped"].iloc[-1]),
                    # The loosest cutoff keeps every site with NCR, so it is the callset before filtering
                    "baseline_precision": rows["precision"].iloc[-1],
                    "baseline_recall": rows["recall"].iloc[-1],
                }
                for method, index in suggest(rows).items():
                    row = rows.loc[index]
                    failed_sites = sites - int(row["sites_passing"])
                    failed_calls = calls - int(row["rec_numerator"])
                    record.update(
                        {
                            f"{method}_cutoff": row["cutoff"],
                            f"{method}_precision": row["precision"],
                            f"{method}_recall": row["recall"],
                            f"{method}_failed_sites": failed_sites,
                            f"{method}_failed_sites_proportion": failed_sites / sites,
                            f"{method}_failed_calls": failed_calls,
                            f"{method}_failed_calls_proportion": failed_calls / calls,
                            f"{method}_delta_precision": relative_change(row["precision"], record["baseline_precision"]),
                            f"{method}_delta_recall": relative_change(row["recall"], record["baseline_recall"]),
                        }
                    )
                records.append(record)
    columns = [
        "allele_type",
        "size_bin",
        "region",
        "input_pass_sites",
        "input_pass_calls",
        "classified_proband_calls",
        "de_novo_proband_calls",
        "skipped_proband_calls",
        "baseline_precision",
        "baseline_recall",
    ]
    columns += [
        f"{method}_{metric}"
        for method in SUGGESTION_STYLES
        for metric in (
            "cutoff",
            "precision",
            "recall",
            "failed_sites",
            "failed_sites_proportion",
            "failed_calls",
            "failed_calls_proportion",
            "delta_precision",
            "delta_recall",
        )
    ]
    return pd.DataFrame(records, columns=columns)


def length_bounds(size_bin):
    """Inclusive length bounds of a size bin label such as 0-0, 50-499 or 500+, as FilterNCRCutoffs reads them."""
    if size_bin.endswith("+"):
        return int(size_bin[:-1]), MAX_LENGTH
    low, high = size_bin.split("-")
    return int(low), int(high)


def recommended_cutoffs(table, regions, method, fallback_cutoff):
    """One cutoff per PASS allele type, size bin and region, from the region, then every region, then the fallback.

    Every stratum with PASS sites gets a row per region holding its sites, so the rows never overlap. Sites missing
    allele_length have no length bounds and get no row.
    """
    present = table[
        (table["allele_type"] != "all")
        & (~table["size_bin"].isin(["all", "missing"]))
        & (table["filter"] == "all")
        & (table["sites_total"] > 0)
    ]
    records, left_out = [], []
    for (allele_type, size_bin), _ in present[present["region"] == "all"].groupby(["allele_type", "size_bin"], sort=False):
        pooled = select(table, allele_type, size_bin)
        pooled_index = suggest(pooled).get(method) if len(pooled) else None
        for region in regions:
            if not (
                (present["allele_type"] == allele_type) & (present["size_bin"] == size_bin) & (present["region"] == region)
            ).any():
                continue
            rows = select(table, allele_type, size_bin, region=region)
            index = suggest(rows).get(method) if len(rows) else None
            if index is not None:
                cutoff, source = rows.loc[index, "cutoff"], "region"
            elif pooled_index is not None:
                cutoff, source = pooled.loc[pooled_index, "cutoff"], "all_regions"
            elif fallback_cutoff is not None:
                cutoff, source = fallback_cutoff, "fallback"
            else:
                left_out.append(f"{allele_type} {size_bin} {region}")
                continue
            min_length, max_length = length_bounds(size_bin)
            records.append(
                {
                    "allele_type": allele_type,
                    "min_length": min_length,
                    "max_length": max_length,
                    "region": region,
                    "cutoff": f"{cutoff:g}",
                    "method": method,
                    "source": source,
                }
            )
    if left_out:
        print(f"No {method} cutoff and no fallback for: {', '.join(left_out)}")
    return pd.DataFrame(records, columns=["allele_type", "min_length", "max_length", "region", "cutoff", "method", "source"])


def relative_change(value, baseline):
    """Signed percentage change from the unfiltered baseline, such as +0.22% or -13.54%; + is an improvement."""
    # Adding 0.0 turns a rounded -0.0 into 0.0, so no change reads +0.00%
    return f"{round((value - baseline) / baseline * 100, 2) + 0.0:+.2f}%"


def classified_calls(rows):
    """Classified and inherited proband calls at the loosest cutoff, where every call at a site with NCR passes."""
    return tuple(int(rows[column].iloc[-1]) for column in ("prec_denominator", "prec_numerator"))


def main():
    args = parse_args()
    apply_style()
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    raw, missing_ncr = load(args.tsv)
    table = derive(raw)
    passed = pass_only(table)
    filters, regions = filter_groups(table), region_groups(passed)

    # Only the FILTER breakdown and the distributions keep non-PASS sites
    for allele_type, size_bins in strata(passed).items():
        name = os.path.join(output_dir, allele_type.upper())
        plot_allele_type(passed, allele_type, size_bins, f"{name}.combined.png")
        plot_breakdown(passed, allele_type, size_bins, "region", regions, f"{name}.by_region.png")
    for allele_type, size_bins in strata(table).items():
        name = os.path.join(output_dir, allele_type.upper())
        plot_breakdown(table, allele_type, size_bins, "filter", filters, f"{name}.by_filter.png")
        for measure in DISTRIBUTION_MEASURES:
            plot_distribution(
                table, allele_type, size_bins, list(regions), filters, measure, f"{name}.distribution_{measure}.png"
            )

    cutoffs = cutoff_table(passed, regions)
    cutoffs.to_csv(os.path.join(output_dir, "cutoffs.tsv"), sep="\t", index=False, float_format="%.6f", na_rep="-")
    print(cutoffs.to_string(index=False))
    overall = missing_ncr[(missing_ncr[KEYS] == "all").all(axis=1)]
    for _, row in overall.iterrows():
        print(f"{row['rec_numerator']} calls at {row['sites_total']} sites without NCR are outside every curve")
    if args.recommended_tsv:
        recommended = recommended_cutoffs(passed, list(regions), args.method, args.fallback_cutoff)
        recommended.to_csv(args.recommended_tsv, sep="\t", index=False)
        print(f"Wrote {len(recommended)} recommended cutoffs to {args.recommended_tsv}")
    print(f"Wrote plots and cutoffs.tsv to {output_dir}")


if __name__ == "__main__":
    main()
