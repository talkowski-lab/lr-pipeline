#!/usr/bin/env python3
"""Sum per-contig EvaluateCallsetAccuracy count tables across contigs and plot precision and recall.

Each input TSV has one row per level (site or call), sample, allele type, length bin, region and matching stage,
with tp_eval, fp, tp_truth and fn counts. The counts are summed across contigs and precision, recall and F1
recomputed, giving --combined_tsv in the same layout; --per_contig_tsv stacks the inputs under a leading contig
column. Allele types, length bins, regions, samples and stages are read from the tables, in the order they first
appear. Figures in --output_dir, from the last matching stage unless noted:
  accuracy_sites.png, accuracy_calls.png - precision and recall per allele type and length bin, one bar pair per
      region, with the callset and truth counts under each length bin.
  accuracy_by_sample.png - call precision and recall per sample and allele type, across every length bin and region.
  accuracy_by_stage.png  - site precision and recall per allele type and length bin after each matching stage.
summary.tsv holds the site and call counts, precision, recall and F1 of every allele type, length bin and region
across samples at the last stage.
"""

import argparse
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
# Categorical slots in fixed order, so a region or stage keeps its colour across figures
CATEGORY_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]
ALL_COLOR = INK_SECONDARY
# Recall bars repeat their precision bar's colour, lightened
RECALL_ALPHA = 0.45

KEYS = ["level", "sample", "allele_type", "length_bin", "region", "match_stage"]
COUNT_COLUMNS = ["tp_eval", "fp", "tp_truth", "fn"]
METRICS = ["precision", "recall", "f1"]
PANEL_HEIGHT = 3.2
GROUP_WIDTH = 0.8
DEFINITIONS = (
    "•  Precision = callset sites (or calls) matching the truth / all callset sites (or calls); recall = truth sites "
    "(or calls) matched by the callset / all truth sites (or calls).\n"
    "•  A call is one sample's non-reference genotype; a call matches when the same sample is non-reference in a "
    "matched record on the other side. A bar is missing where its side has no records."
)
COUNTS_DEFINITION = "\n•  Counts under each length bin are callset / truth sites (or calls) across every region."


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tsvs", nargs="+", required=True, help="Per-contig accuracy TSVs of CountCallsetAccuracy")
    parser.add_argument("--contigs", nargs="+", required=True, help="Contig of each TSV, in the same order")
    parser.add_argument("--combined_tsv", required=True, help="Path for the counts summed across contigs")
    parser.add_argument("--per_contig_tsv", required=True, help="Path for the per-contig rows under a contig column")
    parser.add_argument("--output_dir", required=True, help="Directory for the figures and summary.tsv")
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
            "axes.titlesize": 10,
            "axes.titleweight": "semibold",
            "axes.labelsize": 9,
            "axes.grid": True,
            "axes.grid.axis": "y",
            "axes.axisbelow": True,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelcolor": INK_SECONDARY,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 8,
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "figure.titlesize": 12.5,
            "figure.titleweight": "semibold",
            "figure.dpi": 200,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.3,
        }
    )


def load(path):
    return pd.read_csv(path, sep="\t", dtype={key: str for key in KEYS}, keep_default_na=False)


def with_metrics(table):
    """Precision, recall and F1 from the counts; a stratum with no records on a side leaves its metric undefined."""
    table = table.copy()
    precision = table["tp_eval"] / (table["tp_eval"] + table["fp"]).replace(0, np.nan)
    recall = table["tp_truth"] / (table["tp_truth"] + table["fn"]).replace(0, np.nan)
    table["precision"], table["recall"] = precision, recall
    table["f1"] = 2 * precision * recall / (precision + recall)
    return table


def combine(tables, contigs):
    """Rows of every contig under a contig column, and their counts summed across contigs."""
    per_contig = pd.concat([table[KEYS + COUNT_COLUMNS].assign(contig=contig) for table, contig in zip(tables, contigs)])
    per_contig = per_contig[["contig"] + KEYS + COUNT_COLUMNS]
    combined = per_contig.groupby(KEYS, sort=False, as_index=False)[COUNT_COLUMNS].sum()
    return with_metrics(per_contig), with_metrics(combined)


def write_tsv(table, path):
    table.to_csv(path, sep="\t", index=False, float_format="%.4f", na_rep=".")


def ordered(values, last="all"):
    """Unique values in first-seen order, with the `all` marginal moved to the end or dropped when last is None."""
    unique = [value for value in pd.unique(values) if value != "all"]
    return unique + ([last] if last and "all" in set(values) else [])


def rows_of(table, **selection):
    mask = np.ones(len(table), dtype=bool)
    for column, value in selection.items():
        mask &= table[column].to_numpy() == value
    return table[mask]


def abbreviate(count):
    """A count in thousands or millions to 1 dp, such as 5.2M or 106.7K."""
    for divisor, suffix in ((1e6, "M"), (1e3, "K")):
        if count >= divisor:
            return f"{count / divisor:.1f}{suffix}"
    return f"{int(count):,}"


def series_colors(labels):
    colors = iter(CATEGORY_COLORS)
    return {label: ALL_COLOR if label == "all" else next(colors) for label in labels}


def draw_bars(ax, groups, series, value_of, colors):
    """Precision and recall bars per series in each group, with value labels; undefined values draw no bar."""
    width = GROUP_WIDTH / (2 * len(series))
    for group_index, group in enumerate(groups):
        for series_index, label in enumerate(series):
            for metric_index, metric in enumerate(("precision", "recall")):
                value = value_of(group, label, metric)
                if value is None or not np.isfinite(value):
                    continue
                x = group_index - GROUP_WIDTH / 2 + width * (2 * series_index + metric_index + 0.5)
                alpha = 1 if metric == "precision" else RECALL_ALPHA
                ax.bar(x, value, width * 0.92, color=colors[label], alpha=alpha, linewidth=0, zorder=2)
                ax.text(x, value + 0.01, f"{value:.3f}", ha="center", va="bottom", rotation=90, fontsize=5.5, color=INK)
    ax.set_xticks(range(len(groups)))
    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_ylim(0, 1.18)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_ylabel("Precision / Recall")


def legend_handles(series, colors, names=None):
    handles = []
    for label in series:
        name = (names or {}).get(label, label)
        handles.append(Patch(color=colors[label], label=f"{name} precision"))
        handles.append(Patch(color=colors[label], alpha=RECALL_ALPHA, label=f"{name} recall"))
    return handles


def type_figure(allele_types, groups_of, title, handles, group_width=2.2, definitions=DEFINITIONS):
    """One panel per allele type, group_width inches per group, with the legend and definitions below."""
    widths = [max(len(groups_of(allele_type)), 1) for allele_type in allele_types]
    fig, axes = plt.subplots(
        1,
        len(allele_types),
        figsize=(max(group_width * sum(widths), 6), PANEL_HEIGHT + 1.4),
        gridspec_kw={"width_ratios": widths},
        squeeze=False,
        layout="constrained",
    )
    fig.suptitle(title)
    fig.legend(handles=handles, loc="outside lower center", ncol=min(len(handles), 6))
    fig.text(0.01, -0.02, definitions, ha="left", va="top", fontsize=7, linespacing=1.5, color=INK_SECONDARY)
    return fig, axes[0]


def plot_by_region(table, level, stage, path):
    """Precision and recall per allele type and length bin, one bar pair per region."""
    rows = rows_of(table, level=level, sample="all", match_stage=stage)
    allele_types = ordered(rows["allele_type"], last=None)
    regions = ordered(rows["region"])
    colors = series_colors(regions)
    lookup = rows.set_index(["allele_type", "length_bin", "region"])

    def bins_of(allele_type):
        return ordered(rows_of(rows, allele_type=allele_type)["length_bin"], last=None)

    noun = "Site" if level == "site" else "Call"
    fig, axes = type_figure(
        allele_types,
        bins_of,
        f"{noun} Precision and Recall by Region",
        legend_handles(regions, colors, {"all": "All regions"}),
        definitions=DEFINITIONS + COUNTS_DEFINITION,
    )
    for ax, allele_type in zip(axes, allele_types):
        bins = bins_of(allele_type)

        def value_of(length_bin, region, metric, allele_type=allele_type):
            key = (allele_type, length_bin, region)
            return lookup.loc[key, metric] if key in lookup.index else None

        draw_bars(ax, bins, regions, value_of, colors)
        labels = []
        for length_bin in bins:
            totals = lookup.loc[(allele_type, length_bin, "all")]
            callset, truth = totals["tp_eval"] + totals["fp"], totals["tp_truth"] + totals["fn"]
            labels.append(f"{length_bin} bp\n{abbreviate(callset)} / {abbreviate(truth)}")
        ax.set_xticklabels(labels)
        ax.set_title(allele_type.upper())
    fig.savefig(path)
    plt.close(fig)


def plot_by_sample(table, stage, path):
    """Call precision and recall per sample and allele type, across every length bin and region."""
    rows = rows_of(table, level="call", length_bin="all", region="all", match_stage=stage)
    allele_types = ordered(rows["allele_type"], last=None)
    samples = ordered(rows["sample"])
    colors = {"calls": CATEGORY_COLORS[0]}
    lookup = rows.set_index(["allele_type", "sample"])
    fig, axes = type_figure(
        allele_types,
        lambda allele_type: samples,
        "Call Precision and Recall by Sample",
        legend_handles(["calls"], colors, {"calls": "Call"}),
        group_width=0.9,
    )
    for ax, allele_type in zip(axes, allele_types):

        def value_of(sample, _, metric, allele_type=allele_type):
            key = (allele_type, sample)
            return lookup.loc[key, metric] if key in lookup.index else None

        draw_bars(ax, samples, ["calls"], value_of, colors)
        ax.set_xticklabels(samples, rotation=45, ha="right")
        ax.set_title(allele_type.upper())
    fig.savefig(path)
    plt.close(fig)


def plot_by_stage(table, stages, path):
    """Site precision and recall per allele type and length bin after each cumulative matching stage."""
    rows = rows_of(table, level="site", sample="all", region="all")
    allele_types = ordered(rows["allele_type"], last=None)
    colors = series_colors(stages)
    lookup = rows.set_index(["allele_type", "length_bin", "match_stage"])

    def bins_of(allele_type):
        return ordered(rows_of(rows, allele_type=allele_type)["length_bin"], last=None)

    names = {stage: f"After {stage}" for stage in stages}
    fig, axes = type_figure(
        allele_types, bins_of, "Site Precision and Recall by Matching Stage", legend_handles(stages, colors, names)
    )
    for ax, allele_type in zip(axes, allele_types):
        bins = bins_of(allele_type)

        def value_of(length_bin, stage, metric, allele_type=allele_type):
            key = (allele_type, length_bin, stage)
            return lookup.loc[key, metric] if key in lookup.index else None

        draw_bars(ax, bins, stages, value_of, colors)
        ax.set_xticklabels([f"{length_bin} bp" for length_bin in bins])
        ax.set_title(allele_type.upper())
    fig.savefig(path)
    plt.close(fig)


def summary(table, stage):
    """Site and call counts and metrics side by side per allele type, length bin and region, across samples.

    Every site stratum has calls, since a kept record carries at least one non-reference genotype.
    """
    rows = rows_of(table, sample="all", match_stage=stage)
    wide = None
    for level in ("site", "call"):
        part = rows_of(rows, level=level)[["allele_type", "length_bin", "region"] + COUNT_COLUMNS + METRICS]
        part = part.rename(columns={column: f"{level}_{column}" for column in COUNT_COLUMNS + METRICS})
        wide = part if wide is None else wide.merge(part, on=["allele_type", "length_bin", "region"], how="left")
    return wide


def main():
    args = parse_args()
    if len(args.tsvs) != len(args.contigs):
        raise SystemExit(f"{len(args.tsvs)} TSVs but {len(args.contigs)} contigs")
    apply_style()
    os.makedirs(args.output_dir, exist_ok=True)
    per_contig, combined = combine([load(path) for path in args.tsvs], args.contigs)
    write_tsv(per_contig, args.per_contig_tsv)
    write_tsv(combined, args.combined_tsv)

    stages = ordered(combined["match_stage"], last=None)
    final = stages[-1]
    for level, name in (("site", "accuracy_sites.png"), ("call", "accuracy_calls.png")):
        plot_by_region(combined, level, final, os.path.join(args.output_dir, name))
    plot_by_sample(combined, final, os.path.join(args.output_dir, "accuracy_by_sample.png"))
    plot_by_stage(combined, stages, os.path.join(args.output_dir, "accuracy_by_stage.png"))
    write_tsv(summary(combined, final), os.path.join(args.output_dir, "summary.tsv"))

    overall = rows_of(combined, sample="all", allele_type="all", length_bin="all", region="all", match_stage=final)
    for row in overall.itertuples():
        print(f"{row.level}: precision {row.precision:.4f}, recall {row.recall:.4f}, F1 {row.f1:.4f}")
    print(f"Summed {len(args.tsvs)} contigs into {args.combined_tsv}; wrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
