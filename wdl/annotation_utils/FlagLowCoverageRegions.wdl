version 1.0

import "../utils/Structs.wdl"

workflow FlagLowCoverageRegions {
    input {
        Array[File] mosdepth_bed_files
        Array[String] sample_ids
        File ped
        String prefix

        Int bin_size
        Float median_coverage_lower_threshold = 0.2

        String utils_docker

        RuntimeAttr? runtime_attr_bin_coverage
        RuntimeAttr? runtime_attr_aggregate_coverage
    }

    Array[Pair[File, String]] sample_inputs = zip(mosdepth_bed_files, sample_ids)

    scatter (sample_input in sample_inputs) {
        call BinSampleCoverage {
            input:
                mosdepth_bed = sample_input.left,
                ped = ped,
                sample_id = sample_input.right,
                prefix = "~{prefix}.~{sample_input.right}",
                bin_size = bin_size,
                median_coverage_lower_threshold = median_coverage_lower_threshold,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bin_coverage
        }
    }

    call AggregateLowCoverage {
        input:
            low_coverage_beds = BinSampleCoverage.low_coverage_bed,
            chromosome_coverage = mosdepth_bed_files[0],
            sample_histograms = BinSampleCoverage.coverage_histogram,
            prefix = "~{prefix}.cohort",
            bin_size = bin_size,
            sample_count = length(sample_inputs),
            docker = utils_docker,
            runtime_attr_override = runtime_attr_aggregate_coverage
    }

    output {
        Array[File] binned_coverage_tsvs = BinSampleCoverage.binned_coverage_tsv
        Array[File] low_coverage_beds = BinSampleCoverage.low_coverage_bed
        Array[File] sample_cutoff_tsvs = BinSampleCoverage.cutoff_tsv
        File sample_histograms_tar = AggregateLowCoverage.sample_histograms_tar
        File chromosome_low_coverage_tar = AggregateLowCoverage.chromosome_low_coverage_tar
        File cohort_low_coverage_counts = AggregateLowCoverage.cohort_low_coverage_counts
    }
}

task BinSampleCoverage {
    input {
        File mosdepth_bed
        File ped
        String sample_id
        String prefix
        Int bin_size
        Float median_coverage_lower_threshold
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv
import gzip
import math
import re
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np


INPUT = "~{mosdepth_bed}"
PED = "~{ped}"
SAMPLE_ID = "~{sample_id}"
BIN_SIZE = ~{bin_size}
MEDIAN_COVERAGE_LOWER_THRESHOLD = ~{median_coverage_lower_threshold}
PREFIX = "~{prefix}"
CHR_X = "chrX"
CHR_Y = "chrY"

plt.switch_backend("Agg")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def weighted_median(depth_counts, base_count):
    lower_rank = (base_count + 1) // 2
    upper_rank = (base_count + 2) // 2
    cumulative = 0
    lower_value = None
    upper_value = None
    for depth in sorted(depth_counts):
        cumulative += depth_counts[depth]
        if lower_value is None and cumulative >= lower_rank:
            lower_value = depth
        if cumulative >= upper_rank:
            upper_value = depth
            break
    return (lower_value + upper_value) / 2


def format_depth(depth):
    if depth.is_integer():
        return str(int(depth))
    return f"{depth:.10g}"


def read_sample_sex(path, sample_id):
    sex = None
    with open(path, "r") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"PED line {line_number} has fewer than 5 columns")
            if fields[1] != sample_id:
                continue
            if fields[4] == "1":
                parsed_sex = "male"
            elif fields[4] == "2":
                parsed_sex = "female"
            else:
                raise ValueError(
                    f"Sample {sample_id} has unsupported PED sex code {fields[4]}"
                )
            if sex is not None and sex != parsed_sex:
                raise ValueError(f"Sample {sample_id} has conflicting PED sex entries")
            sex = parsed_sex
    if sex is None:
        raise ValueError(f"Sample {sample_id} is missing from PED")
    return sex


def write_binned_coverage(input_path, output_path):
    coverage_counts = Counter()
    histogram = Counter()
    current_chrom = None
    current_start = 0
    current_end = BIN_SIZE
    covered_bases = 0
    previous_end = None

    def finish_bin(writer, end):
        nonlocal coverage_counts, covered_bases
        if covered_bases == 0:
            return
        median = weighted_median(coverage_counts, covered_bases)
        writer.writerow([current_chrom, current_start, end, format_depth(median)])
        histogram[median] += 1
        coverage_counts = Counter()
        covered_bases = 0

    with open_text(input_path) as source, gzip.open(output_path, "wt", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"Expected at least 4 columns at line {line_number}")
            chrom, start_text, end_text, depth_text = fields[:4]
            if SEX == "female" and chrom == CHR_Y:
                continue
            start = int(start_text)
            end = int(end_text)
            depth = float(depth_text)
            if start < 0 or end <= start or depth < 0:
                raise ValueError(f"Invalid interval at line {line_number}: {line.rstrip()}")

            if chrom != current_chrom:
                if current_chrom is not None:
                    finish_bin(writer, previous_end)
                if start != 0:
                    raise ValueError(f"Chromosome {chrom} starts at {start}, not 0")
                current_chrom = chrom
                current_start = 0
                current_end = BIN_SIZE
                previous_end = 0
            elif start != previous_end:
                raise ValueError(
                    f"Coverage intervals are not contiguous at line {line_number}: "
                    f"expected {chrom}:{previous_end}, found {chrom}:{start}"
                )

            position = start
            while position < end:
                overlap_end = min(end, current_end)
                overlap = overlap_end - position
                coverage_counts[depth] += overlap
                covered_bases += overlap
                position = overlap_end
                if position == current_end:
                    finish_bin(writer, current_end)
                    current_start = current_end
                    current_end += BIN_SIZE
            previous_end = end

        if current_chrom is None:
            raise ValueError("Coverage input contains no records")
        finish_bin(writer, previous_end)
    return histogram


def weighted_quantile(histogram, fraction):
    bin_count = sum(histogram.values())
    target_rank = max(1, math.ceil(fraction * bin_count))
    cumulative = 0
    for depth in sorted(histogram):
        cumulative += histogram[depth]
        if cumulative >= target_rank:
            return depth, bin_count
    raise ValueError("Cannot calculate cutoff from empty histogram")


def write_low_coverage_bins(binned_path, output_path, regular_cutoff, sex_cutoff):
    low_bin_count = 0
    low_histogram = Counter()
    with gzip.open(binned_path, "rt") as source, gzip.open(
        output_path, "wt", newline=""
    ) as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        for line in source:
            chrom, start, end, depth_text = line.rstrip("\n").split("\t")
            depth = float(depth_text)
            cutoff = (
                sex_cutoff
                if SEX == "male" and chrom in {CHR_X, CHR_Y}
                else regular_cutoff
            )
            if depth <= cutoff:
                writer.writerow([chrom, start, end, depth_text, SAMPLE_ID])
                low_bin_count += 1
                low_histogram[depth] += 1
    return low_bin_count, low_histogram


def plot_histogram(
    histogram,
    low_histogram,
    median,
    regular_cutoff,
    sex_cutoff,
    output_path,
):
    q1, bin_count = weighted_quantile(histogram, 0.25)
    q3, _ = weighted_quantile(histogram, 0.75)
    q95, _ = weighted_quantile(histogram, 0.95)
    interquartile_range = q3 - q1
    upper_fence = q3 + 3 * interquartile_range
    maximum_depth = max(histogram)
    display_maximum = min(maximum_depth, max(q95, upper_fence))
    if display_maximum == 0:
        display_maximum = 1

    plotted_depths = [depth for depth in sorted(histogram) if depth <= display_maximum]
    plotted_counts = [histogram[depth] for depth in plotted_depths]
    omitted_count = sum(
        count for depth, count in histogram.items() if depth > display_maximum
    )
    if interquartile_range > 0:
        bin_width = 2 * interquartile_range / bin_count ** (1 / 3)
        desired_bin_count = math.ceil(display_maximum / bin_width)
    else:
        desired_bin_count = math.ceil(math.log2(bin_count) + 1)
    histogram_bin_count = min(
        len(plotted_depths), max(10, min(100, desired_bin_count))
    )
    counts, edges = np.histogram(
        plotted_depths,
        bins=histogram_bin_count,
        range=(0, display_maximum),
        weights=plotted_counts,
    )
    low_depths = [depth for depth in plotted_depths if depth in low_histogram]
    low_counts, _ = np.histogram(
        low_depths,
        bins=edges,
        weights=[low_histogram[depth] for depth in low_depths],
    )
    figure, axis = plt.subplots(figsize=(10, 6))
    axis.axvspan(0, regular_cutoff, color="#fdd0a2", alpha=0.25)
    if SEX == "male":
        axis.axvspan(0, sex_cutoff, color="#f16913", alpha=0.15)
    axis.bar(
        edges[:-1],
        counts,
        width=np.diff(edges),
        align="edge",
        color="#4c78a8",
        edgecolor="none",
    )
    axis.bar(
        edges[:-1],
        low_counts,
        width=np.diff(edges),
        align="edge",
        color="#d95f02",
        edgecolor="none",
    )
    axis.axvline(regular_cutoff, color="#8c2d04", linestyle="--", linewidth=1.5)
    if SEX == "male":
        axis.axvline(sex_cutoff, color="#d94801", linestyle=":", linewidth=1.5)
    axis.axvline(median, color="#525252", linestyle="-.", linewidth=1.2)
    axis.set_xlim(0, display_maximum)
    axis.set_xlabel("Median Coverage")
    axis.set_ylabel("Bin Count")
    axis.set_title(SAMPLE_ID)
    annotation = f"Median: {median:.2f}×\nCoverage cutoff: {regular_cutoff:.2f}×"
    if SEX == "male":
        annotation += f"\nchrX/chrY cutoff: {sex_cutoff:.2f}×"
    else:
        annotation += "\nchrY excluded"
    if omitted_count:
        annotation += (
            f"\n{omitted_count:,} bins above {format_depth(display_maximum)}× not shown"
        )
    axis.text(
        0.98,
        0.95,
        annotation,
        horizontalalignment="right",
        verticalalignment="top",
        transform=axis.transAxes,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


if BIN_SIZE <= 0:
    raise ValueError("bin_size must be greater than 0")
if not 0 < MEDIAN_COVERAGE_LOWER_THRESHOLD <= 1:
    raise ValueError("median_coverage_lower_threshold must be greater than 0 and at most 1")
if not re.fullmatch(r"[A-Za-z0-9._-]+", SAMPLE_ID):
    raise ValueError("sample_id may contain only letters, numbers, '.', '_', and '-'")

SEX = read_sample_sex(PED, SAMPLE_ID)
binned_path = f"{PREFIX}.binned_coverage.tsv.gz"
low_path = f"{PREFIX}.low_coverage.bed.gz"
cutoff_path = f"{PREFIX}.low_coverage_cutoff.tsv"
plot_path = f"{PREFIX}.coverage_histogram.png"

histogram = write_binned_coverage(INPUT, binned_path)
median, bin_count = weighted_quantile(histogram, 0.5)
regular_cutoff = median * MEDIAN_COVERAGE_LOWER_THRESHOLD
sex_cutoff = regular_cutoff / 2
low_bin_count, low_histogram = write_low_coverage_bins(
    binned_path, low_path, regular_cutoff, sex_cutoff
)
with open(cutoff_path, "w", newline="") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "sample_id",
            "sex",
            "median",
            "median_coverage_lower_threshold",
            "regular_cutoff",
            "chrX_cutoff",
            "chrY_cutoff",
            "bin_count",
            "low_bin_count",
        ]
    )
    writer.writerow(
        [
            SAMPLE_ID,
            SEX,
            format_depth(median),
            MEDIAN_COVERAGE_LOWER_THRESHOLD,
            format_depth(regular_cutoff),
            format_depth(sex_cutoff) if SEX == "male" else format_depth(regular_cutoff),
            format_depth(sex_cutoff) if SEX == "male" else "excluded",
            bin_count,
            low_bin_count,
        ]
    )
plot_histogram(
    histogram,
    low_histogram,
    median,
    regular_cutoff,
    sex_cutoff,
    plot_path,
)
PYCODE
    >>>

    output {
        File binned_coverage_tsv = "~{prefix}.binned_coverage.tsv.gz"
        File low_coverage_bed = "~{prefix}.low_coverage.bed.gz"
        File cutoff_tsv = "~{prefix}.low_coverage_cutoff.tsv"
        File coverage_histogram = "~{prefix}.coverage_histogram.png"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(mosdepth_bed, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task AggregateLowCoverage {
    input {
        Array[File] low_coverage_beds
        File chromosome_coverage
        Array[File] sample_histograms
        String prefix
        Int bin_size
        Int sample_count
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir chromosome_plots sample_histograms
        printf '%s\n' ~{sep=' ' low_coverage_beds} | xargs -r gzip -dc \
            | cut -f1-3 \
            | LC_ALL=C sort --parallel=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} -k1,1V -k2,2n -k3,3n \
            | uniq -c \
            | awk 'BEGIN {OFS="\t"; print "chrom", "start", "end", "sample_count"} {print $2, $3, $4, $1}' \
            | gzip > ~{prefix}.low_coverage_counts.tsv.gz

        python3 <<'PYCODE'
import csv
import gzip
import re
from pathlib import Path

import matplotlib.pyplot as plt


COUNTS = "~{prefix}.low_coverage_counts.tsv.gz"
CHROMOSOME_COVERAGE = "~{chromosome_coverage}"
BIN_SIZE = ~{bin_size}
SAMPLE_COUNT = ~{sample_count}
OUTPUT_DIR = Path("chromosome_plots")

plt.switch_backend("Agg")


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def safe_filename(value):
    return re.sub(r"[^A-Za-z0-9._-]", "_", value)


def plot_chromosome(chrom, chrom_end, flagged_bins):
    starts = range(0, chrom_end, BIN_SIZE)
    positions = [
        (start + min(start + BIN_SIZE, chrom_end)) / 2_000_000 for start in starts
    ]
    counts = [flagged_bins.get(start, 0) for start in starts]
    figure, axis = plt.subplots(figsize=(14, 5))
    axis.fill_between(positions, counts, step="mid", color="#4c78a8", alpha=0.8)
    axis.plot(
        positions,
        counts,
        color="#2c5282",
        linewidth=0.6,
        drawstyle="steps-mid",
    )
    axis.set_xlabel("Position")
    axis.set_ylabel("Low Coverage Sample Count")
    axis.set_ylim(0, SAMPLE_COUNT)
    axis.set_title(chrom)
    figure.tight_layout()
    figure.savefig(OUTPUT_DIR / f"{safe_filename(chrom)}.low_coverage.png", dpi=150)
    plt.close(figure)


if SAMPLE_COUNT <= 0:
    raise ValueError("sample_count must be greater than 0")

chromosome_ends = {}
with open_text(CHROMOSOME_COVERAGE) as source:
    for line in source:
        if not line.strip() or line.startswith("#"):
            continue
        chrom, _, end_text, _ = line.rstrip("\n").split("\t")[:4]
        chromosome_ends[chrom] = int(end_text)

flagged_bins = {}
with gzip.open(COUNTS, "rt") as source:
    reader = csv.reader(source, delimiter="\t")
    header = next(reader, None)
    if header != ["chrom", "start", "end", "sample_count"]:
        raise ValueError("Unexpected cohort counts header")
    for chrom, start_text, end_text, count_text in reader:
        if chrom not in chromosome_ends:
            raise ValueError(f"Low-coverage bin uses unknown chromosome {chrom}")
        start = int(start_text)
        expected_end = min(start + BIN_SIZE, chromosome_ends[chrom])
        if int(end_text) != expected_end:
            raise ValueError(f"Low-coverage bin has unexpected coordinates: {chrom}:{start}")
        flagged_bins.setdefault(chrom, {})[start] = int(count_text)

for chrom, chrom_end in chromosome_ends.items():
    plot_chromosome(chrom, chrom_end, flagged_bins.get(chrom, {}))
PYCODE

        cp ~{sep=' ' sample_histograms} sample_histograms/
        tar -czf ~{prefix}.sample_histograms.tar.gz -C sample_histograms .
        tar -czf ~{prefix}.chromosome_low_coverage_plots.tar.gz -C chromosome_plots .
    >>>

    output {
        File sample_histograms_tar = "~{prefix}.sample_histograms.tar.gz"
        File chromosome_low_coverage_tar = "~{prefix}.chromosome_low_coverage_plots.tar.gz"
        File cohort_low_coverage_counts = "~{prefix}.low_coverage_counts.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 4 * ceil(size(low_coverage_beds, "GB")) + ceil(size(chromosome_coverage, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
