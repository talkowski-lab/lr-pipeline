version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SummarizeSingletonCalls {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        Int min_length = 50

        Int? records_per_shard
        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_count
        RuntimeAttr? runtime_attr_merge
    }

    scatter (i in range(length(vcfs))) {
        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.input_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [vcfs[i]]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [vcf_idxs[i]]])

        scatter (j in range(length(shard_vcfs))) {
            call Helpers.SubsetVcfByLength {
                input:
                    vcf = shard_vcfs[j],
                    vcf_idx = shard_vcf_idxs[j],
                    min_length = min_length,
                    prefix = "~{prefix}.input_~{i}.shard_~{j}.subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }

            call CountSingletonShard {
                input:
                    vcf = SubsetVcfByLength.subset_vcf,
                    vcf_idx = SubsetVcfByLength.subset_vcf_idx,
                    prefix = "~{prefix}.input_~{i}.shard_~{j}.counts",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_count
            }
        }
    }

    call MergeSingletonCounts {
        input:
            count_tsvs = flatten(CountSingletonShard.counts_tsv),
            sample_files = flatten(CountSingletonShard.samples),
            prefix = "~{prefix}.singleton_diagnostics",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File singleton_counts_tsv = MergeSingletonCounts.merged_counts_tsv
    }
}

task CountSingletonShard {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf} > ~{prefix}.samples.txt
        bcftools query \
            -f '%INFO/allele_type\t%INFO/allele_length\t%INFO/REGION[\t%GT\t%EV]\n' \
            ~{vcf} > variants.tsv

        python3 <<'PYCODE'
import csv
import re
from collections import defaultdict


INPUT = "variants.tsv"
SAMPLES_INPUT = "~{prefix}.samples.txt"
OUTPUT = "~{prefix}.tsv"


def get_size_range(length):
    if length < 50:
        return "<50"
    if length < 200:
        return "50-199"
    if length < 500:
        return "200-499"
    if length < 5000:
        return "500-4999"
    return "5000+"


def get_count_type(ev):
    callers = set()
    if ev not in {"", "."}:
        callers = {value.split("_(", 1)[0].lower() for value in ev.split(",")}
    if callers - {"hapdiff", "dipcall"}:
        return "alignments"
    if callers:
        return "assemblies"
    return "other"


def get_singleton_type(gt):
    alleles = re.split(r"[/|]", gt)
    if not alleles or "." in alleles:
        return None
    alt_count = sum(int(allele) > 0 for allele in alleles)
    if alt_count == 0:
        return "ac_0"
    if alt_count == 1:
        return "ac_1"
    return "ac_>1"


with open(SAMPLES_INPUT, "r") as handle:
    samples = [line.rstrip("\n") for line in handle]

counts = defaultdict(int)
with open(INPUT, "r") as handle:
    for line_number, line in enumerate(handle, 1):
        fields = line.rstrip("\n").split("\t")
        expected_fields = 3 + 2 * len(samples)
        if len(fields) != expected_fields:
            raise ValueError(
                f"Expected {expected_fields} fields at input line {line_number}, found {len(fields)}"
            )

        variant_type, allele_length, region = fields[:3]
        size_range = get_size_range(abs(int(allele_length)))
        for index, sample in enumerate(samples):
            gt = fields[3 + 2 * index]
            ev = fields[4 + 2 * index]
            singleton_type = get_singleton_type(gt)
            if singleton_type is None:
                continue
            count_type = "other" if singleton_type == "ac_0" else get_count_type(ev)
            key = (sample, variant_type, size_range, region, count_type, singleton_type)
            counts[key] += 1

with open(OUTPUT, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "sample",
        "variant_type",
        "size_range",
        "region",
        "count_type",
        "singleton_type",
        "count",
    ])
    for key in sorted(counts):
        writer.writerow([*key, counts[key]])
PYCODE
    >>>

    output {
        File counts_tsv = "~{prefix}.tsv"
        File samples = "~{prefix}.samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size([vcf, vcf_idx], "GB")) + 10,
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

task MergeSingletonCounts {
    input {
        Array[File] count_tsvs
        Array[File] sample_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv
from collections import defaultdict


COUNT_FILES = [path for path in "~{sep=',' count_tsvs}".split(",") if path]
SAMPLE_FILES = [path for path in "~{sep=',' sample_files}".split(",") if path]
OUTPUT = "~{prefix}.tsv"
SIZE_RANGES = ["<50", "50-199", "200-499", "500-4999", "5000+"]
COUNT_TYPES = ["assemblies", "alignments", "other"]
SINGLETON_TYPES = ["ac_0", "ac_1", "ac_>1"]


samples = set()
for path in SAMPLE_FILES:
    with open(path, "r") as handle:
        samples.update(line.rstrip("\n") for line in handle)

counts = defaultdict(int)
variant_regions = set()
expected_header = [
    "sample",
    "variant_type",
    "size_range",
    "region",
    "count_type",
    "singleton_type",
    "count",
]
for path in COUNT_FILES:
    with open(path, "r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if header != expected_header:
            raise ValueError(f"Unexpected count table header in {path}: {header}")
        for row in reader:
            sample, variant_type, size_range, region, count_type, singleton_type, count = row
            counts[(sample, variant_type, size_range, region, count_type, singleton_type)] += int(count)
            variant_regions.add((variant_type, region))

ordered_strata = [
    (variant_type, size_range, region)
    for variant_type in sorted({value[0] for value in variant_regions})
    for size_range in SIZE_RANGES
    for region in sorted(value[1] for value in variant_regions if value[0] == variant_type)
]
columns = [
    (*stratum, count_type, singleton_type)
    for stratum in ordered_strata
    for count_type in COUNT_TYPES
    for singleton_type in SINGLETON_TYPES
]

with open(OUTPUT, "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    header = ["sample"] + [" - ".join(column) for column in columns]
    writer.writerow(header)
    for sample in sorted(samples):
        values = [counts[(sample, *column)] for column in columns]
        writer.writerow([sample, *values])
PYCODE
    >>>

    output {
        File merged_counts_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 4 * ceil(size(count_tsvs, "GB")) + 10,
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
