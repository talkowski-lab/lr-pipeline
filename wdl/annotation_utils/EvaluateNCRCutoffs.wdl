version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow EvaluateNCRCutoffs {
    meta {
        description: [
            "This utility tabulates precision and recall of a callset across a sweep of `INFO/NCR` (no-call rate) cutoffs, so a cutoff can be chosen from the resulting precision-recall curve. A site passes a cutoff when its `NCR` is at or below it. Precision is estimated from trio transmission: trios whose child and both parents are VCF samples are taken from the PED, and each site where a trio child carries a non-reference allele is one observation, counted as transmitted when either parent also carries a non-reference allele and as de novo when both parents are homozygous reference. Phasing and zygosity are ignored, so `0/1`, `1/1` and `0|1` all count as carriers. Observations where the child or either parent has a missing allele are skipped, since treating missing as non-carrier would inflate the de novo count precisely at high-`NCR` sites. Precision is the transmitted fraction of observations. Recall is the fraction of sites, and separately of trio observations, retained relative to the full set after `args_string_vcf`. Sites without `INFO/NCR` are excluded from every row and their count is printed to the summary task log."
        ]
    }

    parameter_meta {
        vcf: "VCF carrying `INFO/NCR` to evaluate."
        vcf_idx: "Index for `vcf`."
        contig: "Contig to evaluate."
        records_per_shard: "Number of variants to keep within a single shard during evaluation."
        ped: "PED file used to identify the trios whose transmission estimates precision."
        args_string_vcf: "`bcftools view` include expression applied to the VCF before evaluation."
        cutoff_step: "Spacing of the `NCR` cutoffs evaluated from `0` to `1` inclusive."
        ncr_cutoffs_tsv: "One row per cutoff with columns `cutoff`, `n_sites` passing it, `recall_sites`, `n_trio_carrier_obs` (trio child carrier observations at passing sites), `recall_trio_carrier_obs`, `n_transmitted`, `n_denovo` and `precision`, the transmitted fraction of observations (empty when there are none)."
    }

    input {
        File vcf
        File vcf_idx
        String contig
        Int? records_per_shard
        String prefix

        File ped
        String? args_string_vcf
        Float cutoff_step = 0.01

        String utils_docker

        RuntimeAttr? runtime_attr_find_trios
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_count_trio_transmission
        RuntimeAttr? runtime_attr_summarize_ncr_cutoffs
    }

    call Helpers.FindTrios {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            ped = ped,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_find_trios
    }

    if (defined(records_per_shard)) {
        call Helpers.ShardVcfByRecords {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                records_per_shard = select_first([records_per_shard]),
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }
    }

    Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [vcf]])
    Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [vcf_idx]])

    scatter (shard_i in range(length(vcfs_to_process))) {
        call CountTrioTransmission {
            input:
                vcf = vcfs_to_process[shard_i],
                vcf_idx = vcf_idxs_to_process[shard_i],
                contig = contig,
                trio_definitions = FindTrios.trio_definitions,
                trio_sample_ids_file = FindTrios.trio_sample_ids_file,
                include_args = args_string_vcf,
                prefix = "~{prefix}.~{contig}.shard_~{shard_i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_count_trio_transmission
        }
    }

    call SummarizeNCRCutoffs {
        input:
            count_tsvs = CountTrioTransmission.ncr_counts_tsv,
            cutoff_step = cutoff_step,
            prefix = "~{prefix}.ncr_cutoffs",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_summarize_ncr_cutoffs
    }

    output {
        File ncr_cutoffs_tsv = SummarizeNCRCutoffs.ncr_cutoffs_tsv
    }
}

task CountTrioTransmission {
    input {
        File vcf
        File vcf_idx
        String contig
        File trio_definitions
        File trio_sample_ids_file
        String? include_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [[ "~{vcf_idx}" != "~{vcf}.tbi" ]]; then
            ln -sf "~{vcf_idx}" "~{vcf}.tbi"
        fi

        if [[ ! -s ~{trio_definitions} ]]; then
            echo "No trios with all three members in the VCF were found in the PED" >&2
            exit 1
        fi

        # Filter on all samples before subsetting to trio members, so sample-level expressions behave as in bcftools view
        bcftools view \
            -r ~{contig} \
            ~{if defined(include_args) then "-i '~{include_args}'" else ""} \
            -Ou \
            ~{vcf} \
            | bcftools query \
                -H \
                -S ~{trio_sample_ids_file} \
                -f '%INFO/NCR[\t%GT]\n' \
                > genotypes.tsv

        python3 <<'CODE'
import csv
import re
from collections import defaultdict

ALLELE_SPLIT = re.compile(r'[/|]')

def is_called(gt):
    return '.' not in gt

def is_carrier(gt):
    return any(allele != '0' for allele in ALLELE_SPLIT.split(gt))

with open("~{trio_definitions}") as handle:
    trios = [line.rstrip('\n').split('\t')[:3] for line in handle if line.strip()]

# Keyed by NCR exactly as printed, so the summary task can apply any cutoff step without rebinning error
counts = defaultdict(lambda: [0, 0, 0, 0])
with open('genotypes.tsv') as handle:
    header = next(handle).lstrip('#').rstrip('\n').split('\t')
    columns = {re.sub(r'^\s*\[\d+\]|:GT$', '', name): i for i, name in enumerate(header)}
    trio_columns = [(columns[child], columns[father], columns[mother]) for child, father, mother in trios]
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        row = counts[fields[0]]
        row[0] += 1
        for child_col, father_col, mother_col in trio_columns:
            child, father, mother = fields[child_col], fields[father_col], fields[mother_col]
            if not (is_called(child) and is_called(father) and is_called(mother)) or not is_carrier(child):
                continue
            row[1] += 1
            if is_carrier(father) or is_carrier(mother):
                row[2] += 1
            else:
                row[3] += 1

with open("~{prefix}.ncr_counts.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['ncr', 'n_sites', 'n_trio_carrier_obs', 'n_transmitted', 'n_denovo'])
    for ncr, row in counts.items():
        writer.writerow([ncr, *row])

print(f'Counted {sum(row[0] for row in counts.values())} sites across {len(trios)} trios')
CODE
    >>>

    output {
        File ncr_counts_tsv = "~{prefix}.ncr_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
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

task SummarizeNCRCutoffs {
    input {
        Array[File] count_tsvs
        Float cutoff_step
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import csv
from collections import defaultdict

COUNT_COLUMNS = ['n_sites', 'n_trio_carrier_obs', 'n_transmitted', 'n_denovo']
STEP = ~{cutoff_step}
EPSILON = 1e-9

def ratio(numerator, denominator):
    return f'{numerator / denominator:.6f}' if denominator else ''

totals = defaultdict(lambda: [0, 0, 0, 0])
missing_ncr_sites = 0
for path in "~{sep=',' count_tsvs}".split(','):
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['ncr'] == '.':
                missing_ncr_sites += int(row['n_sites'])
                continue
            total = totals[float(row['ncr'])]
            for i, column in enumerate(COUNT_COLUMNS):
                total[i] += int(row[column])

ncr_values = sorted(totals)
all_sites = sum(total[0] for total in totals.values())
all_obs = sum(total[1] for total in totals.values())
cutoffs = [round(i * STEP, 10) for i in range(int(round(1 / STEP)) + 1)]
if cutoffs[-1] < 1:
    cutoffs.append(1.0)

with open("~{prefix}.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow([
        'cutoff', 'n_sites', 'recall_sites', 'n_trio_carrier_obs', 'recall_trio_carrier_obs',
        'n_transmitted', 'n_denovo', 'precision',
    ])

    # Cutoffs and NCR values are both ascending, so one pass accumulates every row
    cumulative = [0, 0, 0, 0]
    next_value = 0
    for cutoff in cutoffs:
        while next_value < len(ncr_values) and ncr_values[next_value] <= cutoff + EPSILON:
            for i, count in enumerate(totals[ncr_values[next_value]]):
                cumulative[i] += count
            next_value += 1
        n_sites, n_obs, n_transmitted, n_denovo = cumulative
        writer.writerow([
            f'{cutoff:g}', n_sites, ratio(n_sites, all_sites), n_obs, ratio(n_obs, all_obs),
            n_transmitted, n_denovo, ratio(n_transmitted, n_obs),
        ])

print(f'Excluded {missing_ncr_sites} sites without INFO/NCR')
CODE
    >>>

    output {
        File ncr_cutoffs_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(count_tsvs, "GB")) + 10,
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
