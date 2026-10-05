version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow EvaluateNCRCutoffs {
    meta {
        description: [
            "This utility counts the precision and recall of a callset across a sweep of `INFO/NCR` (no-call rate) cutoffs, so a cutoff can be chosen from the resulting curve. Only sites passing `args_string_vcf` are considered, and a site passes a cutoff when its `NCR` is at or below it. Both metrics are counted in non-reference calls: each site where a sample has any called non-reference allele is one call, tagged with the site's `NCR`, so phasing and zygosity are ignored. Every row is stratified by `INFO/allele_type` and by size bin of `abs(INFO/allele_length)` over `length_bins`, with `all` marginals for each, and sites missing either field are binned as `missing`. Calls at sites without `INFO/NCR` are kept out of the sweep and reported in a separate row per stratum whose `cutoff` is `missing`.",
            "Recall is the calls of every VCF sample at sites passing the cutoff over all their calls at sites with `INFO/NCR`. Precision uses only the calls of trio probands, from trios whose child and both parents are VCF samples in the PED. A proband call is inherited when a relevant parent has any called non-reference allele at the site, skipped when no relevant parent carries it but a relevant parent has a missing allele, and uninherited otherwise. Both parents are relevant everywhere except chrY, where only the father is. Only the genotypes are consulted, so haploid calls such as a father's `1` on chrX need no special handling. Precision is inherited over inherited plus uninherited proband calls."
        ]
    }

    parameter_meta {
        vcf: "VCF carrying `INFO/NCR`, `INFO/allele_type` and `INFO/allele_length` to evaluate."
        vcf_idx: "Index for `vcf`."
        contig: "Contig to evaluate."
        records_per_shard: "Number of variants to keep within a single shard during evaluation."
        ped: "Six-column PED used to identify trios."
        args_string_vcf: "`bcftools view` include expression applied to the VCF before evaluation; sites failing it enter no count."
        length_bins: "Size-bin edges for stratifying sites by `abs(INFO/allele_length)`."
        cutoff_step: "Spacing of the `NCR` cutoffs evaluated from `0` to `1` inclusive."
        ncr_cutoffs_tsv: "One row per `allele_type`, `size_bin` and `cutoff`, plus a `missing` cutoff row per stratum for sites without `INFO/NCR`, with proband call counts `prec_numerator` (inherited), `prec_denominator` (inherited and uninherited) and `prec_skipped` (skipped), and all-sample call counts `rec_numerator` (calls in the row) and `rec_denominator` (calls at sites with `INFO/NCR`, or the row's own calls in the `missing` row)."
    }

    input {
        File vcf
        File vcf_idx
        String contig
        Int? records_per_shard
        String prefix

        File ped
        String? args_string_vcf
        Array[Int] length_bins = [0, 1, 50, 500]
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
                include_args = args_string_vcf,
                length_bins = length_bins,
                prefix = "~{prefix}.~{contig}.shard_~{shard_i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_count_trio_transmission
        }
    }

    call SummarizeNCRCutoffs {
        input:
            count_tsvs = CountTrioTransmission.ncr_counts_tsv,
            length_bins = length_bins,
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
        String? include_args
        Array[Int] length_bins
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

        bcftools view \
            -r ~{contig} \
            ~{if defined(include_args) then "-i '~{include_args}'" else ""} \
            -Ou \
            ~{vcf} \
            | bcftools query \
                -H \
                -f '%INFO/NCR\t%INFO/allele_type\t%INFO/allele_length[\t%GT]\n' \
                > genotypes.tsv

        python3 <<'CODE'
import csv
import re
from collections import defaultdict

# Any allele index other than 0 contains a digit from 1 to 9, so one search finds a called non-reference allele
NONREF_ALLELE = re.compile(r'[1-9]')
LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f'{start}-{end - 1}' for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f'{LENGTH_BINS[-1]}+']
OUTCOMES = ['inherited', 'uninherited', 'skipped']
# A Y chromosome comes from the father alone, so the mother's genotype there is never consulted
FATHER_ONLY = re.sub(r'^chr', '', "~{contig}", flags=re.IGNORECASE).upper() == 'Y'

def has_nonref(gt):
    return NONREF_ALLELE.search(gt) is not None

def has_missing(gt):
    return '.' in gt

def classify(child, parents):
    if not has_nonref(child):
        return None
    if any(has_nonref(parent) for parent in parents):
        return 'inherited'
    if any(has_missing(parent) for parent in parents):
        return 'skipped'
    return 'uninherited'

def size_bin(length):
    if length == '.':
        return 'missing'
    size = abs(int(length))
    for index, label in enumerate(SIZE_LABELS):
        if index + 1 == len(LENGTH_BINS) or size < LENGTH_BINS[index + 1]:
            return label

with open("~{trio_definitions}") as handle:
    trios = [line.rstrip('\n').split('\t')[:3] for line in handle if line.strip()]

# Keyed by allele type, size bin and NCR as printed; holds one count per proband outcome, then all-sample calls
counts = defaultdict(lambda: [0] * (len(OUTCOMES) + 1))
with open('genotypes.tsv') as handle:
    header = next(handle).lstrip('#').rstrip('\n').split('\t')
    columns = {re.sub(r'^\s*\[\d+\]|:GT$', '', name): i for i, name in enumerate(header)}
    trio_columns = [(columns[child], [columns[father]] if FATHER_ONLY else [columns[father], columns[mother]])
                    for child, father, mother in trios]
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        ncr, allele_type, allele_length = fields[:3]
        allele_type = allele_type.lower() if allele_type != '.' else 'missing'
        row = counts[(allele_type, size_bin(allele_length), ncr)]
        for child_col, parent_cols in trio_columns:
            outcome = classify(fields[child_col], [fields[col] for col in parent_cols])
            if outcome is not None:
                row[OUTCOMES.index(outcome)] += 1
        row[-1] += sum(1 for gt in fields[3:] if has_nonref(gt))

with open("~{prefix}.ncr_counts.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['allele_type', 'size_bin', 'ncr'] + [f'n_{outcome}' for outcome in OUTCOMES] + ['n_calls'])
    for key, row in counts.items():
        writer.writerow([*key, *row])

totals = [sum(row[i] for row in counts.values()) for i in range(len(OUTCOMES) + 1)]
print(f'Counted {len(trios)} trios: ' + ', '.join(f'{n} {outcome}' for n, outcome in zip(totals, OUTCOMES))
      + f'; {totals[-1]} calls across {len(header) - 3} samples')
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
        Array[Int] length_bins
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

COUNT_COLUMNS = ['n_inherited', 'n_uninherited', 'n_skipped', 'n_calls']
LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f'{start}-{end - 1}' for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f'{LENGTH_BINS[-1]}+']
SIZE_ORDER = {label: i for i, label in enumerate(['all'] + SIZE_LABELS + ['missing'])}
STEP = ~{cutoff_step}
EPSILON = 1e-9

def new_counts():
    return [0] * len(COUNT_COLUMNS)

def add(target, counts):
    for i, count in enumerate(counts):
        target[i] += count

def write_row(writer, stratum, cutoff, counts, rec_denominator):
    inherited, uninherited, skipped, calls = counts
    writer.writerow([*stratum, cutoff, inherited, inherited + uninherited, skipped, calls, rec_denominator])

# Each shard row feeds its own stratum and the three marginals over allele type and size bin
by_ncr = defaultdict(lambda: defaultdict(new_counts))
missing_ncr = defaultdict(new_counts)
for path in "~{sep=',' count_tsvs}".split(','):
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            allele_type, size, ncr = row['allele_type'], row['size_bin'], row['ncr']
            counts = [int(row[column]) for column in COUNT_COLUMNS]
            for stratum in [(allele_type, size), (allele_type, 'all'), ('all', size), ('all', 'all')]:
                add(missing_ncr[stratum] if ncr == '.' else by_ncr[stratum][float(ncr)], counts)

cutoffs = [round(i * STEP, 10) for i in range(int(round(1 / STEP)) + 1)]
if cutoffs[-1] < 1:
    cutoffs.append(1.0)

def stratum_order(stratum):
    allele_type, size = stratum
    return (allele_type != 'all', allele_type == 'missing', allele_type, SIZE_ORDER.get(size, len(SIZE_ORDER)))

with open("~{prefix}.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow([
        'allele_type', 'size_bin', 'cutoff',
        'prec_numerator', 'prec_denominator', 'prec_skipped', 'rec_numerator', 'rec_denominator',
    ])
    for stratum in sorted(set(by_ncr) | set(missing_ncr), key=stratum_order):
        totals = by_ncr[stratum]
        ncr_values = sorted(totals)
        rec_denominator = sum(counts[-1] for counts in totals.values())

        # Cutoffs and NCR values are both ascending, so one pass accumulates every row
        cumulative = new_counts()
        next_value = 0
        for cutoff in cutoffs:
            while next_value < len(ncr_values) and ncr_values[next_value] <= cutoff + EPSILON:
                add(cumulative, totals[ncr_values[next_value]])
                next_value += 1
            write_row(writer, stratum, f'{cutoff:g}', cumulative, rec_denominator)

        if stratum in missing_ncr:
            write_row(writer, stratum, 'missing', missing_ncr[stratum], missing_ncr[stratum][-1])

print(f"{missing_ncr[('all', 'all')][-1]} calls at sites without INFO/NCR reported in the missing rows")
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
