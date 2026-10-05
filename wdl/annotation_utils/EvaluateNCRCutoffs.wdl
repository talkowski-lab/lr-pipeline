version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow EvaluateNCRCutoffs {
    meta {
        description: [
            "This utility tabulates de novo rate and recall of a callset across a sweep of `INFO/NCR` (no-call rate) cutoffs, so a cutoff can be chosen from the resulting curve. A site passes a cutoff when it passes `args_string_vcf` and its `NCR` is at or below the cutoff. Every row is stratified by `INFO/allele_type` and by size bin of `abs(INFO/allele_length)` over `length_bins`, with `all` marginals for each, and sites missing either field are binned as `missing`. Recall is the fraction of sites in the stratum retained relative to every site on the contig, before `args_string_vcf` and including sites without `INFO/NCR`, which fail every cutoff. The filter is applied as a soft filter so this denominator is counted in the same pass.",
            "The de novo rate is estimated from trios whose child and both parents are VCF samples in the PED. Each passing site where a trio child has any called non-reference allele is one observation, so phasing, zygosity and partially missing child genotypes such as `./1` do not matter. An observation is transmitted when any called allele of a relevant parent is non-reference, de novo when every relevant parent has at least its expected number of called alleles and all are reference, and skipped otherwise. The de novo rate is de novo over transmitted plus de novo observations, and skipped observations are counted separately by reason.",
            "Relevant parents follow inheritance on GRCh38. On autosomes, the pseudoautosomal regions and any other contig, both parents are relevant with two expected alleles each. On non-PAR chrX, sons inherit from the mother only (two alleles) and daughters from both parents (father one allele, mother two). On non-PAR chrY, sons inherit from the father only (one allele) and daughters have no parental source, so their calls are skipped as `no_parent`. On chrM the mother alone is relevant with one allele. On non-PAR chrX and chrY, children whose PED sex is neither `1` nor `2` are skipped as `unknown_sex`. Records are placed in or out of a pseudoautosomal region by `POS`."
        ]
    }

    parameter_meta {
        vcf: "VCF carrying `INFO/NCR`, `INFO/allele_type` and `INFO/allele_length` to evaluate."
        vcf_idx: "Index for `vcf`."
        contig: "Contig to evaluate."
        records_per_shard: "Number of variants to keep within a single shard during evaluation."
        ped: "Six-column PED used to identify trios and the sex of each child."
        args_string_vcf: "`bcftools filter` include expression marking the sites eligible to pass a cutoff; sites failing it still count toward the recall denominator."
        length_bins: "Size-bin edges for stratifying sites by `abs(INFO/allele_length)`."
        cutoff_step: "Spacing of the `NCR` cutoffs evaluated from `0` to `1` inclusive."
        ncr_cutoffs_tsv: "One row per `allele_type`, `size_bin` and `cutoff`, with `n_sites_unfiltered` in the stratum, `n_sites` passing the cutoff, `recall`, `n_transmitted`, `n_denovo`, `denovo_rate` (empty when there are no classified observations) and the skipped observation counts `n_skipped_parent_missing`, `n_skipped_unknown_sex` and `n_skipped_no_parent`."
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
                trio_sample_ids_file = FindTrios.trio_sample_ids_file,
                ped = ped,
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
        File trio_sample_ids_file
        File ped
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

        # Soft filter on all samples before trio subsetting, so failing sites still reach the recall denominator
        bcftools filter \
            -r ~{contig} \
            ~{if defined(include_args) then "-i '~{include_args}' -s EVAL_EXCLUDED -m+" else ""} \
            -Ou \
            ~{vcf} \
            | bcftools query \
                -H \
                -S ~{trio_sample_ids_file} \
                -f '%POS\t%INFO/NCR\t%INFO/allele_type\t%INFO/allele_length\t%FILTER[\t%GT]\n' \
                > genotypes.tsv

        python3 <<'CODE'
import csv
import re
from collections import Counter, defaultdict

ALLELE_SPLIT = re.compile(r'[/|]')
SEX_BY_CODE = {'1': 'male', '2': 'female'}
# GRCh38 pseudoautosomal regions, 1-based inclusive
PAR_REGIONS = {'X': [(10001, 2781479), (155701383, 156030895)], 'Y': [(10001, 2781479), (56887903, 57217415)]}
LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f'{start}-{end - 1}' for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f'{LENGTH_BINS[-1]}+']
CONTIG = re.sub(r'^chr', '', "~{contig}", flags=re.IGNORECASE).upper()
OUTCOMES = ['transmitted', 'denovo', 'skipped_parent_missing', 'skipped_unknown_sex', 'skipped_no_parent']

def contig_class(pos):
    if CONTIG in ('M', 'MT'):
        return 'chrM'
    if CONTIG not in PAR_REGIONS or any(start <= pos <= end for start, end in PAR_REGIONS[CONTIG]):
        return 'autosome'
    return 'chr' + CONTIG

def expected_parent_alleles(cls, sex):
    """Called alleles each relevant parent needs before a reference genotype rules out transmission."""
    if cls == 'autosome':
        return {'father': 2, 'mother': 2}
    if cls == 'chrM':
        return {'mother': 1}
    if sex is None:
        return None
    if cls == 'chrX':
        return {'mother': 2} if sex == 'male' else {'father': 1, 'mother': 2}
    return {'father': 1} if sex == 'male' else {}

def called_alleles(gt):
    return [allele for allele in ALLELE_SPLIT.split(gt) if allele != '.']

def classify(child, parents, expected):
    if not any(allele != '0' for allele in called_alleles(child)):
        return None
    if expected is None:
        return 'skipped_unknown_sex'
    if not expected:
        return 'skipped_no_parent'
    called = {parent: called_alleles(parents[parent]) for parent in expected}
    if any(allele != '0' for alleles in called.values() for allele in alleles):
        return 'transmitted'
    if all(len(called[parent]) >= n for parent, n in expected.items()):
        return 'denovo'
    return 'skipped_parent_missing'

def size_bin(length):
    if length == '.':
        return 'missing'
    size = abs(int(length))
    for index, label in enumerate(SIZE_LABELS):
        if index + 1 == len(LENGTH_BINS) or size < LENGTH_BINS[index + 1]:
            return label

with open("~{ped}") as handle:
    sex_by_sample = {}
    for line in handle:
        if line.startswith('#') or not line.strip():
            continue
        fields = line.rstrip('\n').split('\t')
        sex_by_sample[fields[1]] = SEX_BY_CODE.get(fields[4]) if len(fields) > 4 else None

with open("~{trio_definitions}") as handle:
    trios = [line.rstrip('\n').split('\t')[:3] for line in handle if line.strip()]

# Keyed by allele type, size bin and NCR as printed, so the summary task can apply any cutoff step exactly
COUNT_HEADER = ['n_sites_unfiltered', 'n_sites'] + [f'n_{outcome}' for outcome in OUTCOMES]
counts = defaultdict(lambda: [0] * len(COUNT_HEADER))
outcome_totals = Counter()
expected_cache = {}
with open('genotypes.tsv') as handle:
    header = next(handle).lstrip('#').rstrip('\n').split('\t')
    columns = {re.sub(r'^\s*\[\d+\]|:GT$', '', name): i for i, name in enumerate(header)}
    trio_columns = [(columns[child], columns[father], columns[mother], sex_by_sample.get(child))
                    for child, father, mother in trios]
    for line in handle:
        fields = line.rstrip('\n').split('\t')
        pos, ncr, allele_type, allele_length, filters = fields[:5]
        allele_type = allele_type.lower() if allele_type != '.' else 'missing'
        row = counts[(allele_type, size_bin(allele_length), ncr)]
        row[0] += 1
        if 'EVAL_EXCLUDED' in filters.split(';'):
            continue
        row[1] += 1
        cls = contig_class(int(pos))
        for child_col, father_col, mother_col, sex in trio_columns:
            key = (cls, sex)
            if key not in expected_cache:
                expected_cache[key] = expected_parent_alleles(cls, sex)
            parents = {'father': fields[father_col], 'mother': fields[mother_col]}
            outcome = classify(fields[child_col], parents, expected_cache[key])
            if outcome is not None:
                row[2 + OUTCOMES.index(outcome)] += 1
                outcome_totals[outcome] += 1

with open("~{prefix}.ncr_counts.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['allele_type', 'size_bin', 'ncr'] + COUNT_HEADER)
    for key, row in counts.items():
        writer.writerow([*key, *row])

print(f'Counted {sum(row[0] for row in counts.values())} sites across {len(trios)} trios: '
      + ', '.join(f'{outcome_totals[o]} {o}' for o in OUTCOMES))
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

SKIP_COLUMNS = ['n_skipped_parent_missing', 'n_skipped_unknown_sex', 'n_skipped_no_parent']
COUNT_COLUMNS = ['n_sites', 'n_transmitted', 'n_denovo'] + SKIP_COLUMNS
LENGTH_BINS = [~{sep=", " length_bins}]
SIZE_LABELS = [f'{start}-{end - 1}' for start, end in zip(LENGTH_BINS, LENGTH_BINS[1:])] + [f'{LENGTH_BINS[-1]}+']
SIZE_ORDER = {label: i for i, label in enumerate(['all'] + SIZE_LABELS + ['missing'])}
STEP = ~{cutoff_step}
EPSILON = 1e-9

def ratio(numerator, denominator):
    return f'{numerator / denominator:.6f}' if denominator else ''

# Each shard row feeds its own stratum and the three marginals over allele type and size bin
totals = defaultdict(lambda: defaultdict(lambda: [0] * len(COUNT_COLUMNS)))
unfiltered = defaultdict(int)
missing_ncr_sites = 0
for path in "~{sep=',' count_tsvs}".split(','):
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            allele_type, size, ncr = row['allele_type'], row['size_bin'], row['ncr']
            strata = [(allele_type, size), (allele_type, 'all'), ('all', size), ('all', 'all')]
            for stratum in strata:
                unfiltered[stratum] += int(row['n_sites_unfiltered'])
            if ncr == '.':
                missing_ncr_sites += int(row['n_sites_unfiltered'])
                continue
            for stratum in strata:
                total = totals[stratum][float(ncr)]
                for i, column in enumerate(COUNT_COLUMNS):
                    total[i] += int(row[column])

cutoffs = [round(i * STEP, 10) for i in range(int(round(1 / STEP)) + 1)]
if cutoffs[-1] < 1:
    cutoffs.append(1.0)

def stratum_order(stratum):
    allele_type, size = stratum
    return (allele_type != 'all', allele_type == 'missing', allele_type, SIZE_ORDER.get(size, len(SIZE_ORDER)))

with open("~{prefix}.tsv", 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow([
        'allele_type', 'size_bin', 'cutoff', 'n_sites_unfiltered', 'n_sites', 'recall',
        'n_transmitted', 'n_denovo', 'denovo_rate', *SKIP_COLUMNS,
    ])
    for stratum in sorted(unfiltered, key=stratum_order):
        by_ncr = totals[stratum]
        ncr_values = sorted(by_ncr)
        n_unfiltered = unfiltered[stratum]

        # Cutoffs and NCR values are both ascending, so one pass accumulates every row
        cumulative = [0] * len(COUNT_COLUMNS)
        next_value = 0
        for cutoff in cutoffs:
            while next_value < len(ncr_values) and ncr_values[next_value] <= cutoff + EPSILON:
                for i, count in enumerate(by_ncr[ncr_values[next_value]]):
                    cumulative[i] += count
                next_value += 1
            n_sites, n_transmitted, n_denovo, *skipped = cumulative
            writer.writerow([
                *stratum, f'{cutoff:g}', n_unfiltered, n_sites, ratio(n_sites, n_unfiltered),
                n_transmitted, n_denovo, ratio(n_denovo, n_transmitted + n_denovo), *skipped,
            ])

print(f'{missing_ncr_sites} sites without INFO/NCR count toward recall denominators but pass no cutoff')
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
