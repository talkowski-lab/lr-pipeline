version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FilterNCRCutoffs {
    meta {
        description: [
            "This utility removes sites that carry no informative genotypes and flags the sites that are largely uncalled. When `remove_no_calls` is set, every variant whose samples carry no alternate allele - an allele count of zero - is removed. When `ncr_cutoffs_tsv` is given, every surviving variant is matched to a row by its lowercased `INFO/allele_type`, its absolute `INFO/allele_length` within the row's inclusive length bounds and its `INFO/REGION`, where a row's `ALL` region matches any region. A matched variant whose no-call rate is above the row's cutoff is given the `HIGH_NCR` FILTER value, whose header line the workflow adds; those variants are flagged rather than removed. A cutoff keeps the sites at or below it, as an `EvaluateNCRCutoffs` cutoff does, so a cutoff read off its curves keeps the same sites here. Variants matching no row are left untouched, so the table need not cover every stratum. When `ncr_filter_pass_only` is set, only variants whose FILTER is PASS are matched, and `ncr_cutoffs_tsv` must be given; otherwise variants of every FILTER are matched.",
            "A genotype counts as a carrier when any of its alleles is alternate, so a partially called genotype such as './1' keeps the site from being an allele count of zero. Carriers are counted from the genotypes rather than read from `INFO/AC`, so an allele count left stale by an upstream step does not affect which sites are removed.",
            "Variants are classified the way `EvaluateNCRCutoffs` strata are, so cutoffs chosen from its output apply to the same variants. Rows whose strata overlap are rejected when the table is read, and a variant checked against the table without `INFO/REGION` stops the task.",
            "The no-call rate comes from `INFO/NCR`, which is a single value at biallelic and multiallelic sites alike and is used as it stands. A matched variant that carries no `INFO/NCR` falls back to counting the proportion of alleles without a call, which is what that field holds, and every such variant ID is printed to the task log. It optionally shards by record count and outputs the filtered VCF, whose header, when `clean_vcf_header` is set, is cleaned last of all to keep only the primary contigs and the INFO, FILTER, FORMAT and ALT definitions that records use, dropping tool provenance lines such as command lines and program versions."
        ]
    }

    parameter_meta {
        vcf: "Cohort VCF to filter."
        vcf_idx: "Index for the cohort VCF."
        remove_no_calls: "Whether to remove sites whose samples carry no alternate allele."
        ncr_cutoffs_tsv: "Headed TSV of no-call rate cutoffs with columns `allele_type`, `min_length`, `max_length`, `region` and `cutoff`. Lengths are inclusive bounds on the absolute `INFO/allele_length`, `region` is `US`, `RM`, `SD`, `SR` or `ALL` for any region, and a matched site is given the `HIGH_NCR` FILTER value when its no-call rate is above `cutoff`. Rows may cover only some strata but may not overlap. When absent, no site is flagged."
        ncr_filter_pass_only: "Whether to check only sites whose FILTER is PASS against `ncr_cutoffs_tsv`, which must then be given. When false, sites of every FILTER are checked."
        records_per_shard: "Number of variants per shard. When set, variants are processed in parallel shards and concatenated."
        clean_vcf_header: "Whether to clean the header of the filtered VCF, last of all, keeping only the primary contigs and the definitions that records use."
        filtered_vcf: "VCF with the sites carrying no alternate allele removed and the high no-call rate sites flagged."
        filtered_vcf_idx: "Index for `filtered_vcf`."
    }

    input {
        File vcf
        File vcf_idx
        String prefix

        Boolean remove_no_calls
        File? ncr_cutoffs_tsv
        Boolean ncr_filter_pass_only

        Int? records_per_shard
        Boolean clean_vcf_header

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_clean_header
    }

    if (defined(records_per_shard)) {
        call Helpers.ShardVcfByRecords {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                records_per_shard = select_first([records_per_shard]),
                prefix = "~{prefix}.sharded",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }
    }

    Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [vcf]])
    Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [vcf_idx]])

    scatter (i in range(length(vcfs_to_process))) {
        call FilterNCRCutoffsShard {
            input:
                vcf = vcfs_to_process[i],
                vcf_idx = vcf_idxs_to_process[i],
                remove_no_calls = remove_no_calls,
                ncr_cutoffs_tsv = ncr_cutoffs_tsv,
                ncr_filter_pass_only = ncr_filter_pass_only,
                prefix = "~{prefix}.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter
        }
    }

    if (defined(records_per_shard)) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = FilterNCRCutoffsShard.filtered_vcf,
                vcf_idxs = FilterNCRCutoffsShard.filtered_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.ncr_filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }
    }

    File unclean_vcf = select_first([ConcatVcfs.concat_vcf, FilterNCRCutoffsShard.filtered_vcf[0]])
    File unclean_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, FilterNCRCutoffsShard.filtered_vcf_idx[0]])

    if (clean_vcf_header) {
        call Helpers.CleanVcfHeader {
            input:
                vcf = unclean_vcf,
                vcf_idx = unclean_vcf_idx,
                prefix = "~{prefix}.cleaned",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_clean_header
        }
    }

    output {
        File filtered_vcf = select_first([CleanVcfHeader.cleaned_vcf, unclean_vcf])
        File filtered_vcf_idx = select_first([CleanVcfHeader.cleaned_vcf_idx, unclean_vcf_idx])
    }
}

task FilterNCRCutoffsShard {
    input {
        File vcf
        File vcf_idx
        Boolean remove_no_calls
        File? ncr_cutoffs_tsv
        Boolean ncr_filter_pass_only
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Drop sites with no alternate allele carrier and flag the rest against their stratum's no-call rate cutoff
        python3 <<'PYCODE'
import csv
import sys

import pysam

VCF = "~{vcf}"
REMOVE_NO_CALLS = ~{if remove_no_calls then "True" else "False"}
CUTOFFS_TSV = "~{default='' ncr_cutoffs_tsv}"
PASS_ONLY = ~{if ncr_filter_pass_only then "True" else "False"}
OUTPUT_VCF = "~{prefix}.vcf.gz"
FILTER_NAME = "HIGH_NCR"
FILTER_DESCRIPTION = "Variant whose no-call rate is above the cutoff for its allele type, length and region."
REGIONS = {"US", "RM", "SD", "SR", "ALL"}
# EvaluateNCRCutoffs keeps sites with NCR <= cutoff + EPSILON, so the same sites are kept here
EPSILON = 1e-9

def read_cutoffs(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"allele_type", "min_length", "max_length", "region", "cutoff"}
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            raise ValueError("ncr_cutoffs_tsv must contain allele_type, min_length, max_length, region and cutoff columns")
        cutoffs = []
        for row in reader:
            region = row["region"].upper()
            if region not in REGIONS:
                raise ValueError(f"ncr_cutoffs_tsv region must be one of {sorted(REGIONS)}, not {row['region']}")
            cutoffs.append((row["allele_type"].lower(), int(row["min_length"]), int(row["max_length"]), region,
                            float(row["cutoff"])))
    # A variant matching two rows would have no single cutoff
    for i, (type_a, min_a, max_a, region_a, _) in enumerate(cutoffs):
        for type_b, min_b, max_b, region_b, _ in cutoffs[i + 1:]:
            same_region = region_a == region_b or "ALL" in (region_a, region_b)
            if type_a == type_b and min_a <= max_b and min_b <= max_a and same_region:
                raise ValueError(f"ncr_cutoffs_tsv rows overlap: {type_a} {min_a}-{max_a} {region_a} "
                                 f"and {type_b} {min_b}-{max_b} {region_b}")
    return cutoffs

# Classify as EvaluateNCRCutoffs does: lowercased type, absolute length, and the record's region
def find_cutoff(record):
    allele_type = record.info.get("allele_type")
    allele_length = record.info.get("allele_length")
    if allele_type is None or allele_length is None:
        return None
    region = record.info.get("REGION")
    if region is None:
        sys.exit(f"INFO/REGION is missing at {record.chrom}:{record.pos}")
    length = abs(int(allele_length))
    for row_type, min_length, max_length, row_region, cutoff in CUTOFFS:
        if row_type == allele_type.lower() and min_length <= length <= max_length and row_region in ("ALL", region):
            return cutoff
    return None

if PASS_ONLY and not CUTOFFS_TSV:
    sys.exit("ncr_filter_pass_only is set, but no ncr_cutoffs_tsv was given to apply to PASS sites")
CUTOFFS = read_cutoffs(CUTOFFS_TSV) if CUTOFFS_TSV else []

def has_alt_call(record):
    for sample in record.samples.values():
        gt = sample.get("GT")
        if gt is None:
            continue
        if any(allele is not None and allele > 0 for allele in gt):
            return True
    return False

# Count the no-call rate the way INFO/NCR defines it, as the proportion of alleles without a call
def count_no_call_rate(record):
    alleles = 0
    no_calls = 0
    for sample in record.samples.values():
        gt = sample.get("GT")
        if gt is None:
            continue
        alleles += len(gt)
        no_calls += sum(1 for allele in gt if allele is None)
    if alleles == 0:
        return 1.0
    return no_calls / alleles

counted_variant_ids = []

with pysam.VariantFile(VCF) as source:
    header = source.header.copy()
    if CUTOFFS and FILTER_NAME not in header.filters:
        header.add_line(f'##FILTER=<ID={FILTER_NAME},Description="{FILTER_DESCRIPTION}">')

    with pysam.VariantFile(OUTPUT_VCF, "wz", header=header) as destination:
        for record in source:
            if REMOVE_NO_CALLS and not has_alt_call(record):
                continue
            record.translate(header)
            # Checking FILTER first spares non-PASS records the REGION check when only PASS is filtered
            cutoff = find_cutoff(record) if CUTOFFS and (not PASS_ONLY or "PASS" in record.filter) else None
            if cutoff is not None:
                no_call_rate = record.info.get("NCR")
                if no_call_rate is None:
                    counted_variant_ids.append(record.id or f"{record.chrom}:{record.pos}")
                    no_call_rate = count_no_call_rate(record)
                # NCR is a 32-bit float, so round to the 6 digits bcftools prints for EvaluateNCRCutoffs to keep sites at the cutoff
                if float(f"{no_call_rate:.6g}") > cutoff + EPSILON and FILTER_NAME not in record.filter:
                    if "PASS" in record.filter:
                        record.filter.clear()
                    record.filter.add(FILTER_NAME)
            destination.write(record)

if counted_variant_ids:
    print(f"Counted the no-call rate for {len(counted_variant_ids)} variants carrying no INFO/NCR:")
    for variant_id in counted_variant_ids:
        print(variant_id)
PYCODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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
