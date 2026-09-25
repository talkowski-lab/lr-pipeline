version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FilterLowCallSites {
    meta {
        description: [
            "This utility removes sites that carry no informative genotypes and flags the sites that are largely uncalled. When `remove_no_calls` is set, every variant whose samples carry no alternate allele - an allele count of zero - is removed. When `min_ncr_filter` is set to a non-negative value, every surviving variant whose no-call rate reaches that value is given the `HIGH_NCR` FILTER value, whose header line the workflow adds; those variants are flagged rather than removed.",
            "A genotype counts as a carrier when any of its alleles is alternate, so a partially called genotype such as './1' keeps the site from being an allele count of zero. Carriers are counted from the genotypes rather than read from `INFO/AC`, so an allele count left stale by an upstream step does not affect which sites are removed.",
            "The no-call rate comes from `INFO/NCR`, which is a single value at biallelic and multiallelic sites alike and is used as it stands. A variant that carries no `INFO/NCR` falls back to counting the proportion of alleles without a call, which is what that field holds, and every such variant ID is printed to the task log. It optionally shards by record count and outputs the filtered VCF."
        ]
    }

    parameter_meta {
        vcf: "Cohort VCF to filter."
        vcf_idx: "Index for the cohort VCF."
        remove_no_calls: "Whether to remove sites whose samples carry no alternate allele."
        min_ncr_filter: "No-call rate at or above which a site is given the `HIGH_NCR` FILTER value. A negative value leaves sites unflagged."
        records_per_shard: "Number of variants per shard. When set, variants are processed in parallel shards and concatenated."
        filtered_vcf: "VCF with the sites carrying no alternate allele removed and the high no-call rate sites flagged."
        filtered_vcf_idx: "Index for `filtered_vcf`."
    }

    input {
        File vcf
        File vcf_idx
        String prefix

        Boolean remove_no_calls
        Float min_ncr_filter = -1

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat_vcfs
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
        call FilterLowCallSitesShard {
            input:
                vcf = vcfs_to_process[i],
                vcf_idx = vcf_idxs_to_process[i],
                remove_no_calls = remove_no_calls,
                min_ncr_filter = min_ncr_filter,
                prefix = "~{prefix}.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter
        }
    }

    if (defined(records_per_shard)) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = FilterLowCallSitesShard.filtered_vcf,
                vcf_idxs = FilterLowCallSitesShard.filtered_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.low_call_sites_filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }
    }

    output {
        File filtered_vcf = select_first([ConcatVcfs.concat_vcf, FilterLowCallSitesShard.filtered_vcf[0]])
        File filtered_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, FilterLowCallSitesShard.filtered_vcf_idx[0]])
    }
}

task FilterLowCallSitesShard {
    input {
        File vcf
        File vcf_idx
        Boolean remove_no_calls
        Float min_ncr_filter
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Drop sites with no alternate allele carrier and flag the rest against the no-call rate threshold
        python3 <<'PYCODE'
import pysam

VCF = "~{vcf}"
REMOVE_NO_CALLS = ~{if remove_no_calls then "True" else "False"}
MIN_NCR = ~{min_ncr_filter}
OUTPUT_VCF = "~{prefix}.vcf.gz"
FILTER_NAME = "HIGH_NCR"
FILTER_DESCRIPTION = f"Variant whose no-call rate is at or above {MIN_NCR}."

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
    if MIN_NCR >= 0 and FILTER_NAME not in header.filters:
        header.add_line(f'##FILTER=<ID={FILTER_NAME},Description="{FILTER_DESCRIPTION}">')

    with pysam.VariantFile(OUTPUT_VCF, "wz", header=header) as destination:
        for record in source:
            if REMOVE_NO_CALLS and not has_alt_call(record):
                continue
            record.translate(header)
            if MIN_NCR >= 0:
                no_call_rate = record.info.get("NCR")
                if no_call_rate is None:
                    counted_variant_ids.append(record.id or f"{record.chrom}:{record.pos}")
                    no_call_rate = count_no_call_rate(record)
                if no_call_rate >= MIN_NCR and FILTER_NAME not in record.filter:
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
