version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FilterLowCoverageRegions {
    input {
        File vcf
        File vcf_idx
        File low_coverage_regions_bed
        Array[String] contigs
        String prefix

        Float min_region_coverage_cutoff
        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_bed
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }
        }

        call Helpers.SubsetBedToContig {
            input:
                bed = low_coverage_regions_bed,
                contig = contig,
                prefix = "~{prefix}.~{contig}.low_coverage_regions",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_bed
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (s in range(length(shard_vcfs))) {
            call AddLowCoverageRegionFilter {
                input:
                    vcf = shard_vcfs[s],
                    vcf_idx = shard_vcf_idxs[s],
                    low_coverage_regions_bed = SubsetBedToContig.subset_bed,
                    min_region_coverage_cutoff = min_region_coverage_cutoff,
                    prefix = "~{prefix}.~{contig}.shard_~{s}.low_coverage_filtered",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_filter
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatShards {
                input:
                    vcfs = AddLowCoverageRegionFilter.filtered_vcf,
                    vcf_idxs = AddLowCoverageRegionFilter.filtered_vcf_idx,
                    allow_overlaps = false,
                    naive = true,
                    prefix = "~{prefix}.~{contig}.low_coverage_filtered",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File contig_filtered_vcf = select_first([ConcatShards.concat_vcf, AddLowCoverageRegionFilter.filtered_vcf[0]])
        File contig_filtered_vcf_idx = select_first([ConcatShards.concat_vcf_idx, AddLowCoverageRegionFilter.filtered_vcf_idx[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = contig_filtered_vcf,
                vcf_idxs = contig_filtered_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.low_coverage_filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File low_coverage_region_filtered_vcf = select_first([ConcatVcfs.concat_vcf, contig_filtered_vcf[0]])
        File low_coverage_region_filtered_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, contig_filtered_vcf_idx[0]])
    }
}

task AddLowCoverageRegionFilter {
    input {
        File vcf
        File vcf_idx
        File low_coverage_regions_bed
        Float min_region_coverage_cutoff
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view -H "~{vcf}" \
        | awk 'BEGIN { OFS="\t" } {
            start = $2 - 1
            print $1, start, start + length($4), NR
        }' > records.bed

        sort -k1,1 -k2,2n "~{low_coverage_regions_bed}" > sorted_low_coverage_regions.bed
        bedtools merge -i sorted_low_coverage_regions.bed > merged_low_coverage_regions.bed

        bedtools coverage \
            -a records.bed \
            -b merged_low_coverage_regions.bed \
        | awk -v cutoff="~{min_region_coverage_cutoff}" '$6 / $7 >= cutoff { print $4 }' > flagged_record_numbers.txt

        python3 <<'PYTHON'
import pysam

filter_name = "LOW_COVERAGE_REGION"
description = "Variant region determined to have low coverage relative to entire genome."

with open("flagged_record_numbers.txt") as flagged_file:
    flagged_record_numbers = {int(line) for line in flagged_file if line.strip()}

with pysam.VariantFile("~{vcf}") as source:
    header = source.header.copy()
    if filter_name not in header.filters:
        header.add_line(f'##FILTER=<ID={filter_name},Description="{description}">')

    with pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header) as destination:
        for record_number, record in enumerate(source, start=1):
            record.translate(header)
            if record_number in flagged_record_numbers and filter_name not in record.filter:
                if "PASS" in record.filter:
                    record.filter.clear()
                record.filter.add(filter_name)
            destination.write(record)
PYTHON

        tabix -f -p vcf "~{prefix}.vcf.gz"
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size([vcf, vcf_idx, low_coverage_regions_bed], "GB")) + 10,
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
