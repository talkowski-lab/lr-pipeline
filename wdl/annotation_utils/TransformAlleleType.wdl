version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TransformAlleleType {
    input {
        File vcf
        File vcf_idx
        String prefix

        Int? records_per_shard

        Int dup_breakpoint_window = 10
        Float dup_size_similarity = 0.9
        Int min_dup_size = 50

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_transform
        RuntimeAttr? runtime_attr_concat_shards
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
        call TransformAlleleTypes {
            input:
                vcf = vcfs_to_process[i],
                vcf_idx = vcf_idxs_to_process[i],
                dup_breakpoint_window = dup_breakpoint_window,
                dup_size_similarity = dup_size_similarity,
                min_dup_size = min_dup_size,
                prefix = "~{prefix}.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_transform
        }
    }

    if (defined(records_per_shard)) {
        call Helpers.ConcatVcfs as ConcatShards {
            input:
                vcfs = TransformAlleleTypes.transformed_vcf,
                vcf_idxs = TransformAlleleTypes.transformed_vcf_idx,
                allow_overlaps = false,
                naive = false,
                prefix = "~{prefix}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_shards
        }
    }

    output {
        File transformed_vcf = select_first([ConcatShards.concat_vcf, TransformAlleleTypes.transformed_vcf[0]])
        File transformed_vcf_idx = select_first([ConcatShards.concat_vcf_idx, TransformAlleleTypes.transformed_vcf_idx[0]])
    }
}

task TransformAlleleTypes {
    input {
        File vcf
        File vcf_idx
        Int dup_breakpoint_window
        Float dup_size_similarity
        Int min_dup_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import re
import pysam

dup_breakpoint_window = ~{dup_breakpoint_window}
dup_size_similarity = ~{dup_size_similarity}
min_dup_size = ~{min_dup_size}

INS_SUBTYPES = {"complex_dup", "dup_interspersed", "inv_dup", "alu_ins", "line_ins", "sva_ins", "numt"}
DEL_SUBTYPES = {"alu_del", "line_del", "sva_del"}


def parse_origin(origin):
    m = re.match(r'^(.+):(\d+)-(\d+)_?[+-]$', origin)
    if m is None:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))


def passes_criteria(pos, dup_length, origin_start, origin_end):
    origin_length = origin_end - origin_start
    larger = max(dup_length, origin_length)
    if larger == 0 or min(dup_length, origin_length) / larger < dup_size_similarity:
        return False
    if origin_start <= pos <= origin_end:
        return True
    return abs(origin_start - pos) <= dup_breakpoint_window or abs(origin_end - pos) <= dup_breakpoint_window


def is_tandem(record):
    allele_length = record.info.get("allele_length")
    if isinstance(allele_length, (list, tuple)):
        allele_length = allele_length[0]
    if allele_length is None or abs(int(allele_length)) < min_dup_size:
        return False

    origin = record.info.get("ORIGIN")
    if origin is None:
        return False
    origin = origin if isinstance(origin, str) else ",".join(origin)
    if "," in origin:
        return False

    parsed = parse_origin(origin)
    if parsed is None:
        return False

    origin_chrom, origin_start, origin_end = parsed
    dup_length = abs(int(allele_length))

    return passes_criteria(record.pos, dup_length, origin_start, origin_end)


vcf_in = pysam.VariantFile("~{vcf}")
header = vcf_in.header.copy()
header.info.add("allele_subtype", 1, "String", "Original allele_type subtype for reclassified variants: tandem_dup, dup_interspersed, complex_dup, inv_dup, alu_ins, line_ins, sva_ins, numt, alu_del, line_del or sva_del.")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)

for record in vcf_in:
    record.translate(vcf_out.header)

    allele_type = record.info.get("allele_type")
    if isinstance(allele_type, (list, tuple)):
        allele_type = allele_type[0]

    if allele_type == "dup":
        record.info["allele_subtype"] = "tandem_dup"
        if not is_tandem(record):
            record.info["allele_type"] = "ins"
    elif allele_type in INS_SUBTYPES:
        record.info["allele_subtype"] = allele_type
        record.info["allele_type"] = "ins"
    elif allele_type in DEL_SUBTYPES:
        record.info["allele_subtype"] = allele_type
        record.info["allele_type"] = "del"

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File transformed_vcf = "~{prefix}.vcf.gz"
        File transformed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(size(vcf, "GB")) + 4,
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
