version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow DropGenotypes {
    input {
        File vcf
        File vcf_idx
        String prefix

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_strip
        RuntimeAttr? runtime_attr_concat
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
        call Helpers.StripGenotypes as StripShardGenotypes {
            input:
                vcf = vcfs_to_process[i],
                vcf_idx = vcf_idxs_to_process[i],
                prefix = "~{prefix}.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_strip
        }
    }

    if (defined(records_per_shard)) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = StripShardGenotypes.stripped_vcf,
                vcf_idxs = StripShardGenotypes.stripped_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File dropped_vcf = select_first([ConcatVcfs.concat_vcf, StripShardGenotypes.stripped_vcf[0]])
        File dropped_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, StripShardGenotypes.stripped_vcf_idx[0]])
    }
}
