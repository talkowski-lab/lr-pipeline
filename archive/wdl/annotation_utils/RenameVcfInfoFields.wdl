version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow RenameVcfInfoFields {
    input {
        File vcf
        File vcf_idx
        Array[String] current_info_strings
        Array[String] replace_info_strings
        Array[String] replace_info_descriptions
        String prefix

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_rename
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
        call RenameInfoFields as RenameShardInfoFields {
            input:
                vcf = vcfs_to_process[i],
                vcf_idx = vcf_idxs_to_process[i],
                current_info_strings = current_info_strings,
                replace_info_strings = replace_info_strings,
                replace_info_descriptions = replace_info_descriptions,
                prefix = "~{prefix}.shard_~{i}.renamed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename
        }
    }

    if (defined(records_per_shard)) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = RenameShardInfoFields.renamed_vcf,
                vcf_idxs = RenameShardInfoFields.renamed_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.renamed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File renamed_vcf = select_first([ConcatVcfs.concat_vcf, RenameShardInfoFields.renamed_vcf[0]])
        File renamed_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, RenameShardInfoFields.renamed_vcf_idx[0]])
    }
}

task RenameInfoFields {
    input {
        File vcf
        File vcf_idx
        Array[String] current_info_strings
        Array[String] replace_info_strings
        Array[String] replace_info_descriptions
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        CURRENT_FILE="~{write_lines(current_info_strings)}"
        REPLACE_FILE="~{write_lines(replace_info_strings)}"
        DESC_FILE="~{write_lines(replace_info_descriptions)}"

        n_current=$(wc -l < "${CURRENT_FILE}")
        n_replace=$(wc -l < "${REPLACE_FILE}")
        n_desc=$(wc -l < "${DESC_FILE}")
        if [ "${n_current}" != "${n_replace}" ] || [ "${n_current}" != "${n_desc}" ]; then
            echo "ERROR: current_info_strings (${n_current}), replace_info_strings (${n_replace}), and replace_info_descriptions (${n_desc}) must have the same length" >&2
            exit 1
        fi

        paste \
            <(awk 'NF{print "INFO/"$0}' "${CURRENT_FILE}") \
            <(awk 'NF{print $0}'        "${REPLACE_FILE}") \
            > rename_annots.tsv

        bcftools annotate \
            --rename-annots rename_annots.tsv \
            -Oz -o renamed.vcf.gz \
            ~{vcf}

        bcftools view -h renamed.vcf.gz > new_header.txt

        paste "${REPLACE_FILE}" "${DESC_FILE}" | while IFS=$'\t' read -r new_id new_desc; do
            [ -z "${new_id}" ] && continue
            esc_desc=$(printf '%s' "${new_desc}" | sed -e 's/[\\&|]/\\&/g')
            sed -i -E "s|(^##INFO=<ID=${new_id},[^>]*Description=\")[^\"]*(\")|\1${esc_desc}\2|" new_header.txt
        done

        bcftools reheader -h new_header.txt -o ~{prefix}.vcf.gz renamed.vcf.gz
        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
