version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow AnnotateTruvariRemap {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        Array[File] ref_bwa_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        String type_field = "allele_type"
        String type_ins = "ins"

        Int min_length
        Int max_length
        Int mm2_threshold
        Float cov_threshold

        String remap_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_ins_remap
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByArgs {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                include_args = "INFO/~{type_field}=\"~{type_ins}\"",
                extra_args = if single_contig then "" else "--regions ~{contig}",
                prefix = "~{prefix}.~{contig}.ins",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetVcfByArgs.subset_vcf,
                    vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.ins",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [SubsetVcfByArgs.subset_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [SubsetVcfByArgs.subset_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call InsRemap {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    min_length = min_length,
                    max_length = max_length,
                    mm2_threshold = mm2_threshold,
                    cov_threshold = cov_threshold,
                    ref_fa = ref_fa,
                    ref_bwa_idx = ref_bwa_idx,
                    prefix = "~{prefix}.~{contig}.remapped.shard_~{i}",
                    docker = remap_docker,
                    runtime_attr_override = runtime_attr_ins_remap
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = InsRemap.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.remap_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, InsRemap.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.remap_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_remap = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task InsRemap {
    input {
        File vcf
        File vcf_idx
        Int min_length
        Int max_length
        Int mm2_threshold
        Float cov_threshold
        File ref_fa
        Array[File] ref_bwa_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))

        mkdir ref_files
        mv ~{ref_fa} ref_files/
        for x in ~{sep=' ' ref_bwa_idx}; do
            mv $x ref_files/
        done

        ref_fa_basename=$(basename ~{ref_fa})

        truvari anno remap \
            -r ref_files/$ref_fa_basename \
            --min-length ~{min_length} \
            --max-length ~{max_length} \
            --mm2-threshold ~{mm2_threshold} \
            --cov-threshold ~{cov_threshold} \
            --threads ${N_THREADS} \
            -o ~{prefix}.vcf.gz \
            ~{vcf}

        bcftools query \
            -i 'INFO/remap_classification!="."' \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/remap_classification\t%INFO/remap_coords\t%INFO/remap_ori\t%INFO/remap_perc\n' \
            ~{prefix}.vcf.gz \
            > ~{prefix}.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(ref_fa, "GB") + size(ref_bwa_idx, "GB")) + 25,
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
