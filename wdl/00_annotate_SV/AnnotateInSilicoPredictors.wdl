version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateInSilicoPredictors {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht
        String annotate_in_silico_predictors_script = "https://raw.githubusercontent.com/talkowski-lab/lr-annotation/main/scripts/miscellaneous/annotate_insilico_predictors.py"
        String genome_build = "GRCh38"

        String hail_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_vcf
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

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.in_silico_predictors",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call AnnotateInSilicoPredictorsTask {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    cadd_ht = cadd_ht,
                    pangolin_ht = pangolin_ht,
                    phylop_ht = phylop_ht,
                    revel_ht = revel_ht,
                    spliceai_ht = spliceai_ht,
                    annotate_in_silico_predictors_script = annotate_in_silico_predictors_script,
                    genome_build = genome_build,
                    prefix = "~{prefix}.~{contig}.in_silico_predictors.shard_~{i}",
                    docker = hail_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = AnnotateInSilicoPredictorsTask.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.in_silico_predictors",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, AnnotateInSilicoPredictorsTask.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.in_silico_predictors",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    output {
        File annotations_tsv_insilico = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task AnnotateInSilicoPredictorsTask {
    input {
        File vcf
        File vcf_idx
        String cadd_ht
        String pangolin_ht
        String phylop_ht
        String revel_ht
        String spliceai_ht
        String annotate_in_silico_predictors_script
        String genome_build = "GRCh38"
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        curl ~{annotate_in_silico_predictors_script} > annotate.py

        python3 annotate.py \
            --build ~{genome_build} \
            --cores ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb])} \
            --cadd_ht ~{cadd_ht} \
            --pangolin_ht ~{pangolin_ht} \
            --phylop_ht ~{phylop_ht} \
            --revel_ht ~{revel_ht} \
            --spliceai_ht ~{spliceai_ht} \
            --vcf ~{vcf} \
            --output_tsv ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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