version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSVAnnotate {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Int min_length
        File coding_gtf
        File noncoding_bed

        String utils_docker
        String gatk_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_convert_symbolic
        RuntimeAttr? runtime_attr_annotate_func
        RuntimeAttr? runtime_attr_concat_unannotated
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_annotated
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_revert_symbolic
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByLength {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                min_length = min_length,
                extra_args = if single_contig then "-G" else "-G --regions ~{contig}",
                prefix = "~{prefix}.~{contig}.subset",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetVcfByLength.subset_vcf,
                    vcf_idx = SubsetVcfByLength.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [SubsetVcfByLength.subset_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [SubsetVcfByLength.subset_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call Helpers.ConvertToSymbolic {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    move_all_dups = true,
                    prefix = "~{prefix}.~{contig}.converted.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_convert_symbolic
            }

            call AnnotateFunctionalConsequences {
                input:
                    vcf = ConvertToSymbolic.processed_vcf,
                    vcf_idx = ConvertToSymbolic.processed_vcf_idx,
                    noncoding_bed = noncoding_bed,
                    coding_gtf = coding_gtf,
                    prefix = "~{prefix}.~{contig}.functionally_annotated.shard_~{i}",
                    docker = gatk_docker,
                    runtime_attr_override = runtime_attr_annotate_func
            }

            call Helpers.RevertSymbolicAlleles {
                input:
                    annotated_vcf = AnnotateFunctionalConsequences.anno_vcf,
                    annotated_vcf_idx = AnnotateFunctionalConsequences.anno_vcf_idx,
                    original_vcf = vcfs_to_process[i],
                    original_vcf_idx = vcf_idxs_to_process[i],
                    prefix = "~{prefix}.~{contig}.reverted.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_revert_symbolic
            }

            call Helpers.ExtractVcfAnnotations {
                input:
                    vcf = RevertSymbolicAlleles.reverted_vcf,
                    vcf_idx = RevertSymbolicAlleles.reverted_vcf_idx,
                    original_vcf = vcfs_to_process[i],
                    original_vcf_idx = vcf_idxs_to_process[i],
                    prefix = "~{prefix}.~{contig}.shard_~{i}",
                    docker = utils_docker
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = ExtractVcfAnnotations.annotations_tsv,
                    sort_output = true,
                    prefix = "~{prefix}.~{contig}.svannotate_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }

            call Helpers.MergeHeaderLines as MergeShardHeaders {
                input:
                    header_files = ExtractVcfAnnotations.annotations_header,
                    prefix = "~{prefix}.~{contig}.svannotate_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, ExtractVcfAnnotations.annotations_tsv[0]])
        File final_annotations_header = select_first([MergeShardHeaders.merged_header, ExtractVcfAnnotations.annotations_header[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.svannotate_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_annotated
        }

        call Helpers.MergeHeaderLines {
            input:
                header_files = final_annotations_header,
                prefix = "~{prefix}.svannotate_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    output {
        File annotations_tsv_svannotate = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
        File annotations_header_svannotate = select_first([MergeHeaderLines.merged_header, final_annotations_header[0]])
    }
}

task AnnotateFunctionalConsequences {
    input {
        File vcf
        File vcf_idx
        File noncoding_bed
        File coding_gtf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
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
    Int java_mem_mb = 1000 * ceil(0.8 * select_first([runtime_attr.mem_gb, default_attr.mem_gb]))

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{java_mem_mb}m" SVAnnotate \
            -V ~{vcf} \
            --non-coding-bed ~{noncoding_bed} \
            --protein-coding-gtf ~{coding_gtf} \
            -O ~{prefix}.vcf.gz
    >>>

    output {
        File anno_vcf = "~{prefix}.vcf.gz"
        File anno_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

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
