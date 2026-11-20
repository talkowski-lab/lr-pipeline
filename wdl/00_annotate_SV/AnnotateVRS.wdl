version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateVRS {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        File seqrepo_tar

        String utils_docker
        String vrs_docker

        RuntimeAttr? runtime_attr_drop_fields
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate_vrs
        RuntimeAttr? runtime_attr_extract
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

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.vrs",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call Helpers.DropVcfFields {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    drop_fields = "INFO/VRS_Allele_IDs,INFO/VRS_Error,INFO/VRS_Starts,INFO/VRS_Ends,INFO/VRS_States,INFO/VRS_Lengths,INFO/VRS_RepeatSubunitLengths",
                    prefix = "~{prefix}.~{contig}.vrs_cleared.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_drop_fields
            }

            call AnnotateVcfWithVRS {
                input:
                    vcf = DropVcfFields.dropped_vcf,
                    vcf_idx = DropVcfFields.dropped_vcf_idx,
                    seqrepo_tar = seqrepo_tar,
                    prefix = "~{prefix}.~{contig}.vrs.shard_~{i}",
                    docker = vrs_docker,
                    runtime_attr_override = runtime_attr_annotate_vrs
            }

            call ExtractVRSAnnotations {
                input:
                    vcf = AnnotateVcfWithVRS.annotated_vcf,
                    vcf_idx = AnnotateVcfWithVRS.annotated_vcf_idx,
                    prefix = "~{prefix}.~{contig}.vrs_annotations.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = ExtractVRSAnnotations.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.vrs_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, ExtractVRSAnnotations.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.vrs_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_vrs = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task AnnotateVcfWithVRS {
    input {
        File vcf
        File vcf_idx
        File seqrepo_tar
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir seqrepo_data
        tar -xzf ~{seqrepo_tar} -C seqrepo_data --strip-components 1

        SEQREPO_PATH=$(find $(pwd)/seqrepo_data -name "aliases.sqlite3" | xargs dirname)

        vrs-annotate vcf \
            --dataproxy-uri="seqrepo+file://${SEQREPO_PATH}" \
            --vcf-out ~{prefix}.vcf.gz \
            --vrs-attributes \
            ~{vcf}

        rm -rf seqrepo_data

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(seqrepo_tar, "GB")) + 20,
        boot_disk_gb: 50,
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

task ExtractVRSAnnotations {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/VRS_Allele_IDs\t%INFO/VRS_Error\t%INFO/VRS_Starts\t%INFO/VRS_Ends\t%INFO/VRS_States\t%INFO/VRS_Lengths\t%INFO/VRS_RepeatSubunitLengths\n' \
            ~{vcf} \
            > ~{prefix}.annotations.tsv
    >>>

    output {
        File annotations_tsv = "~{prefix}.annotations.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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
