version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateIndelTRs {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Int min_tandem_repeat_length = 9
        Int min_repeats = 3
        Int min_repeat_unit_length = 1

        String stranalysis_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByArgs {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                exclude_args = "INFO/allele_type=\"trv\" || INFO/TR_ENVELOPED=1",
                extra_args = if single_contig then "" else "--regions " + contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetVcfByArgs.subset_vcf,
                    vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.tr_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [SubsetVcfByArgs.subset_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [SubsetVcfByArgs.subset_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call RunFilterVcfToTRs {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}.~{contig}.tr_annotations.shard_~{i}",
                    min_tandem_repeat_length = min_tandem_repeat_length,
                    min_repeats = min_repeats,
                    min_repeat_unit_length = min_repeat_unit_length,
                    docker = stranalysis_docker,
                    runtime_attr_override = runtime_attr_filter
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = RunFilterVcfToTRs.tr_annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.tr_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, RunFilterVcfToTRs.tr_annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.tr_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_trs = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task RunFilterVcfToTRs {
    input {
        File vcf
        File vcf_idx
        File ref_fa
        File ref_fai
        Int min_tandem_repeat_length
        Int min_repeats
        Int min_repeat_unit_length
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        python3 -m str_analysis.filter_vcf_to_tandem_repeats catalog \
            -R ~{ref_fa} \
            --output-prefix ~{prefix} \
            --min-tandem-repeat-length ~{min_tandem_repeat_length} \
            --min-repeats ~{min_repeats} \
            --min-repeat-unit-length ~{min_repeat_unit_length} \
            --write-vcf \
            --trf-executable-path $(which trf) \
            --trf-threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{vcf}

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t1\n' \
            ~{prefix}.tandem_repeats.vcf.gz \
            > ~{prefix}.tsv
    >>>

    output {
        File tr_annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(ref_fa, "GB")) + 20,
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
