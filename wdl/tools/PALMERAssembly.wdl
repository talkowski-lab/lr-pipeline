version 1.0

import "../utils/Helpers.wdl"
import "../utils/PALMERTasks.wdl"
import "../utils/Structs.wdl"

workflow PALMERAssembly {
    input {
        File? bam_pat
        File? bai_pat
        File? bam_mat
        File? bai_mat
        Array[File]? override_palmer_calls_pat
        Array[File]? override_palmer_tsd_files_pat
        Array[File]? override_palmer_calls_mat
        Array[File]? override_palmer_tsd_files_mat
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        String sample
        String mode
        Array[String] mei_types
        Array[String]? truvari_collapse_params

        String annotate_palmer_docker
        String palmer_docker
        String utils_docker

        RuntimeAttr? runtime_attr_split_bam
        RuntimeAttr? runtime_attr_run_palmer
        RuntimeAttr? runtime_attr_merge_palmer_outputs
        RuntimeAttr? runtime_attr_palmer_to_vcf
        RuntimeAttr? runtime_attr_truvari_collapse
        RuntimeAttr? runtime_attr_concat
    }

    scatter (idx in range(length(mei_types))) {
        String mei_type = mei_types[idx]
        String collapse_params = if defined(truvari_collapse_params) then select_first([truvari_collapse_params])[idx] else "--pctsize 0.9 --pctovl 0.9 --pctseq 0.9 --refdist 500"

        # Paternal haplotype
        if (!defined(override_palmer_calls_pat)) {
            call Helpers.SplitBam as SplitBamPat {
                input:
                    bam = select_first([bam_pat]),
                    bai = select_first([bai_pat]),
                    prefix = "~{prefix}.pat",
                    contigs = contigs,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_split_bam
            }

            scatter (i in range(length(SplitBamPat.bams))) {
                call PALMERTasks.RunPALMERShard as RunPALMERShardPat {
                    input:
                        bam = SplitBamPat.bams[i],
                        bai = SplitBamPat.bais[i],
                        mode = mode,
                        mei_type = mei_type,
                        ref_fa = ref_fa,
                        prefix = "~{prefix}.shard_{i}.pat",
                        docker = palmer_docker,
                        runtime_attr_override = runtime_attr_run_palmer
                }
            }

            call PALMERTasks.MergePALMEROutputs as MergePALMEROutputsPat {
                input:
                    calls_shards = RunPALMERShardPat.calls_shard,
                    tsd_reads_shards = RunPALMERShardPat.tsd_reads_shard,
                    mei_type = mei_type,
                    prefix = "~{prefix}.merged_pat",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_palmer_outputs
            }
        }

        File calls_file_pat = if defined (override_palmer_calls_pat) then select_first([override_palmer_calls_pat])[idx] else select_first([MergePALMEROutputsPat.calls])
        File tsd_file_pat = if defined (override_palmer_tsd_files_pat) then select_first([override_palmer_tsd_files_pat])[idx] else select_first([MergePALMEROutputsPat.tsd_reads])

        call PALMERTasks.ConvertPALMERToVcf as ConvertPALMERToVcfPat {
            input:
                palmer_calls = calls_file_pat,
                palmer_tsd_reads = tsd_file_pat,
                mei_type = mei_type,
                sample = sample,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                haplotype = "1|0",
                prefix = "~{prefix}.~{mei_type}.pat",
                docker = annotate_palmer_docker,
                runtime_attr_override = runtime_attr_palmer_to_vcf
        }

        # Maternal haplotype
        if (!defined(override_palmer_calls_mat)) {
            call Helpers.SplitBam as SplitBamMat {
                input:
                    bam = select_first([bam_mat]),
                    bai = select_first([bai_mat]),
                    prefix = "~{prefix}.mat",
                    contigs = contigs,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_split_bam
            }

            scatter (i in range(length(SplitBamMat.bams))) {
                call PALMERTasks.RunPALMERShard as RunPALMERShardMat {
                    input:
                        bam = SplitBamMat.bams[i],
                        bai = SplitBamMat.bais[i],
                        mode = mode,
                        mei_type = mei_type,
                        ref_fa = ref_fa,
                        prefix = "~{prefix}.shard_{i}.mat",
                        docker = palmer_docker,
                        runtime_attr_override = runtime_attr_run_palmer
                }
            }

            call PALMERTasks.MergePALMEROutputs as MergePALMEROutputsMat {
                input:
                    calls_shards = RunPALMERShardMat.calls_shard,
                    tsd_reads_shards = RunPALMERShardMat.tsd_reads_shard,
                    mei_type = mei_type,
                    prefix = "~{prefix}.merged_mat",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_palmer_outputs
            }
        }

        File calls_file_mat = if defined (override_palmer_calls_mat) then select_first([override_palmer_calls_mat])[idx] else select_first([MergePALMEROutputsMat.calls])
        File tsd_file_mat = if defined (override_palmer_tsd_files_mat) then select_first([override_palmer_tsd_files_mat])[idx] else select_first([MergePALMEROutputsMat.tsd_reads])

        call PALMERTasks.ConvertPALMERToVcf as ConvertPALMERToVcfMat {
            input:
                palmer_calls = calls_file_mat,
                palmer_tsd_reads = tsd_file_mat,
                mei_type = mei_type,
                sample = sample,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                haplotype = "0|1",
                prefix = "~{prefix}.~{mei_type}.mat",
                docker = annotate_palmer_docker,
                runtime_attr_override = runtime_attr_palmer_to_vcf
        }

        # Merge haplotypes
        call TruvariCollapse {
            input:
                vcf_pat = ConvertPALMERToVcfPat.vcf,
                vcf_pat_idx = ConvertPALMERToVcfPat.vcf_idx,
                vcf_mat = ConvertPALMERToVcfMat.vcf,
                vcf_mat_idx = ConvertPALMERToVcfMat.vcf_idx,
                collapse_params = collapse_params,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{mei_type}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_truvari_collapse
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = TruvariCollapse.diploid_vcf,
            vcf_idxs = TruvariCollapse.diploid_vcf_idx,
            allow_overlaps = true,
            naive = false,
            prefix = "~{prefix}.concat",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        Array[File] palmer_pat_calls = calls_file_pat
        Array[File] palmer_pat_tsd_reads = tsd_file_pat
        Array[File] palmer_pat_vcfs = ConvertPALMERToVcfPat.vcf
        Array[File] palmer_pat_vcf_idxs = ConvertPALMERToVcfPat.vcf_idx

        Array[File] palmer_mat_calls = calls_file_mat
        Array[File] palmer_mat_tsd_reads = tsd_file_mat
        Array[File] palmer_mat_vcfs = ConvertPALMERToVcfMat.vcf
        Array[File] palmer_mat_vcf_idxs = ConvertPALMERToVcfMat.vcf_idx

        Array[File] palmer_diploid_vcfs = TruvariCollapse.diploid_vcf
        Array[File] palmer_diploid_vcf_idxs = TruvariCollapse.diploid_vcf_idx

        File palmer_combined_vcf = ConcatVcfs.concat_vcf
        File palmer_combined_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task TruvariCollapse {
    input {
        File vcf_pat
        File vcf_pat_idx
        File vcf_mat
        File vcf_mat_idx
        String collapse_params
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat \
            -a \
            ~{vcf_pat} \
            ~{vcf_mat} \
        | bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o combined.vcf.gz

        tabix combined.vcf.gz

        truvari collapse \
            --reference ~{ref_fa} \
            -i combined.vcf.gz \
            -o ~{prefix}.merged.vcf.gz \
            -c ~{prefix}.collapsed.vcf.gz \
            --hap \
            ~{collapse_params}

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o ~{prefix}.merged.sorted.vcf.gz \
            ~{prefix}.merged.vcf.gz

        tabix -p vcf ~{prefix}.merged.sorted.vcf.gz
    >>>

    output {
        File diploid_vcf = "~{prefix}.merged.sorted.vcf.gz"
        File diploid_vcf_idx = "~{prefix}.merged.sorted.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf_pat, "GB") + size(vcf_mat, "GB")) + 10,
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
