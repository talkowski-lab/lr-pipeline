version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PALMERDiploid {
    input {
        File? bam
        File? bai
        Array[File]? override_palmer_calls
        Array[File]? override_palmer_tsd_files
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        String sample
        String mode
        Array[String] mei_types

        String annotate_palmer_docker
        String palmer_docker
        String utils_docker

        RuntimeAttr? runtime_attr_split_bam
        RuntimeAttr? runtime_attr_run_palmer
        RuntimeAttr? runtime_attr_merge_palmer_outputs
        RuntimeAttr? runtime_attr_palmer_to_vcf
        RuntimeAttr? runtime_attr_concat
    }

    scatter (idx in range(length(mei_types))) {
        String mei_type = mei_types[idx]

        if (!defined(override_palmer_calls)) {
            call Helpers.SplitBam {
                input:
                    bam = select_first([bam]),
                    bai = select_first([bai]),
                    prefix = "~{prefix}.split",
                    contigs = contigs,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_split_bam
            }

            scatter (i in range(length(SplitBam.bams))) {
                call Helpers.RunPALMERShard {
                    input:
                        bam = SplitBam.bams[i],
                        bai = SplitBam.bais[i],
                        mode = mode,
                        mei_type = mei_type,
                        ref_fa = ref_fa,
                        prefix = "~{prefix}.shard_{i}",
                        docker = palmer_docker,
                        runtime_attr_override = runtime_attr_run_palmer
                }
            }

            call Helpers.MergePALMEROutputs {
                input:
                    calls_shards = RunPALMERShard.calls_shard,
                    tsd_reads_shards = RunPALMERShard.tsd_reads_shard,
                    mei_type = mei_type,
                    prefix = "~{prefix}.merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_palmer_outputs
            }
        }

        if (defined(override_palmer_calls)) {
            call AddMeiTypeColumn as AddMeiTypeToCallsOverride {
                input:
                    input_file = select_first([override_palmer_calls])[idx],
                    mei_type = mei_type,
                    file_type = "calls",
                    prefix = "~{prefix}.mei_type_added",
                    docker = utils_docker
            }
        }

        if (defined(override_palmer_tsd_files)) {
            call AddMeiTypeColumn as AddMeiTypeToTsdOverride {
                input:
                    input_file = select_first([override_palmer_tsd_files])[idx],
                    mei_type = mei_type,
                    file_type = "tsd_reads",
                    prefix = "~{prefix}.mei_type_added",
                    docker = utils_docker
            }
        }

        File calls_file = if defined (override_palmer_calls) then select_first([AddMeiTypeToCallsOverride.output_file]) else select_first([MergePALMEROutputs.calls])
        File tsd_file = if defined (override_palmer_tsd_files) then select_first([AddMeiTypeToTsdOverride.output_file]) else select_first([MergePALMEROutputs.tsd_reads])

        call Helpers.ConvertPALMERToVcf {
            input:
                palmer_calls = calls_file,
                palmer_tsd_reads = tsd_file,
                mei_type = mei_type,
                sample = sample,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                haplotype = "1/1",
                prefix = "~{prefix}.~{mei_type}",
                docker = annotate_palmer_docker,
                runtime_attr_override = runtime_attr_palmer_to_vcf
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = ConvertPALMERToVcf.vcf,
            vcf_idxs = ConvertPALMERToVcf.vcf_idx,
            allow_overlaps = true,
            naive = false,
            prefix = "~{prefix}.concat",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        Array[File] palmer_calls = calls_file
        Array[File] palmer_tsd_reads = tsd_file
        Array[File] palmer_diploid_vcfs = ConvertPALMERToVcf.vcf
        Array[File] palmer_diploid_vcf_idxs = ConvertPALMERToVcf.vcf_idx

        File palmer_combined_vcf = ConcatVcfs.concat_vcf
        File palmer_combined_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task AddMeiTypeColumn {
    input {
        File input_file
        String mei_type
        String file_type
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        sed "s/$/\t~{mei_type}/" ~{input_file} > ~{prefix}_~{mei_type}_~{file_type}.txt
    >>>

    output {
        File output_file = "~{prefix}_~{mei_type}_~{file_type}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(input_file, "GB")) + 10,
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
