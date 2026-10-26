version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PALMERAssembly {
    meta {
        description: [
            "This workflow runs PALMER on a pair of aligned assembly haplotypes in order to generate MEI calls. It then convets the raw PALMER calls generated into a VCF, merges calls across the haplotypes to create a diploid VCF per haplotype and then finally integrates these into a final VCF containing multiple MEI types."
        ]
    }

    parameter_meta {
        bam_hap2: "Aligned haplotype 2 assembly."
        bai_hap2: "Index for `bam_hap2`."
        bam_hap1: "Aligned haplotype 1 assembly."
        bai_hap1: "Index for `bam_hap1`."
        override_palmer_calls_hap2: "Optional precomputed PALMER calls for haplotype 2, causing the workflow to bypass execution."
        override_palmer_tsd_files_hap2: "Optional precomputed PALMER TSD files for haplotype 2, causing the workflow to bypass execution."
        override_palmer_calls_hap1: "Optional precomputed PALMER calls for haplotype 1, causing the workflow to bypass execution."
        override_palmer_tsd_files_hap1: "Optional precomputed PALMER TSD files for haplotype 1, causing the workflow to bypass execution."
        ref_fa: "From references."
        ref_fai: "From references."
        contigs: "Contigs to run PALMER on."
        sample: "ID of the sample being processed."
        mode: "PALMER run mode."
        mei_types: "MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`."
        truvari_collapse_params: "Per-MEI-type Truvari parameters used when merging calls across haplotypes."
        palmer_hap2_calls: "Raw PALMER calls for haplotype 2, per MEI type."
        palmer_hap2_tsd_reads: "PALMER TSD reads for haplotype 2, per MEI type."
        palmer_hap2_vcfs: "Haplotype 2 PALMER VCFs, per MEI type."
        palmer_hap2_vcf_idxs: "Indexes for the haplotype 2 VCFs."
        palmer_hap1_calls: "Raw PALMER calls for haplotype 1, per MEI type."
        palmer_hap1_tsd_reads: "PALMER TSD reads for haplotype 1, per MEI type."
        palmer_hap1_vcfs: "Haplotype 1 PALMER VCFs, per MEI type."
        palmer_hap1_vcf_idxs: "Indexes for the haplotype 1 VCFs."
        palmer_diploid_vcfs: "Diploid PALMER VCFs merged across haplotypes, per MEI type."
        palmer_diploid_vcf_idxs: "Indexes for the diploid VCFs."
        palmer_combined_vcf: "Final VCF combining all MEI types."
        palmer_combined_vcf_idx: "Index for the combined VCF."
    }

    input {
        File? bam_hap2
        File? bai_hap2
        File? bam_hap1
        File? bai_hap1
        Array[File]? override_palmer_calls_hap2
        Array[File]? override_palmer_tsd_files_hap2
        Array[File]? override_palmer_calls_hap1
        Array[File]? override_palmer_tsd_files_hap1
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        String sample
        String mode
        Array[String] mei_types
        Array[String]? truvari_collapse_params

        String palmer_docker
        String utils_docker

        RuntimeAttr? runtime_attr_split_bam
        RuntimeAttr? runtime_attr_run_palmer
        RuntimeAttr? runtime_attr_merge_palmer_outputs
        RuntimeAttr? runtime_attr_palmer_to_vcf
        RuntimeAttr? runtime_attr_truvari_collapse
        RuntimeAttr? runtime_attr_concat
    }

    # Cromwell rejects an undefined optional entering a nested scatter, so these cross as one-element arrays
    Array[RuntimeAttr?] runtime_attr_run_palmer_wrapped = [runtime_attr_run_palmer]

    scatter (idx in range(length(mei_types))) {
        String mei_type = mei_types[idx]
        String collapse_params = if defined(truvari_collapse_params) then select_first([truvari_collapse_params])[idx] else "--pctsize 0.9 --pctovl 0.9 --pctseq 0.9 --refdist 500"

        # Call MEIs on haplotype 2, placed first in the genotype
        if (!defined(override_palmer_calls_hap2)) {
            call Helpers.SplitBam as SplitBamHap2 {
                input:
                    bam = select_first([bam_hap2]),
                    bai = select_first([bai_hap2]),
                    prefix = "~{prefix}.hap2",
                    contigs = contigs,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_split_bam
            }

            scatter (i in range(length(SplitBamHap2.bams))) {
                call Helpers.RunPALMERShard as RunPALMERShardHap2 {
                    input:
                        bam = SplitBamHap2.bams[i],
                        bai = SplitBamHap2.bais[i],
                        mode = mode,
                        mei_type = mei_type,
                        ref_fa = ref_fa,
                        prefix = "~{prefix}.shard_{i}.hap2",
                        docker = palmer_docker,
                        runtime_attr_override = runtime_attr_run_palmer_wrapped[0]
                }
            }

            call Helpers.MergePALMEROutputs as MergePALMEROutputsHap2 {
                input:
                    calls_shards = RunPALMERShardHap2.calls_shard,
                    tsd_reads_shards = RunPALMERShardHap2.tsd_reads_shard,
                    mei_type = mei_type,
                    prefix = "~{prefix}.merged_hap2",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_palmer_outputs
            }
        }

        File calls_file_hap2 = if defined (override_palmer_calls_hap2) then select_first([override_palmer_calls_hap2])[idx] else select_first([MergePALMEROutputsHap2.calls])
        File tsd_file_hap2 = if defined (override_palmer_tsd_files_hap2) then select_first([override_palmer_tsd_files_hap2])[idx] else select_first([MergePALMEROutputsHap2.tsd_reads])

        call Helpers.ConvertPALMERToVcf as ConvertPALMERToVcfHap2 {
            input:
                palmer_calls = calls_file_hap2,
                palmer_tsd_reads = tsd_file_hap2,
                mei_type = mei_type,
                sample = sample,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                haplotype = "1|0",
                prefix = "~{prefix}.~{mei_type}.hap2",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_palmer_to_vcf
        }

        # Call MEIs on haplotype 1, placed second in the genotype
        if (!defined(override_palmer_calls_hap1)) {
            call Helpers.SplitBam as SplitBamHap1 {
                input:
                    bam = select_first([bam_hap1]),
                    bai = select_first([bai_hap1]),
                    prefix = "~{prefix}.hap1",
                    contigs = contigs,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_split_bam
            }

            scatter (i in range(length(SplitBamHap1.bams))) {
                call Helpers.RunPALMERShard as RunPALMERShardHap1 {
                    input:
                        bam = SplitBamHap1.bams[i],
                        bai = SplitBamHap1.bais[i],
                        mode = mode,
                        mei_type = mei_type,
                        ref_fa = ref_fa,
                        prefix = "~{prefix}.shard_{i}.hap1",
                        docker = palmer_docker,
                        runtime_attr_override = runtime_attr_run_palmer_wrapped[0]
                }
            }

            call Helpers.MergePALMEROutputs as MergePALMEROutputsHap1 {
                input:
                    calls_shards = RunPALMERShardHap1.calls_shard,
                    tsd_reads_shards = RunPALMERShardHap1.tsd_reads_shard,
                    mei_type = mei_type,
                    prefix = "~{prefix}.merged_hap1",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_palmer_outputs
            }
        }

        File calls_file_hap1 = if defined (override_palmer_calls_hap1) then select_first([override_palmer_calls_hap1])[idx] else select_first([MergePALMEROutputsHap1.calls])
        File tsd_file_hap1 = if defined (override_palmer_tsd_files_hap1) then select_first([override_palmer_tsd_files_hap1])[idx] else select_first([MergePALMEROutputsHap1.tsd_reads])

        call Helpers.ConvertPALMERToVcf as ConvertPALMERToVcfHap1 {
            input:
                palmer_calls = calls_file_hap1,
                palmer_tsd_reads = tsd_file_hap1,
                mei_type = mei_type,
                sample = sample,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                haplotype = "0|1",
                prefix = "~{prefix}.~{mei_type}.hap1",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_palmer_to_vcf
        }

        # Collapse the two haplotype callsets into a diploid callset
        call TruvariCollapse {
            input:
                vcf_hap2 = ConvertPALMERToVcfHap2.vcf,
                vcf_hap2_idx = ConvertPALMERToVcfHap2.vcf_idx,
                vcf_hap1 = ConvertPALMERToVcfHap1.vcf,
                vcf_hap1_idx = ConvertPALMERToVcfHap1.vcf_idx,
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
        Array[File] palmer_hap2_calls = calls_file_hap2
        Array[File] palmer_hap2_tsd_reads = tsd_file_hap2
        Array[File] palmer_hap2_vcfs = ConvertPALMERToVcfHap2.vcf
        Array[File] palmer_hap2_vcf_idxs = ConvertPALMERToVcfHap2.vcf_idx

        Array[File] palmer_hap1_calls = calls_file_hap1
        Array[File] palmer_hap1_tsd_reads = tsd_file_hap1
        Array[File] palmer_hap1_vcfs = ConvertPALMERToVcfHap1.vcf
        Array[File] palmer_hap1_vcf_idxs = ConvertPALMERToVcfHap1.vcf_idx

        Array[File] palmer_diploid_vcfs = TruvariCollapse.diploid_vcf
        Array[File] palmer_diploid_vcf_idxs = TruvariCollapse.diploid_vcf_idx

        File palmer_combined_vcf = ConcatVcfs.concat_vcf
        File palmer_combined_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task TruvariCollapse {
    input {
        File vcf_hap2
        File vcf_hap2_idx
        File vcf_hap1
        File vcf_hap1_idx
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
            ~{vcf_hap2} \
            ~{vcf_hap1} \
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
        disk_gb: 3 * ceil(size(vcf_hap2, "GB") + size(vcf_hap1, "GB")) + 10,
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
