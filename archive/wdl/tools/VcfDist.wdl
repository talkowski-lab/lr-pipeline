version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow VcfDist {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        File ref_fa
        Array[String] contigs
        String prefix

        File? bed_regions
        String? mode
        Float? threshold
        String? vcfdist_args

        String utils_docker
        String vcfdist_docker

        RuntimeAttr? runtime_attr_subset_eval
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_vcfdist
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig as SubsetVcfEval {
            input:
                vcf = vcf_eval,
                vcf_idx = vcf_eval_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_eval
        }

        call Helpers.SubsetVcfToContig as SubsetVcfTruth {
            input:
                vcf = vcf_truth,
                vcf_idx = vcf_truth_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_truth
        }

        call RunVcfDist {
            input:
                vcf_eval = SubsetVcfEval.subset_vcf,
                vcf_eval_idx = SubsetVcfEval.subset_vcf_idx,
                vcf_truth = SubsetVcfTruth.subset_vcf,
                vcf_truth_idx = SubsetVcfTruth.subset_vcf_idx,
                ref_fa = ref_fa,
                bed_regions = bed_regions,
                mode = mode,
                threshold = threshold,
                vcfdist_args = vcfdist_args,
                prefix = "~{prefix}.~{contig}",
                docker = vcfdist_docker,
                runtime_attr_override = runtime_attr_vcfdist
        }
    }

    output {
        Array[File] vcfdist_phasing_summary_tsv = RunVcfDist.phasing_summary_tsv
        Array[File] vcfdist_switchflips_tsv = RunVcfDist.switchflips_tsv
        Array[File] vcfdist_precision_recall_tsv = RunVcfDist.precision_recall_tsv
        Array[File] vcfdist_precision_recall_summary_tsv = RunVcfDist.precision_recall_summary_tsv
        Array[File] vcfdist_phase_blocks_tsv = RunVcfDist.phase_blocks_tsv
        Array[File] vcfdist_superclusters_tsv = RunVcfDist.superclusters_tsv
        Array[File] vcfdist_query_tsv = RunVcfDist.query_tsv
        Array[File] vcfdist_truth_tsv = RunVcfDist.truth_tsv
        Array[File] vcfdist_summary_vcf = RunVcfDist.summary_vcf
    }
}

task RunVcfDist {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        File ref_fa
        File? bed_regions
        String? mode
        Float? threshold
        String? vcfdist_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        
        mkdir -p results
        
        vcfdist \
            ~{vcf_eval} \
            ~{vcf_truth} \
            ~{ref_fa} \
            -p "results/~{prefix}." \
            ~{if defined(bed_regions) then "-b " + bed_regions else ""} \
            ~{if defined(mode) then "-m " + mode else ""} \
            ~{if defined(threshold) then "-t " + threshold else ""} \
            ~{if defined(vcfdist_args) then vcfdist_args else ""}
    >>>

    output {
        File phasing_summary_tsv = "results/~{prefix}.phasing-summary.tsv"
        File switchflips_tsv = "results/~{prefix}.switchflips.tsv"
        File precision_recall_tsv = "results/~{prefix}.precision-recall.tsv"
        File precision_recall_summary_tsv = "results/~{prefix}.precision-recall-summary.tsv"
        File phase_blocks_tsv = "results/~{prefix}.phase-blocks.tsv"
        File superclusters_tsv = "results/~{prefix}.superclusters.tsv"
        File query_tsv = "results/~{prefix}.query.tsv"
        File truth_tsv = "results/~{prefix}.truth.tsv"
        File summary_vcf = "results/~{prefix}.summary.vcf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 32,
        disk_gb: 2 * ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) + 10,
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