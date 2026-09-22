version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow Whatshap {
    meta {
        description: [
            "This tool haplotags a sample's BAM against a phased VCF using WhatsHap (https://github.com/whatshap/whatshap), per contig, then merges the tagged reads into a single BAM. It outputs the haplotagged BAM and per-contig haplotag read lists."
        ]
    }

    parameter_meta {
        bam: "Aligned reads to haplotag."
        bai: "Index for `bam`."
        phased_vcf: "Phased VCF used to assign haplotypes."
        phased_vcf_idx: "Index for `phased_vcf`."
        ref_fa: "From references."
        ref_fai: "From references."
        contigs: "Contigs to process."
        extra_args: "Additional arguments passed to WhatsHap."
        haplotagged_bam: "Haplotagged BAM."
        haplotagged_bai: "Index for the haplotagged BAM."
        haplotag_lists: "Per-contig haplotag read assignments."
    }

    input {
        File bam
        File bai
        File phased_vcf
        File phased_vcf_idx
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        String? extra_args

        String utils_docker
        String whatshap_docker

        RuntimeAttr? runtime_attr_subset_bam
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_haplotag
        RuntimeAttr? runtime_attr_merge_bams
    }

    scatter (contig in contigs) {
        call Helpers.SubsetBamToContig {
            input:
                bam = bam,
                bai = bai,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_bam
        }

        call Helpers.SubsetVcfToContig {
            input:
                vcf = phased_vcf,
                vcf_idx = phased_vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}.phased",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Haplotag {
            input:
                bam = SubsetBamToContig.contig_bam,
                bai = SubsetBamToContig.contig_bai,
                phased_vcf = SubsetVcfToContig.subset_vcf,
                phased_vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                extra_args = extra_args,
                prefix = "~{prefix}.~{contig}.haplotagged",
                docker = whatshap_docker,
                runtime_attr_override = runtime_attr_haplotag
        }
    }

    call Helpers.MergeBams {
        input:
            bams = Haplotag.haplotagged_bam,
            bais = Haplotag.haplotagged_bai,
            prefix = "~{prefix}.haplotagged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_bams
    }

    output {
        File haplotagged_bam = MergeBams.merged_bam
        File haplotagged_bai = MergeBams.merged_bam_idx
        Array[File] haplotag_lists = Haplotag.haplotag_list
    }
}

task Haplotag {
    input {
        File bam
        File bai
        File phased_vcf
        File phased_vcf_idx
        File ref_fa
        File ref_fai
        String? extra_args
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        whatshap haplotag \
            -o ~{prefix}.out.bam \
            --reference ~{ref_fa} \
            --output-threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --output-haplotag-list ~{prefix}.haplotag_list.tsv \
            ~{if defined(extra_args) then extra_args else ""} \
            ~{phased_vcf} \
            ~{bam}

        samtools sort \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            -o ~{prefix}.bam \
            ~{prefix}.out.bam

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File haplotagged_bam = "~{prefix}.bam"
        File haplotagged_bai = "~{prefix}.bam.bai"
        File haplotag_list = "~{prefix}.haplotag_list.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 5 * ceil(size(bam, "GB") + size(ref_fa, "GB")) + 25,
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
