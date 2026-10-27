version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MethylationProfiling {
    meta {
        description: [
            "This tool generates a sample's CpG methylation pileups using pb-CpG-tools (https://github.com/PacificBiosciences/pb-CpG-tools), producing combined and per-haplotype methylation BED tracks.",
            "The reads are first haplotagged per contig with WhatsHap (https://github.com/whatshap/whatshap) against the sample's phased heterozygous biallelic SNVs, so haplotype 1 follows the first GT allele; run on the sample's backbone-phased VCF, that keeps each haplotype on one parent across the contig. The haplotagged contigs are merged before the pileup."
        ]
    }

    parameter_meta {
        bam: "Aligned reads."
        bai: "Index for `bam`."
        phased_vcf: "The sample's phased VCF, ideally from backbone phasing."
        phased_vcf_idx: "Index for `phased_vcf`."
        contigs: "Contigs to haplotag and profile."
        ref_fa: "From references."
        ref_fai: "From references."
        cpg_combined_bed: "Combined methylation pileup BED."
        cpg_combined_bed_idx: "Index for the combined BED."
        cpg_hap1_bed: "Haplotype 1 methylation pileup BED."
        cpg_hap1_bed_idx: "Index for the haplotype 1 BED."
        cpg_hap2_bed: "Haplotype 2 methylation pileup BED."
        cpg_hap2_bed_idx: "Index for the haplotype 2 BED."
    }

    input {
        File bam
        File bai
        File phased_vcf
        File phased_vcf_idx
        Array[String] contigs
        File ref_fa
        File ref_fai
        String prefix

        String utils_docker
        String whatshap_docker
        String cpg_docker

        RuntimeAttr? runtime_attr_subset_bam
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_haplotag
        RuntimeAttr? runtime_attr_merge_bams
        RuntimeAttr? runtime_attr_cpg_pileup
    }

    scatter (contig in contigs) {
        call Helpers.SubsetBamToRegions {
            input:
                bam = bam,
                bai = bai,
                region_file = write_lines([contig]),
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_bam
        }

        call Helpers.SubsetVcfByArgs {
            input:
                vcf = phased_vcf,
                vcf_idx = phased_vcf_idx,
                include_args = 'GT="het"',
                extra_args = "--regions ~{contig} -v snps -m2 -M2",
                prefix = "~{prefix}.~{contig}.phased_snvs",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Haplotag {
            input:
                bam = SubsetBamToRegions.subset_bam,
                bai = SubsetBamToRegions.subset_bai,
                phased_vcf = SubsetVcfByArgs.subset_vcf,
                phased_vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
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

    call CpgPileup {
        input:
            bam = MergeBams.merged_bam,
            bai = MergeBams.merged_bam_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = cpg_docker,
            runtime_attr_override = runtime_attr_cpg_pileup
    }

    output {
        File cpg_combined_bed = CpgPileup.pileup_combined_bed
        File cpg_combined_bed_idx = CpgPileup.pileup_combined_bed_idx
        File cpg_hap1_bed = CpgPileup.pileup_hap1_bed
        File cpg_hap1_bed_idx = CpgPileup.pileup_hap1_bed_idx
        File cpg_hap2_bed = CpgPileup.pileup_hap2_bed
        File cpg_hap2_bed_idx = CpgPileup.pileup_hap2_bed_idx
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
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # The VCF holds one sample, so read group sample names need not match it; existing HP/PS tags are replaced
        whatshap haplotag \
            -o ~{prefix}.bam \
            --reference ~{ref_fa} \
            --ignore-read-groups \
            --output-threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{phased_vcf} \
            ~{bam}

        samtools index \
            -@ ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix}.bam
    >>>

    output {
        File haplotagged_bam = "~{prefix}.bam"
        File haplotagged_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 8,
        disk_gb: 3 * ceil(size(bam, "GB") + size(ref_fa, "GB")) + 25,
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

task CpgPileup {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        aligned_bam_to_cpg_scores \
            --bam ~{bam} \
            --ref ~{ref_fa} \
            --output-prefix ~{prefix} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --pileup-mode model \
            --min-coverage 4 \
            --min-mapq 10
    >>>

    output {
        File pileup_combined_bed = "~{prefix}.combined.bed.gz"
        File pileup_combined_bed_idx = "~{prefix}.combined.bed.gz.tbi"
        File pileup_hap1_bed = "~{prefix}.hap1.bed.gz"
        File pileup_hap1_bed_idx = "~{prefix}.hap1.bed.gz.tbi"
        File pileup_hap2_bed = "~{prefix}.hap2.bed.gz"
        File pileup_hap2_bed_idx = "~{prefix}.hap2.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 12,
        disk_gb: 3 * ceil(size(bam, "GB") + size(ref_fa, "GB")) + 25,
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
