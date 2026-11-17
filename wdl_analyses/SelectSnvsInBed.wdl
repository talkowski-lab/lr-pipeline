version 1.0

import "Structs.wdl"

# Select SNVs from one or more VCF/BCF files that fall inside a BED of regions (e.g. high-confidence exome regions for
# PCA). Each input file (e.g. one per chromosome) is processed in parallel, then the per-file outputs are concatenated
# into one VCF. All input files must share the same samples in the same order. Records are kept when they overlap a BED
# interval and are SNVs. With biallelic_only, multiallelic sites are dropped (including those with spanning-deletion
# alleles); with pass_only, only FILTER=PASS records are kept (set false for pre-filtered files whose FILTER is ".").
# Genotypes are kept unchanged. Run the same workflow on case and control VCFs so both sides share the same site set.

workflow SelectSnvsInBed {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        File regions_bed
        String prefix

        Boolean pass_only = true
        Boolean biallelic_only = true

        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

        RuntimeAttr? runtime_attr_select_snvs_in_bed
        RuntimeAttr? runtime_attr_concat_snvs
    }

    scatter (i in range(length(vcfs))) {
        call SelectSnvs {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                regions_bed = regions_bed,
                pass_only = pass_only,
                biallelic_only = biallelic_only,
                prefix = "~{prefix}.shard_~{i}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_select_snvs_in_bed
        }
    }

    call ConcatSnvs {
        input:
            vcfs = SelectSnvs.snv_vcf,
            vcf_idxs = SelectSnvs.snv_vcf_idx,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_snvs
    }

    output {
        File snv_vcf = ConcatSnvs.snv_vcf
        File snv_vcf_idx = ConcatSnvs.snv_vcf_idx
        File snv_count = ConcatSnvs.snv_count
    }
}

task SelectSnvs {
    input {
        File vcf
        File vcf_idx
        File regions_bed
        Boolean pass_only
        Boolean biallelic_only
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String pass_arg = if pass_only then "-f PASS" else ""
    String biallelic_arg = if biallelic_only then "-m2 -M2" else ""

    command <<<
        set -euo pipefail

        bcftools view \
            -R ~{regions_bed} \
            -v snps \
            ~{pass_arg} \
            ~{biallelic_arg} \
            --threads 2 \
            -Oz -o ~{prefix}.snvs_in_bed.vcf.gz \
            ~{vcf}
        bcftools index -t ~{prefix}.snvs_in_bed.vcf.gz
    >>>

    output {
        File snv_vcf = "~{prefix}.snvs_in_bed.vcf.gz"
        File snv_vcf_idx = "~{prefix}.snvs_in_bed.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
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

task ConcatSnvs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat \
            -a \
            --file-list ~{write_lines(vcfs)} \
            --threads 2 \
            -Oz -o ~{prefix}.snvs_in_bed.vcf.gz
        bcftools index -t ~{prefix}.snvs_in_bed.vcf.gz
        bcftools index -n ~{prefix}.snvs_in_bed.vcf.gz > ~{prefix}.snvs_in_bed.count.txt
    >>>

    output {
        File snv_vcf = "~{prefix}.snvs_in_bed.vcf.gz"
        File snv_vcf_idx = "~{prefix}.snvs_in_bed.vcf.gz.tbi"
        File snv_count = "~{prefix}.snvs_in_bed.count.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
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
