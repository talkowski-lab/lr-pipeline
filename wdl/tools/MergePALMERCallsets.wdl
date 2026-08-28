version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergePALMERCallsets {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call MergeVCFs {
            input:
                vcfs = vcfs,
                vcf_idxs = vcf_idxs,
                contig = contig,
                prefix = "~{prefix}.{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_vcfs
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = MergeVCFs.merged_vcf,
            vcf_idxs = MergeVCFs.merged_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.palmer_merged",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File palmer_merged_vcf = ConcatVcfs.concat_vcf
        File palmer_merged_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task MergeVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir inputs

        vcf_files=(~{sep=" " vcfs})
        vcf_idxs=(~{sep=" " vcf_idxs})

        for i in ${!vcf_files[@]}; do
            ln -s ${vcf_files[$i]} inputs/input_$i.vcf.gz
            ln -s ${vcf_idxs[$i]} inputs/input_$i.vcf.gz.tbi
        done

        bcftools merge \
            -r ~{contig} \
            --missing-to-ref \
            -Oz -o ~{prefix}.merged.vcf.gz \
            inputs/*.vcf.gz

        tabix -p vcf ~{prefix}.merged.vcf.gz

        bcftools annotate \
            -r ~{contig} \
            --set-id '%INFO/ME_TYPE\_%CHROM\_%POS\_%INFO/allele_length' \
            -Oz -o ~{prefix}.vcf.gz \
            ~{prefix}.merged.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{prefix}.vcf.gz"
        File merged_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcfs, "GB")) + 10,
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
