version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow CountVariantClasses {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        Int max_indel = 49

        File count_script

        String utils_docker

        RuntimeAttr? runtime_attr_count
        RuntimeAttr? runtime_attr_sum
    }

    scatter (i in range(length(vcfs))) {
        call CountPerVcf {
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxs[i],
                max_indel = max_indel,
                count_script = count_script,
                prefix = "~{prefix}.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_count
        }
    }

    call SumCounts {
        input:
            counts = CountPerVcf.per_contig,
            count_script = count_script,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_sum
    }

    output {
        File per_contig_counts = SumCounts.per_contig
        File total_counts = SumCounts.total
    }
}

task CountPerVcf {
    input {
        File vcf
        File vcf_idx
        Int max_indel
        File count_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{count_script} \
            --vcfs ~{vcf} \
            --max-indel ~{max_indel} \
            --out-per-contig ~{prefix}.per_contig.tsv \
            --out-total ~{prefix}.total.tsv
    >>>

    output {
        File per_contig = "~{prefix}.per_contig.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) + 10,
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

task SumCounts {
    input {
        Array[File] counts
        File count_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 ~{count_script} \
            --counts ~{sep=" " counts} \
            --out-per-contig ~{prefix}.variant_classes.per_contig.tsv \
            --out-total ~{prefix}.variant_classes.total.tsv
    >>>

    output {
        File per_contig = "~{prefix}.variant_classes.per_contig.tsv"
        File total = "~{prefix}.variant_classes.total.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10,
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
