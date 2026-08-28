version 1.0

import "utils/Structs.wdl"

workflow ExtractBamRegion {
    input {
        File bam
        File bai
        Int start
        Int end
        String chrom
        String gatk_docker
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_extract_region
    }

    call ExtractRegion {
        input:
            bam = bam,
            bai = bai,
            chrom = chrom,
            start = start,
            end = end,
            docker = gatk_docker,
            runtime_attr_override = runtime_attr_extract_region
    }

    output {
        File regional_bam = ExtractRegion.regional_bam
        File regional_bai = ExtractRegion.regional_bai
    }
}

task ExtractRegion {
    input {
        File bam
        File bai
        String chrom
        Int start
        Int end
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam")

    command <<<
        set -euo pipefail

        gatk PrintReads \
           -I ~{bam} \
           -L ~{chrom}:~{start}-~{end} \
           -O "~{prefix}.~{chrom}_~{start}_~{end}.bam"

        samtools index "~{prefix}.~{chrom}_~{start}_~{end}.bam"
    >>>

    output {
        File regional_bam = "~{prefix}.~{chrom}_~{start}_~{end}.bam"
        File regional_bai = "~{prefix}.~{chrom}_~{start}_~{end}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB")) + 25,
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
