version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AlignQueriesBWA {
    input {
        File queries_fa
        String prefix

        String bwa_args = "-a -Y"

        File ref_fa
        Array[File] ref_bwa_idx

        String bwa_docker

        RuntimeAttr? runtime_attr_bwa_mem
    }

    call BwaMemQueries {
        input:
            queries_fa = queries_fa,
            bwa_args = bwa_args,
            ref_fa = ref_fa,
            ref_bwa_idx = ref_bwa_idx,
            prefix = prefix,
            docker = bwa_docker,
            runtime_attr_override = runtime_attr_bwa_mem
    }

    output {
        File alignments_sam = BwaMemQueries.alignments_sam
    }
}

task BwaMemQueries {
    input {
        File queries_fa
        String bwa_args
        File ref_fa
        Array[File] ref_bwa_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir ref_files
        mv ~{ref_fa} ref_files/
        for x in ~{sep=' ' ref_bwa_idx}; do
            mv "$x" ref_files/
        done

        bwa mem -t "$(nproc)" ~{bwa_args} ref_files/$(basename ~{ref_fa}) ~{queries_fa} \
            | samtools view -h -F 4 - \
            | gzip > ~{prefix}.bwa.sam.gz
    >>>

    output {
        File alignments_sam = "~{prefix}.bwa.sam.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(ref_fa, "GB") + size(ref_bwa_idx, "GB") + size(queries_fa, "GB")) + 10,
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
