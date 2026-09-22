version 1.0

import "../utils/Structs.wdl"

workflow ConcatenateMosDepth {
    meta {
        description: [
            "This utility concatenates a sample's per-contig `MosDepth` per-base coverage BED files into a single indexed BED. It outputs the combined per-base coverage BED and its index."
        ]
    }

    parameter_meta {
        mosdepth_bed_files: "Per-contig mosdepth per-base coverage BED files to concatenate."
        mosdepth_per_base_combined: "Combined per-base coverage BED."
        mosdepth_per_base_combined_idx: "Index for the combined BED."
    }

    input {
        Array[File] mosdepth_bed_files
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_concatenate
    }

    call ConcatenateBeds {
        input:
            bed_files = mosdepth_bed_files,
            prefix = "~{prefix}.combined",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concatenate
    }

    output {
        File mosdepth_per_base_combined = ConcatenateBeds.concatenated_bed
        File mosdepth_per_base_combined_idx = ConcatenateBeds.concatenated_bed_idx
    }
}

task ConcatenateBeds {
    input {
        Array[File] bed_files
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat ~{sep=' ' bed_files} | bgzip > ~{prefix}.bed.gz

        tabix -s 1 -b 2 -e 3 ~{prefix}.bed.gz
    >>>

    output {
        File concatenated_bed = "~{prefix}.bed.gz"
        File concatenated_bed_idx = "~{prefix}.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(bed_files, "GB")) + 10,
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
