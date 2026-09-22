# Derived from broadinstitute/long-read-pipelines:
# https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/tasks/Alignment/AlignReads.wdl

version 1.0

import "../utils/Structs.wdl"

workflow MinimapReadAlignment {
    meta {
        description: [
            "This tool aligns a sample's unaligned long reads to a reference using Minimap2 (https://github.com/lh3/minimap2). Every unaligned BAM for the sample is converted to FASTQ, streamed through Minimap2 in a single pass and coordinate-sorted into one indexed BAM.",
            "Base modification tags are carried across from the unaligned BAM, since `samtools fastq` drops all tags by default and downstream methylation profiling needs them. Assemblies are aligned by `MinimapAlignment` instead."
        ]
    }

    parameter_meta {
        bams: "Unaligned BAMs for the sample, one per SMRT cell."
        sample_id: "ID of the sample being aligned, used for the read group ID and sample name."
        map_preset: "Minimap2 preset passed to '-x'."
        tags_to_preserve: "SAM tags carried over from the unaligned BAMs into the aligned BAM."
        ref_fa: "From references."
        ref_fai: "From references."
        aligned_bam: "Coordinate-sorted aligned reads."
        aligned_bai: "Index for the aligned reads."
    }

    input {
        Array[File] bams
        String prefix

        String sample_id
        String map_preset = "map-hifi"
        Array[String] tags_to_preserve = ["MM", "ML"]

        File ref_fa
        File ref_fai

        String minimap2_docker

        RuntimeAttr? runtime_attr_align_with_minimap2
    }

    call AlignWithMinimap2 {
        input:
            bams = bams,
            sample_id = sample_id,
            map_preset = map_preset,
            tags_to_preserve = tags_to_preserve,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = "~{prefix}.~{sample_id}",
            docker = minimap2_docker,
            runtime_attr_override = runtime_attr_align_with_minimap2
    }

    output {
        File aligned_bam = AlignWithMinimap2.aligned_bam
        File aligned_bai = AlignWithMinimap2.aligned_bai
    }
}

task AlignWithMinimap2 {
    input {
        Array[File] bams
        String sample_id
        String map_preset
        Array[String] tags_to_preserve
        File ref_fa
        File ref_fai
        String prefix
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Stream every unaligned BAM through Minimap2 in one pass, preserving the requested tags
        for bam in ~{sep=' ' bams}; do
            samtools fastq -@ 4 -T ~{sep=',' tags_to_preserve} "$bam"
        done \
            | minimap2 \
                -ayYL \
                --MD \
                --eqx \
                -x ~{map_preset} \
                -R "@RG\tID:~{sample_id}\tSM:~{sample_id}" \
                -t ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
                ~{ref_fa} \
                - \
            | samtools sort -@ 4 -m 2G --no-PG -o ~{prefix}.bam

        samtools index -@ 4 ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 32,
        mem_gb: 64,
        disk_gb: 4 * ceil(size(bams, "GB")) + ceil(size(ref_fa, "GB")) + 25,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
