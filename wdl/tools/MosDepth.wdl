version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow MosDepth {
    input {
        File bam
        File bai
        Array[String] contigs
        String prefix

        Boolean single_contig
        Boolean stream_mode
        Boolean fast_mode

        Int? bin_size

        File? ref_fa
        File? ref_fai

        String mosdepth_docker
        String mosdepthstream_docker
        String utils_docker

        RuntimeAttr? runtime_attr_run_mosdepth
        RuntimeAttr? runtime_attr_run_mosdepth_stream
        RuntimeAttr? runtime_attr_split_bam
    }

    if (single_contig) {
        call RunMosDepth as RunMosDepthAllContigs {
            input:
                bam = bam,
                bai = bai,
                bin_size = bin_size,
                fast_mode = fast_mode,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.coverage",
                docker = mosdepth_docker,
                runtime_attr_override = runtime_attr_run_mosdepth
        }
    }

    if (!single_contig && stream_mode) {
        scatter (contig in contigs) {
            call RunMosDepthStream as RunMosDepthPerContigStream {
                input:
                    bam = bam,
                    bai = bai,
                    contig = contig,
                    bin_size = bin_size,
                    fast_mode = fast_mode,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}.~{contig}.coverage",
                    docker = mosdepthstream_docker,
                    runtime_attr_override = runtime_attr_run_mosdepth_stream
            }
        }
    }

    if (!single_contig && !stream_mode) {
        call Helpers.SplitBamByContig {
            input:
                bam = bam,
                bai = bai,
                contigs = contigs,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_split_bam
        }

        scatter (i in range(length(contigs))) {
            call RunMosDepth as RunMosDepthPerContig {
                input:
                    bam = SplitBamByContig.contig_bams[i],
                    bai = SplitBamByContig.contig_bais[i],
                    contig = contigs[i],
                    bin_size = bin_size,
                    fast_mode = fast_mode,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}.~{contigs[i]}.coverage",
                    docker = mosdepth_docker,
                    runtime_attr_override = runtime_attr_run_mosdepth
            }
        }
    }

    output {
        Array[File] mosdepth_dist = flatten([select_all([RunMosDepthAllContigs.dist]), flatten(select_all([RunMosDepthPerContigStream.dist])), flatten(select_all([RunMosDepthPerContig.dist]))])
        Array[File] mosdepth_summary = flatten([select_all([RunMosDepthAllContigs.summary]), flatten(select_all([RunMosDepthPerContigStream.summary])), flatten(select_all([RunMosDepthPerContig.summary]))])
        Array[File] mosdepth_per_base = flatten([select_all([RunMosDepthAllContigs.per_base]), select_all(flatten(select_all([RunMosDepthPerContigStream.per_base]))), select_all(flatten(select_all([RunMosDepthPerContig.per_base])))])
        Array[File] mosdepth_per_base_csi = flatten([select_all([RunMosDepthAllContigs.per_base_csi]), select_all(flatten(select_all([RunMosDepthPerContigStream.per_base_csi]))), select_all(flatten(select_all([RunMosDepthPerContig.per_base_csi])))])
        Array[File] mosdepth_regions_bed = flatten([select_all([RunMosDepthAllContigs.regions_bed]), select_all(flatten(select_all([RunMosDepthPerContigStream.regions_bed]))), select_all(flatten(select_all([RunMosDepthPerContig.regions_bed])))])
        Array[File] mosdepth_regions_bed_csi = flatten([select_all([RunMosDepthAllContigs.regions_bed_csi]), select_all(flatten(select_all([RunMosDepthPerContigStream.regions_bed_csi]))), select_all(flatten(select_all([RunMosDepthPerContig.regions_bed_csi])))])
    }
}

task RunMosDepth {
    input {
        File bam
        File bai
        String? contig
        Int? bin_size
        Boolean fast_mode
        File? ref_fa
        File? ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mosdepth \
            ~{if defined(contig) then "-c " + contig else ""} \
            ~{if defined(bin_size) then "--by " + bin_size + " --no-per-base" else ""} \
            ~{if defined(ref_fa) then "--fasta " + ref_fa else ""} \
            ~{if fast_mode then "--fast-mode" else ""} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix} \
            ~{bam}
    >>>

    output {
        File dist = "~{prefix}.mosdepth.global.dist.txt"
        File summary = "~{prefix}.mosdepth.summary.txt"
        File? per_base = "~{prefix}.per-base.bed.gz"
        File? per_base_csi = "~{prefix}.per-base.bed.gz.csi"
        File? regions_bed = "~{prefix}.regions.bed.gz"
        File? regions_bed_csi = "~{prefix}.regions.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 0,
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

task RunMosDepthStream {
    input {
        File bam
        File bai
        String contig
        Int? bin_size
        Boolean fast_mode
        File? ref_fa
        File? ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: { localization_optional: true }
        bai: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        mosdepth \
            -c ~{contig} \
            ~{if defined(bin_size) then "--by " + bin_size + " --no-per-base" else ""} \
            ~{if defined(ref_fa) then "--fasta " + ref_fa else ""} \
            ~{if fast_mode then "--fast-mode" else ""} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            ~{prefix} \
            ~{bam}
    >>>

    output {
        File dist = "~{prefix}.mosdepth.global.dist.txt"
        File summary = "~{prefix}.mosdepth.summary.txt"
        File? per_base = "~{prefix}.per-base.bed.gz"
        File? per_base_csi = "~{prefix}.per-base.bed.gz.csi"
        File? regions_bed = "~{prefix}.regions.bed.gz"
        File? regions_bed_csi = "~{prefix}.regions.bed.gz.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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
