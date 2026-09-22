# Derived from broadinstitute/long-read-pipelines:
# https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/tasks/VariantCalling/PBSV.wdl

version 1.0

import "../utils/Structs.wdl"

workflow PBSV {
    meta {
        description: [
            "This tool calls structural variants from a sample's aligned long reads using pbsv (https://github.com/PacificBiosciences/pbsv). Signatures of structural variation are discovered from the alignments and then genotyped into a bgzipped, indexed VCF.",
            "Supplying a tandem repeat BED lets pbsv collapse the alignment noise inside repeats, which reduces false calls at those loci."
        ]
    }

    parameter_meta {
        bam: "Aligned reads for the sample."
        bai: "Index for the aligned reads."
        is_hifi: "Whether the reads are HiFi, which enables the pbsv optimisations for low-error reads."
        ref_fa: "From references."
        ref_fai: "From references."
        tandem_repeat_bed: "Tandem repeat intervals used to suppress alignment noise inside repeats."
        pbsv_vcf: "Structural variant calls for the sample."
        pbsv_vcf_idx: "Index for the structural variant calls."
        pbsv_svsig: "Structural variant signatures discovered from the alignments."
    }

    input {
        File bam
        File bai
        String prefix

        Boolean is_hifi = true

        File ref_fa
        File ref_fai
        File? tandem_repeat_bed

        String pbsv_docker

        RuntimeAttr? runtime_attr_discover_signatures
        RuntimeAttr? runtime_attr_call_svs
    }

    call DiscoverSignatures {
        input:
            bam = bam,
            bai = bai,
            tandem_repeat_bed = tandem_repeat_bed,
            prefix = prefix,
            docker = pbsv_docker,
            runtime_attr_override = runtime_attr_discover_signatures
    }

    call CallSVs {
        input:
            svsig = DiscoverSignatures.svsig,
            is_hifi = is_hifi,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = pbsv_docker,
            runtime_attr_override = runtime_attr_call_svs
    }

    output {
        File pbsv_vcf = CallSVs.vcf
        File pbsv_vcf_idx = CallSVs.vcf_idx
        File pbsv_svsig = DiscoverSignatures.svsig
    }
}

task DiscoverSignatures {
    input {
        File bam
        File bai
        File? tandem_repeat_bed
        String prefix
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        pbsv discover \
            ~{if defined(tandem_repeat_bed) then "--tandem-repeats " + tandem_repeat_bed else ""} \
            ~{bam} \
            ~{prefix}.svsig.gz
    >>>

    output {
        File svsig = "~{prefix}.svsig.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: ceil(size(bam, "GB")) + 25,
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

task CallSVs {
    input {
        File svsig
        Boolean is_hifi
        File ref_fa
        File ref_fai
        String prefix
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        pbsv call \
            -j ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --log-level INFO \
            ~{true='--hifi' false='' is_hifi} \
            ~{ref_fa} \
            ~{svsig} \
            ~{prefix}.pbsv.vcf

        bgzip -@ 4 ~{prefix}.pbsv.vcf
        tabix -p vcf ~{prefix}.pbsv.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.pbsv.vcf.gz"
        File vcf_idx = "~{prefix}.pbsv.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 64,
        disk_gb: ceil(size(svsig, "GB") + size(ref_fa, "GB")) + 25,
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
