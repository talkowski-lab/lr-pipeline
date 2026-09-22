# Derived from broadinstitute/long-read-pipelines:
# https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/tasks/VariantCalling/Sniffles2.wdl

version 1.0

import "../utils/Structs.wdl"

workflow Sniffles {
    meta {
        description: [
            "This tool calls structural variants from a sample's aligned long reads using Sniffles2 (https://github.com/fritzsedlazeck/Sniffles). It emits both a bgzipped, indexed single-sample VCF and the sample's SNF file.",
            "The SNF file holds the sample's raw structural variant candidates and is what Sniffles2 population mode re-genotypes across a cohort, so it is retained even though this pipeline merges callsets by other means."
        ]
    }

    parameter_meta {
        bam: "Aligned reads for the sample."
        bai: "Index for the aligned reads."
        sample_id: "ID of the sample being called, written to the VCF sample column."
        min_sv_len: "Minimum structural variant length in base pairs to report."
        ref_fa: "From references."
        ref_fai: "From references."
        tandem_repeat_bed: "Tandem repeat intervals used to suppress alignment noise inside repeats."
        sniffles_vcf: "Structural variant calls for the sample."
        sniffles_vcf_idx: "Index for the structural variant calls."
        sniffles_snf: "Structural variant candidates for the sample, for later population-mode calling."
    }

    input {
        File bam
        File bai
        String prefix

        String sample_id
        Int min_sv_len = 50

        File ref_fa
        File ref_fai
        File? tandem_repeat_bed

        String sniffles_docker

        RuntimeAttr? runtime_attr_run_sniffles
    }

    call RunSniffles {
        input:
            bam = bam,
            bai = bai,
            sample_id = sample_id,
            min_sv_len = min_sv_len,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            tandem_repeat_bed = tandem_repeat_bed,
            prefix = prefix,
            docker = sniffles_docker,
            runtime_attr_override = runtime_attr_run_sniffles
    }

    output {
        File sniffles_vcf = RunSniffles.vcf
        File sniffles_vcf_idx = RunSniffles.vcf_idx
        File sniffles_snf = RunSniffles.snf
    }
}

task RunSniffles {
    input {
        File bam
        File bai
        String sample_id
        Int min_sv_len
        File ref_fa
        File ref_fai
        File? tandem_repeat_bed
        String prefix
        String docker

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # The reference is only needed so that deletion alleles carry their sequence rather than a symbolic ALT
        sniffles \
            --input ~{bam} \
            --reference ~{ref_fa} \
            ~{if defined(tandem_repeat_bed) then "--tandem-repeats " + tandem_repeat_bed else ""} \
            --sample-id ~{sample_id} \
            --minsvlen ~{min_sv_len} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --vcf ~{prefix}.sniffles.vcf.gz \
            --snf ~{prefix}.sniffles.snf

        tabix -f -p vcf ~{prefix}.sniffles.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.sniffles.vcf.gz"
        File vcf_idx = "~{prefix}.sniffles.vcf.gz.tbi"
        File snf = "~{prefix}.sniffles.snf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 32,
        disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB")) + 25,
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
