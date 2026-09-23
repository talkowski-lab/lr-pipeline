version 1.0

import "utils/Structs.wdl"

workflow PALMERToVcf {
    meta {
        description: [
            "This utility converts a sample's PALMER mobile-element calls into a VCF. Each mobile-element type is converted separately, and the per-type records are concatenated and sorted into one indexed VCF."
        ]
    }

    parameter_meta {
        PALMER_calls: "PALMER call files, one per entry in `mei_types`."
        mei_types: "Mobile-element type for each entry in `PALMER_calls`, in the same order."
        sample: "ID of the sample being processed."
        ref_fai: "From references."
        PALMER_combined_vcf: "PALMER calls converted to VCF."
        PALMER_combined_vcf_idx: "Index for `PALMER_combined_vcf`."
    }

    input {
        Array[File] PALMER_calls
        Array[String] mei_types

        String prefix
        String sample

        File ref_fai

        String pipeline_docker
        
        RuntimeAttr? runtime_attr_palmer_to_vcf
        RuntimeAttr? runtime_attr_concat_sort_vcfs
    }

    scatter(i in range(length(PALMER_calls))) {
        call ConvertPALMERToVcf {
            input:
                PALMER_calls = PALMER_calls[i],
                mei_type = mei_types[i],
                ref_fai = ref_fai,
                sample = sample,
                docker = pipeline_docker,
                runtime_attr_override = runtime_attr_palmer_to_vcf
        }
    }

    call ConcatSortVcfs {
        input:
            vcfs = ConvertPALMERToVcf.vcf,
            vcf_idxs = ConvertPALMERToVcf.vcf_idx,
            docker = pipeline_docker,
            prefix = prefix,
            runtime_attr_override = runtime_attr_concat_sort_vcfs
    }

    output {
        File PALMER_combined_vcf = ConcatSortVcfs.vcf
        File PALMER_combined_vcf_idx = ConcatSortVcfs.vcf_idx
    }
}

task ConvertPALMERToVcf {
    input {
        File PALMER_calls
        String mei_type
        String sample
        File ref_fai
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python /opt/gnomad-lr/scripts/mei/PALMER_to_vcf.py \
            --palmer_calls ~{PALMER_calls} \
            --mei_type ~{mei_type} \
            --sample ~{sample} \
            --ref_fai ~{ref_fai} \
            | bcftools sort -Oz \
            > ~{sample}.PALMER_calls.~{mei_type}.vcf.gz
        
        tabix ~{sample}.PALMER_calls.~{mei_type}.vcf.gz
    >>>

    output {
        File vcf = "~{sample}.PALMER_calls.~{mei_type}.vcf.gz"
        File vcf_idx = "~{sample}.PALMER_calls.~{mei_type}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 10*ceil(size(PALMER_calls, "GB") + size(ref_fai, "GB")) + 20,
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

task ConcatSortVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat -a ~{sep=" " vcfs} | bcftools sort -Oz > ~{prefix}.vcf.gz
        tabix ~{prefix}.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 5*ceil(size(vcfs, "GB")) + 5,
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
