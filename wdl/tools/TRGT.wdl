version 1.0

import "../utils/Structs.wdl"

workflow TRGT {
    meta {
        description: [
            "This workflow leverages TRGT (https://github.com/PacificBiosciences/trgt) in order to genotype short-tandem repeats."
        ]
    }

    parameter_meta {
        bam: "Aligned reads."
        bai: "Index for aligned reads."
        sample_id: "ID of the sample being genotyped."
        sex: "Sex of sample (one of `M` or `F`)."
        catalog_name: "Name of the repeat catalog used, included in the output VCF filename."
        ref_fa: "From references."
        ref_fai: "From references."
        repeat_catalog_trgt: "From references."
        trgt_vcf: "TRGT tandem-repeat genotype VCF."
        trgt_vcf_idx: "Index for the TRGT VCF."
    }

    input {
        File bam
        File bai
        String prefix

        String sample_id
        String sex
        String catalog_name

        File ref_fa
        File ref_fai
        File repeat_catalog_trgt

        String trgt_docker

        RuntimeAttr? runtime_attr_process_with_trgt
    }

    Boolean is_female = 'F' == sex

    call ProcessWithTRGT {
        input:
            bam = bam,
            bai = bai,
            is_female = is_female,
            catalog_name = catalog_name,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            repeat_catalog_trgt = repeat_catalog_trgt,
            prefix = "~{prefix}.~{sample_id}",
            docker = trgt_docker,
            runtime_attr_override = runtime_attr_process_with_trgt
    }

    output {
        File trgt_vcf = ProcessWithTRGT.trgt_output_vcf
        File trgt_vcf_idx = ProcessWithTRGT.trgt_output_vcf_idx
    }
}

task ProcessWithTRGT {
    input {
        File bam
        File bai
        Boolean is_female
        Boolean verbose = false
        String catalog_name
        File ref_fa
        File ref_fai
        File repeat_catalog_trgt
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String vcf_out_name = prefix + "_trgt." + catalog_name
    String karyotype = if(is_female) then "XX" else "XY"

    command <<<
        set -euo pipefail

        time \
        trgt \
            ~{true='--verbose' false=' ' verbose} \
            genotype \
            --genome ~{ref_fa} \
            --repeats ~{repeat_catalog_trgt} \
            --reads ~{bam} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --output-prefix ~{vcf_out_name} \
            --karyotype ~{karyotype}

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Ob -o ~{vcf_out_name}.sorted.vcf.gz \
            ~{vcf_out_name}.vcf.gz

        mv ~{vcf_out_name}.sorted.vcf.gz ~{vcf_out_name}.vcf.gz

        tabix -p vcf ~{vcf_out_name}.vcf.gz
    >>>

    output {
        File trgt_output_vcf = "~{vcf_out_name}.vcf.gz"
        File trgt_output_vcf_idx = "~{vcf_out_name}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 6,
        disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB") + size(repeat_catalog_trgt, "GB")) + 25,
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
