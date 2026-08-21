version 1.0

import "../utils/Structs.wdl"

workflow HiFiCNV {
    input {
        File bam
        File bai
        String prefix

        String sex

        File ref_fa
        File ref_fai
        File exclude_bed
        File expected_cn_male
        File expected_cn_female

        File? maf
        String? cov_regex
        Boolean disable_vcf_filters = false

        String hificnv_docker

        RuntimeAttr? runtime_attr_run_hificnv
    }

    File sex_specific_cn = if sex == 'M' then expected_cn_male else expected_cn_female

    call RunHiFiCNV {
        input:
            bam = bam,
            bai = bai,
            prefix = prefix,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            exclude_bed = exclude_bed,
            sex_specific_cn = sex_specific_cn,
            maf = maf,
            cov_regex = cov_regex,
            disable_vcf_filters = disable_vcf_filters,
            docker = hificnv_docker,
            runtime_attr_override = runtime_attr_run_hificnv
    }

    output {
        File hificnv_vcf = RunHiFiCNV.hificnv_vcf
        File hificnv_vcf_idx = RunHiFiCNV.hificnv_vcf_idx
        File hificnv_bedgraph = RunHiFiCNV.hificnv_bedgraph
        File hificnv_depth_bw = RunHiFiCNV.hificnv_depth_bw
        File hificnv_log = RunHiFiCNV.hificnv_log
    }
}

task RunHiFiCNV {
    input {
        File bam
        File bai
        String prefix
        File ref_fa
        File ref_fai
        File exclude_bed
        File sex_specific_cn
        File? maf
        String? cov_regex
        Boolean disable_vcf_filters
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        hificnv \
            --bam ~{bam} \
            --ref ~{ref_fa} \
            --exclude ~{exclude_bed} \
            --expected-cn ~{sex_specific_cn} \
            ~{if defined(maf) then "--maf " + maf else ""} \
            ~{if defined(cov_regex) then "--cov-regex " + cov_regex else ""} \
            ~{if disable_vcf_filters then "--disable-vcf-filters" else ""} \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --output-prefix ~{prefix}

        mv ~{prefix}.*.vcf.gz ~{prefix}.vcf.gz
        mv ~{prefix}.*.copynum.bedgraph ~{prefix}.copynum.bedgraph
        mv ~{prefix}.*.depth.bw ~{prefix}.depth.bw

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File hificnv_vcf = "~{prefix}.vcf.gz"
        File hificnv_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File hificnv_bedgraph = "~{prefix}.copynum.bedgraph"
        File hificnv_depth_bw = "~{prefix}.depth.bw"
        File hificnv_log = "~{prefix}.log"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 6,
        disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB")) + 20,
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
