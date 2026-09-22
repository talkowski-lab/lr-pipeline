version 1.0

import "../utils/Structs.wdl"

workflow HiFiCNV {
    meta {
        description: [
            "This tool runs PacBio HiFiCNV (https://github.com/PacificBiosciences/HiFiCNV) on a sample's aligned HiFi BAM to call copy number variants from read depth. It outputs the CNV VCF, a copy-number bedgraph, a depth BigWig track and the tool's log."
        ]
    }

    parameter_meta {
        bam: "Aligned reads for the sample."
        bai: "Index for `bam`."
        sex: "Sex of sample (one of `M` or `F`), used to select the matching expected-CN file."
        ref_fa: "Reference sequences FASTA file."
        ref_fai: "Index for `ref_fa`."
        exclude_bed: "Regions to exclude from CNV calling (e.g. centromeres)."
        exclude_bed_idx: "Index for `exclude_bed`."
        expected_cn_male: "PAR regions and expected copy numbers for sex chromosomes, male."
        expected_cn_female: "PAR regions and expected copy numbers for sex chromosomes, female."
        maf: "Optional minor-allele-frequency track passed to HiFiCNV as `--maf`."
        cov_regex: "Optional regular expression passed to HiFiCNV as `--cov-regex`, selecting the contigs used to estimate expected coverage."
        disable_vcf_filters: "Whether to pass `--disable-vcf-filters`, emitting every call rather than only those HiFiCNV would keep."
        hificnv_vcf: "CNV calls VCF."
        hificnv_vcf_idx: "Index for the CNV calls VCF."
        hificnv_bedgraph: "Per-window copy number bedgraph."
        hificnv_depth_bw: "Depth BigWig track."
        hificnv_log: "HiFiCNV log file."
    }

    input {
        File bam
        File bai
        String prefix

        String sex

        File ref_fa
        File ref_fai
        File exclude_bed
        File exclude_bed_idx
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
            exclude_bed_idx = exclude_bed_idx,
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
        File exclude_bed_idx
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
