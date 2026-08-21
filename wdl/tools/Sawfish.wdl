version 1.0

import "../utils/Structs.wdl"

workflow Sawfish {
    input {
        File bam
        File bai
        String prefix

        String sex

        File ref_fa
        File ref_fai
        File expected_cn_male
        File expected_cn_female

        File? cnv_excluded_regions
        File? maf

        Int min_sv_size = 35
        Int min_sv_mapq = 5
        Boolean fast_cnv_mode = false
        Boolean disable_cnv = false
        Boolean treat_single_copy_as_haploid = false
        Boolean report_supporting_reads = false

        String sawfish_docker

        RuntimeAttr? runtime_attr_run_sawfish
    }

    File sex_specific_cn = if sex == 'M' then expected_cn_male else expected_cn_female

    call RunSawfish {
        input:
            bam = bam,
            bai = bai,
            prefix = prefix,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            expected_cn = sex_specific_cn,
            cnv_excluded_regions = cnv_excluded_regions,
            maf = maf,
            min_sv_size = min_sv_size,
            min_sv_mapq = min_sv_mapq,
            fast_cnv_mode = fast_cnv_mode,
            disable_cnv = disable_cnv,
            treat_single_copy_as_haploid = treat_single_copy_as_haploid,
            report_supporting_reads = report_supporting_reads,
            docker = sawfish_docker,
            runtime_attr_override = runtime_attr_run_sawfish
    }

    output {
        File sawfish_vcf = RunSawfish.sawfish_vcf
        File sawfish_vcf_idx = RunSawfish.sawfish_vcf_idx
        File sawfish_bedgraph = RunSawfish.sawfish_bedgraph
        File sawfish_depth_bw = RunSawfish.sawfish_depth_bw
        File sawfish_log = RunSawfish.sawfish_log
        File? sawfish_supporting_reads = RunSawfish.sawfish_supporting_reads
    }
}

task RunSawfish {
    input {
        File bam
        File bai
        String prefix
        File ref_fa
        File ref_fai
        File expected_cn
        File? cnv_excluded_regions
        File? maf
        Int min_sv_size
        Int min_sv_mapq
        Boolean fast_cnv_mode
        Boolean disable_cnv
        Boolean treat_single_copy_as_haploid
        Boolean report_supporting_reads
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Discover SV candidates and build the depth/CNV segmentation for the sample
        sawfish discover \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --ref ~{ref_fa} \
            --bam ~{bam} \
            --expected-cn ~{expected_cn} \
            ~{if defined(cnv_excluded_regions) then "--cnv-excluded-regions " + cnv_excluded_regions else ""} \
            ~{if defined(maf) then "--maf " + maf else ""} \
            --min-indel-size ~{min_sv_size} \
            --min-sv-mapq ~{min_sv_mapq} \
            ~{if fast_cnv_mode then "--fast-cnv-mode" else ""} \
            ~{if disable_cnv then "--disable-cnv" else ""} \
            --output-dir discover_dir

        # Merge and genotype SVs/CNVs to produce the final VCF
        sawfish joint-call \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --sample discover_dir \
            --min-sv-mapq ~{min_sv_mapq} \
            ~{if treat_single_copy_as_haploid then "--treat-single-copy-as-haploid" else ""} \
            ~{if report_supporting_reads then "--report-supporting-reads" else ""} \
            --output-dir joint_call_dir

        mv joint_call_dir/genotyped.sv.vcf.gz ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz

        cp joint_call_dir/samples/sample*/copynum.bedgraph ~{prefix}.copynum.bedgraph
        cp joint_call_dir/samples/sample*/depth.bw ~{prefix}.depth.bw
        cp joint_call_dir/sawfish.log ~{prefix}.sawfish.log

        if [ -f joint_call_dir/supporting_reads.json.gz ]; then
            cp joint_call_dir/supporting_reads.json.gz ~{prefix}.supporting_reads.json.gz
        fi
    >>>

    output {
        File sawfish_vcf = "~{prefix}.vcf.gz"
        File sawfish_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File sawfish_bedgraph = "~{prefix}.copynum.bedgraph"
        File sawfish_depth_bw = "~{prefix}.depth.bw"
        File sawfish_log = "~{prefix}.sawfish.log"
        File? sawfish_supporting_reads = "~{prefix}.supporting_reads.json.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 16,
        mem_gb: 64,
        disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB")) + 30,
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
