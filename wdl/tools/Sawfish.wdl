version 1.0

import "../utils/Structs.wdl"

workflow Sawfish {
    input {
        Array[File] bams
        Array[File] bais
        Array[String] sexes
        Array[String] sample_ids
        String prefix

        File ref_fa
        File ref_fai
        File expected_cn_male
        File expected_cn_female

        File exclude_bed
        File exclude_bed_idx

        Int min_sv_size = 35
        Int min_sv_mapq = 5
        Boolean fast_cnv_mode = false
        Boolean disable_cnv = false
        Boolean treat_single_copy_as_haploid = false
        Boolean report_supporting_reads = false

        String sawfish_docker

        RuntimeAttr? runtime_attr_discover
        RuntimeAttr? runtime_attr_joint_call
    }

    scatter (i in range(length(bams))) {
        File expected_cn = if sexes[i] == 'M' then expected_cn_male else expected_cn_female

        call Discover {
            input:
                bam = bams[i],
                bai = bais[i],
                sample_id = sample_ids[i],
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                expected_cn = expected_cn,
                exclude_bed = exclude_bed,
                exclude_bed_idx = exclude_bed_idx,
                min_sv_size = min_sv_size,
                min_sv_mapq = min_sv_mapq,
                fast_cnv_mode = fast_cnv_mode,
                disable_cnv = disable_cnv,
                docker = sawfish_docker,
                runtime_attr_override = runtime_attr_discover
        }
    }

    call JointCall {
        input:
            discover_tars = Discover.discover_tar,
            bams = bams,
            bais = bais,
            sample_ids = sample_ids,
            prefix = prefix,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            min_sv_mapq = min_sv_mapq,
            treat_single_copy_as_haploid = treat_single_copy_as_haploid,
            report_supporting_reads = report_supporting_reads,
            docker = sawfish_docker,
            runtime_attr_override = runtime_attr_joint_call
    }

    output {
        File sawfish_vcf = JointCall.sawfish_vcf
        File sawfish_vcf_idx = JointCall.sawfish_vcf_idx
        Array[File] sawfish_bedgraphs = JointCall.sawfish_bedgraphs
        Array[File] sawfish_depth_bws = JointCall.sawfish_depth_bws
        File sawfish_log = JointCall.sawfish_log
        File? sawfish_supporting_reads = JointCall.sawfish_supporting_reads
    }
}

task Discover {
    input {
        File bam
        File bai
        String sample_id
        File ref_fa
        File ref_fai
        File expected_cn
        File exclude_bed
        File exclude_bed_idx
        Int min_sv_size
        Int min_sv_mapq
        Boolean fast_cnv_mode
        Boolean disable_cnv
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        sawfish discover \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --ref ~{ref_fa} \
            --bam ~{bam} \
            --expected-cn ~{expected_cn} \
            ~{if defined(exclude_bed) then "--cnv-excluded-regions " + exclude_bed else ""} \
            --min-indel-size ~{min_sv_size} \
            --min-sv-mapq ~{min_sv_mapq} \
            ~{if fast_cnv_mode then "--fast-cnv-mode" else ""} \
            ~{if disable_cnv then "--disable-cnv" else ""} \
            --output-dir ~{sample_id}_discover

        tar czf ~{sample_id}_discover.tar.gz ~{sample_id}_discover
    >>>

    output {
        File discover_tar = "~{sample_id}_discover.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 4,
        disk_gb: ceil(size(bam, "GB") + size(ref_fa, "GB")) + 30,
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

task JointCall {
    input {
        Array[File] discover_tars
        Array[File] bams
        Array[File] bais
        Array[String] sample_ids
        String prefix
        File ref_fa
        File ref_fai
        Int min_sv_mapq
        Boolean treat_single_copy_as_haploid
        Boolean report_supporting_reads
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        tars=(~{sep=" " discover_tars})
        bams=(~{sep=" " bams})
        bais=(~{sep=" " bais})
        ids=(~{sep=" " sample_ids})

        mkdir -p aligned
        : > samples.csv
        for i in "${!tars[@]}"; do
            tar xzf "${tars[$i]}"
            ln -s "${bams[$i]}" "aligned/${ids[$i]}.bam"
            ln -s "${bais[$i]}" "aligned/${ids[$i]}.bam.bai"
            echo "${ids[$i]}_discover, ${PWD}/aligned/${ids[$i]}.bam" >> samples.csv
        done
        
        sawfish joint-call \
            --threads ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])} \
            --ref ~{ref_fa} \
            --sample-csv samples.csv \
            --min-sv-mapq ~{min_sv_mapq} \
            ~{if treat_single_copy_as_haploid then "--treat-single-copy-as-haploid" else ""} \
            ~{if report_supporting_reads then "--report-supporting-reads" else ""} \
            --output-dir joint_call_dir

        mv joint_call_dir/genotyped.sv.vcf.gz ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz

        cp joint_call_dir/sawfish.log ~{prefix}.sawfish.log

        if [ -f joint_call_dir/supporting_reads.json.gz ]; then
            cp joint_call_dir/supporting_reads.json.gz ~{prefix}.supporting_reads.json.gz
        fi
    >>>

    output {
        File sawfish_vcf = "~{prefix}.vcf.gz"
        File sawfish_vcf_idx = "~{prefix}.vcf.gz.tbi"
        Array[File] sawfish_bedgraphs = glob("joint_call_dir/samples/sample*/copynum.bedgraph")
        Array[File] sawfish_depth_bws = glob("joint_call_dir/samples/sample*/depth.bw")
        File sawfish_log = "~{prefix}.sawfish.log"
        File? sawfish_supporting_reads = "~{prefix}.supporting_reads.json.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 4,
        disk_gb: ceil(size(bams, "GB") + size(discover_tars, "GB") + size(ref_fa, "GB")) + 30,
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
