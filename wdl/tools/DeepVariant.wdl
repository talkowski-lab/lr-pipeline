# Derived from broadinstitute/long-read-pipelines:
# https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/tasks/VariantCalling/DeepVariant.wdl [source branch: sh_update_wgs_callers_outliers]

version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow DeepVariant {
    meta {
        description: "Call and merge small variants with DeepVariant from region-sharded BAMs."
    }

    parameter_meta {
        how_to_shard_wg_for_calling: "Shard ID paired with its BAM and BAI."
        shard_region_files: "Map from each non-alts shard ID to a file containing one DeepVariant region literal per line. Regions must be disjoint across shards."
        ref_fa: "Reference FASTA file."
        ref_fai: "FASTA index for ref_fa."
        prefix: "Prefix for merged VCF, gVCF, and monitoring outputs."
        model_type: "DeepVariant model type, for example PACBIO or ONT_R104."
        haploid_contigs: "Optional comma-separated haploid contigs."
        par_regions_bed: "Optional pseudoautosomal-region BED file."
        threads: "Requested CPU count and DeepVariant internal shard count."
        memory: "Requested task memory in GiB."
        use_gpu: "Run GPU DeepVariant tasks instead of CPU tasks."
        zones: "Space-separated Google Cloud zones for task placement."
        deepvariant_docker: "Docker image for CPU DeepVariant tasks."
        deepvariant_gpu_docker: "Docker image for GPU DeepVariant tasks."
        utils_docker: "Docker image containing bcftools and tabix for VCF merging."
        resource_visualization_docker: "Docker image containing /opt/plot.resources.R."
        runtime_attr_run_deepvariant: "Override runtime attributes for CPU DeepVariant tasks."
        runtime_attr_run_deepvariant_gpu: "Override runtime attributes for GPU DeepVariant tasks."
        runtime_attr_visualize_resource_usage: "Override runtime attributes for resource-usage plotting tasks."
        runtime_attr_merge_gvcfs: "Override runtime attributes for gVCF merge task."
        runtime_attr_merge_vcfs: "Override runtime attributes for VCF merge task."
    }

    input {
        Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling
        Map[String, File] shard_region_files
        File ref_fa
        File ref_fai
        String prefix
        String model_type
        String? haploid_contigs
        File? par_regions_bed
        Int threads
        Int memory
        Boolean use_gpu = false
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        String deepvariant_docker
        String deepvariant_gpu_docker
        String utils_docker
        String resource_visualization_docker

        RuntimeAttr? runtime_attr_run_deepvariant
        RuntimeAttr? runtime_attr_run_deepvariant_gpu
        RuntimeAttr? runtime_attr_visualize_resource_usage
        RuntimeAttr? runtime_attr_merge_gvcfs
        RuntimeAttr? runtime_attr_merge_vcfs
    }

    scatter (shard in how_to_shard_wg_for_calling) {
        String shard_id = shard.left
        File shard_bam = shard.right.left
        File shard_bai = shard.right.right

        if (shard_id != "alts") {
            File shard_region_file = shard_region_files[shard_id]
            Array[String] shard_regions = read_lines(shard_region_file)

            Boolean is_t2t_chr_x = length(how_to_shard_wg_for_calling) < 20 && shard_id == "18_X"
            Boolean is_t2t_shard_8 = length(how_to_shard_wg_for_calling) < 20 && shard_id == "9_15"
            Boolean is_grch38_chr_1 = length(how_to_shard_wg_for_calling) > 20 && shard_id == "1-p"
            Boolean is_grch38_shard_3 = length(how_to_shard_wg_for_calling) > 20 && shard_id == "11-q_17-q"
            Boolean needs_extra_memory = is_t2t_chr_x || is_t2t_shard_8 || is_grch38_chr_1 || is_grch38_shard_3
            Int cpu_memory = if (needs_extra_memory) then 48 else memory

            if (!use_gpu) {
                call RunDeepVariant as RunCpu {
                    input:
                        bam = shard_bam,
                        bai = shard_bai,
                        ref_fa = ref_fa,
                        ref_fai = ref_fai,
                        regions = shard_regions,
                        haploid_contigs = haploid_contigs,
                        par_regions_bed = par_regions_bed,
                        model_type = model_type,
                        threads = threads,
                        memory = cpu_memory,
                        zones = zones,
                        prefix = "~{prefix}.~{shard_id}",
                        docker = deepvariant_docker,
                        runtime_attr_override = runtime_attr_run_deepvariant
                }
            }

            if (use_gpu) {
                call RunDeepVariantGpu as RunGpu {
                    input:
                        bam = shard_bam,
                        bai = shard_bai,
                        ref_fa = ref_fa,
                        ref_fai = ref_fai,
                        regions = shard_regions,
                        haploid_contigs = haploid_contigs,
                        par_regions_bed = par_regions_bed,
                        model_type = model_type,
                        threads = threads,
                        memory = memory,
                        zones = zones,
                        prefix = "~{prefix}.~{shard_id}",
                        docker = deepvariant_gpu_docker,
                        runtime_attr_override = runtime_attr_run_deepvariant_gpu
                }
            }

            File shard_vcf = select_first([RunCpu.vcf, RunGpu.vcf])
            File shard_vcf_idx = select_first([RunCpu.vcf_idx, RunGpu.vcf_idx])
            File shard_gvcf = select_first([RunCpu.gvcf, RunGpu.gvcf])
            File shard_gvcf_idx = select_first([RunCpu.gvcf_idx, RunGpu.gvcf_idx])
            File resource_usage_log = select_first([RunCpu.resource_usage_log, RunGpu.resource_usage_log])
            File visual_report = select_first([RunCpu.visual_report, RunGpu.visual_report])

            call VisualizeResourceUsage {
                input:
                    resource_log = resource_usage_log,
                    output_pdf_name = "~{prefix}.deepvariant.resource-usage.~{shard_id}.pdf",
                    plot_title = "DeepVariant on ~{prefix}, shard ~{shard_id}",
                    docker = resource_visualization_docker,
                    runtime_attr_override = runtime_attr_visualize_resource_usage
            }
        }
    }

    call Helpers.ConcatVcfs as ConcatGvcfs {
        input:
            vcfs = select_all(shard_gvcf),
            vcf_idxs = select_all(shard_gvcf_idx),
            allow_overlaps = true,
            naive = false,
            sort_output = true,
            no_version = false,
            no_address = false,
            prefix = "~{prefix}.deepvariant.g",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_gvcfs
    }

    call Helpers.ConcatVcfs as ConcatVcfs {
        input:
            vcfs = select_all(shard_vcf),
            vcf_idxs = select_all(shard_vcf_idx),
            allow_overlaps = true,
            naive = false,
            sort_output = true,
            no_version = false,
            no_address = false,
            prefix = "~{prefix}.deepvariant",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    output {
        File gvcf = ConcatGvcfs.concat_vcf
        File gvcf_idx = ConcatGvcfs.concat_vcf_idx
        File vcf = ConcatVcfs.concat_vcf
        File vcf_idx = ConcatVcfs.concat_vcf_idx
        Array[File] resource_usage_logs = select_all(resource_usage_log)
        Array[File] resource_usage_visualizations = select_all(VisualizeResourceUsage.plot_pdf)
        Array[File] visual_reports = select_all(visual_report)
    }
}

task RunDeepVariant {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        Array[String] regions
        String? haploid_contigs
        File? par_regions_bed
        String model_type
        Int threads
        Int memory
        String zones
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String output_root = "/mnt/disks/cromwell_root/dv_output"
    Int bam_size = ceil(size(bam, "GB"))
    Int inflation_factor = if (bam_size > 100) then 10 else 5
    Int disk_size = if (inflation_factor * bam_size > 20) then inflation_factor * bam_size else 20

    command <<<
        set -euxo pipefail

        num_core=$(grep -c '^processor' /proc/cpuinfo)
        mkdir -p "~{output_root}"

        export MONITOR_MOUNT_POINT="/mnt/disks/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        monitor_job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=~{model_type} \
            --ref=~{ref_fa} \
            ~{true='--haploid_contigs ' false='' defined(haploid_contigs)}~{select_first([haploid_contigs, ""])} \
            ~{true='--par_regions_bed ' false='' defined(par_regions_bed)}~{select_first([par_regions_bed, ""])} \
            --reads=~{bam} \
            --regions "~{sep=' ' regions}" \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}" || cat resources.log

        if ps -p "${monitor_job_id}" > /dev/null; then kill "${monitor_job_id}"; fi
    >>>

    output {
        File resource_usage_log = "resources.log"
        File vcf = "~{output_root}/~{prefix}.vcf.gz"
        File vcf_idx = "~{output_root}/~{prefix}.vcf.gz.tbi"
        File gvcf = "~{output_root}/~{prefix}.g.vcf.gz"
        File gvcf_idx = "~{output_root}/~{prefix}.g.vcf.gz.tbi"
        File visual_report = "~{output_root}/~{prefix}.visual_report.html"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: threads,
        mem_gb: memory,
        disk_gb: disk_size,
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
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
        zones: zones
    }
}

task RunDeepVariantGpu {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        Array[String] regions
        String? haploid_contigs
        File? par_regions_bed
        String model_type
        Int threads
        Int memory
        String zones
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String output_root = "/mnt/disks/cromwell_root/dv_output"
    Int bam_size = ceil(size(bam, "GB"))
    Int inflation_factor = if (bam_size > 100) then 10 else 5
    Int disk_size = if (inflation_factor * bam_size > 100) then inflation_factor * bam_size else 100
    Int gpu_cpu = if (threads > 12) then 12 else threads
    Int gpu_memory = if (memory > 64) then 64 else memory

    command <<<
        set -euxo pipefail

        num_core=$(grep -c '^processor' /proc/cpuinfo)
        mkdir -p "~{output_root}"

        export MONITOR_MOUNT_POINT="/mnt/disks/cromwell_root/"
        bash vm_local_monitoring_script.sh &> resources.log &
        monitor_job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')
        gpustat -a -i 1 &> gpu.usages.log &
        gpu_monitor_job_id=$(ps -aux | grep -F 'gpustat' | head -1 | awk '{print $2}')

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=~{model_type} \
            --ref=~{ref_fa} \
            ~{true='--haploid_contigs ' false='' defined(haploid_contigs)}~{select_first([haploid_contigs, ""])} \
            ~{true='--par_regions_bed ' false='' defined(par_regions_bed)}~{select_first([par_regions_bed, ""])} \
            --reads=~{bam} \
            --regions "~{sep=' ' regions}" \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}" || cat resources.log

        if ps -p "${monitor_job_id}" > /dev/null; then kill "${monitor_job_id}"; fi
        if ps -p "${gpu_monitor_job_id}" > /dev/null; then kill "${gpu_monitor_job_id}"; fi
    >>>

    output {
        File resource_usage_log = "resources.log"
        File gpu_usage_log = "gpu.usages.log"
        File vcf = "~{output_root}/~{prefix}.vcf.gz"
        File vcf_idx = "~{output_root}/~{prefix}.vcf.gz.tbi"
        File gvcf = "~{output_root}/~{prefix}.g.vcf.gz"
        File gvcf_idx = "~{output_root}/~{prefix}.g.vcf.gz.tbi"
        File visual_report = "~{output_root}/~{prefix}.visual_report.html"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: gpu_cpu,
        mem_gb: gpu_memory,
        disk_gb: disk_size,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
        zones: zones
        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
    }
}

task VisualizeResourceUsage {
    input {
        File resource_log
        String output_pdf_name
        String plot_title
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        /opt/plot.resources.R "~{resource_log}" "~{output_pdf_name}" "~{plot_title}"
    >>>

    output {
        File plot_pdf = "~{output_pdf_name}"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
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
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}
