# Derived from wdl/tasks/VariantCalling/DeepVariant.wdl of broadinstitute/long-read-pipelines (branch sh_update_wgs_callers_outliers)
# https://github.com/broadinstitute/long-read-pipelines/blob/main/wdl/tasks/VariantCalling/DeepVariant.wdl

version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow DeepVariant {
    meta {
        description: [
            "This tool follows the DeepVariant path of the linked source workflow: it reads the source reference-bundle and small-variant-options JSON files, subsets the input BAM into its size-balanced shards, calls DeepVariant, then merges the VCFs and gVCFs. The same region list controls BAM subsetting and DeepVariant `--regions`, preventing duplicate off-shard zero-depth gVCF blocks. It performs no explicit GCS-directory copy; map its declared outputs directly to attributes in the main Terra entity table."
        ]
    }

    parameter_meta {
        bam: "Aligned whole-genome BAM and index."
        bai: "Aligned whole-genome BAM and index."
        sex: "Biological sex; `M` enables the source reference bundle's haploid-contig and PAR settings."
        model_for_dv_andor_pepper: "DeepVariant model, such as `PACBIO` or `ONT_R104`."
        ref_bundle_json_file: "Source-compatible reference bundle JSON. Its size-balanced shard manifests are required."
        small_variant_calling_options_json: "Source-compatible small-variant options JSON, providing DeepVariant threads, memory and haploid contigs."
        gcp_zones: "Placement zones."
        gvcf: "Merged gVCF."
        gvcf_idx: "Index for gvcf."
        vcf: "Merged VCF."
        vcf_idx: "Index for vcf."
        visual_reports: "Per-shard diagnostic outputs."
    }

    input {
        File bam
        File bai
        String sex
        String prefix
        String model_for_dv_andor_pepper
        File ref_bundle_json_file
        File small_variant_calling_options_json
        Array[String] gcp_zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]

        String deepvariant_docker
        String utils_docker

        RuntimeAttr? runtime_attr_collapse_zones
        RuntimeAttr? runtime_attr_subset_bam
        RuntimeAttr? runtime_attr_run_deepvariant
        RuntimeAttr? runtime_attr_merge_gvcfs
        RuntimeAttr? runtime_attr_merge_vcfs
    }

    HumanReferenceBundle ref_bundle = read_json(ref_bundle_json_file)
    SmallVarJobConfig small_variant_options = read_json(small_variant_calling_options_json)

    call CollapseArrayOfStrings as CollapseZones {
        input:
            values = gcp_zones,
            delimiter = " ",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_collapse_zones
    }
    String zones = CollapseZones.collapsed

    File shard_ids_file = select_first([ref_bundle.size_balanced_scatter_interval_ids])
    File shard_region_locators_file = select_first([ref_bundle.size_balanced_scatter_intervallists_locators])
    Array[String] shard_ids = read_lines(shard_ids_file)
    Array[String] shard_region_paths = read_lines(shard_region_locators_file)
    Array[Pair[String, String]] shard_plan = zip(shard_ids, shard_region_paths)

    scatter (shard in shard_plan) {
        String shard_id = shard.left

        if (shard_id != "alts") {
            File shard_region_file = shard.right
            Array[String] shard_regions = read_lines(shard_region_file)

            Boolean is_t2t_chr_x = length(shard_plan) < 20 && shard_id == "18_X"
            Boolean is_t2t_shard_8 = length(shard_plan) < 20 && shard_id == "9_15"
            Boolean is_grch38_chr_1 = length(shard_plan) > 20 && shard_id == "1-p"
            Boolean is_grch38_shard_3 = length(shard_plan) > 20 && shard_id == "11-q_17-q"
            Boolean needs_extra_memory = is_t2t_chr_x || is_t2t_shard_8 || is_grch38_chr_1 || is_grch38_shard_3
            Int cpu_memory = if (needs_extra_memory) then 48 else small_variant_options.dv_memory

            call Helpers.SubsetBamToRegions as SubsetBam {
                input:
                    bam = bam,
                    bai = bai,
                    region_file = shard_region_file,
                    prefix = "~{prefix}.~{shard_id}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_bam
            }

            call RunDeepVariant {
                input:
                    bam = SubsetBam.subset_bam,
                    bai = SubsetBam.subset_bai,
                    ref_fa = ref_bundle.fasta,
                    ref_fai = ref_bundle.fai,
                    regions = shard_regions,
                    sex = sex,
                    haploid_contigs = small_variant_options.haploid_contigs,
                    par_regions_bed = ref_bundle.PAR_bed,
                    model_type = model_for_dv_andor_pepper,
                    threads = small_variant_options.dv_threads,
                    memory = cpu_memory,
                    zones = zones,
                    prefix = "~{prefix}.~{shard_id}",
                    docker = deepvariant_docker,
                    runtime_attr_override = runtime_attr_run_deepvariant
            }
        }
    }

    call MergeAndSortVcfs as MergeGvcfs {
        input:
            vcfs = select_all(RunDeepVariant.gvcf),
            ref_fai = ref_bundle.fai,
            prefix = "~{prefix}.deepvariant.g",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_gvcfs
    }

    call MergeAndSortVcfs as MergeVcfs {
        input:
            vcfs = select_all(RunDeepVariant.vcf),
            ref_fai = ref_bundle.fai,
            prefix = "~{prefix}.deepvariant",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_vcfs
    }

    output {
        File gvcf = MergeGvcfs.vcf
        File gvcf_idx = MergeGvcfs.vcf_idx
        File vcf = MergeVcfs.vcf
        File vcf_idx = MergeVcfs.vcf_idx
        Array[File] visual_reports = select_all(RunDeepVariant.visual_report)
    }
}

task RunDeepVariant {
    input {
        File bam
        File bai
        File ref_fa
        File ref_fai
        Array[String] regions
        String sex
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
    Boolean use_haploid_contigs = sex == "M" && defined(haploid_contigs)
    Boolean use_par_regions_bed = sex == "M" && defined(par_regions_bed)

    command <<<
        set -euo pipefail

        num_core=$(grep -c '^processor' /proc/cpuinfo)
        mkdir -p "~{output_root}"

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=~{model_type} \
            --ref=~{ref_fa} \
            ~{true='--haploid_contigs ' false='' use_haploid_contigs}~{select_first([haploid_contigs, ""])} \
            ~{true='--par_regions_bed ' false='' use_par_regions_bed}~{select_first([par_regions_bed, ""])} \
            --reads=~{bam} \
            --regions "~{sep=' ' regions}" \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}"
    >>>

    output {
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

task CollapseArrayOfStrings {
    input {
        Array[String] values
        String delimiter
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        tr '\n' "~{delimiter}" < "~{write_lines(values)}" > result.txt
    >>>

    output {
        String collapsed = read_string("result.txt")
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

task MergeAndSortVcfs {
    input {
        Array[File] vcfs
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int input_size = ceil(size(vcfs, "GB"))
    Int disk_size = if (input_size > 100) then 5 * input_size else 375
    Int cores = 8
    Int memory = 48

    command <<<
        set -euo pipefail

        printf '%s\n' ~{sep=' ' vcfs} > input_vcfs.txt
        bcftools concat --naive --threads ~{cores - 1} -f input_vcfs.txt --output-type v -o concatenated.vcf.gz
        bcftools reheader --fai ~{ref_fai} -o reheadered.vcf.gz concatenated.vcf.gz
        bcftools sort --temp-dir sort_tmp --output-type z -o ~{prefix}.vcf.gz reheadered.vcf.gz
        bcftools index --tbi --force ~{prefix}.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: cores,
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
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: docker
    }
}
