version 1.0

import "Structs.wdl"
import "Helpers.wdl"

workflow DepthClustering {
    input {
        File depth_vcf
        File depth_vcf_idx
        File ploidy_table
        String prefix
        String variant_prefix

        File contig_list
        File? contig_subset_list

        File ref_fa
        File ref_fai
        File ref_dict

        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        # SVCluster
        Boolean fast_mode = true
        String clustering_algorithm = "SINGLE_LINKAGE"
        Boolean? enable_cnv
        Boolean? default_no_call
        Boolean? omit_members
        String? breakpoint_summary_strategy
        Float? defrag_padding_fraction
        Float? defrag_sample_overlap
        Float depth_sample_overlap = 0
        Float depth_interval_overlap = 0.8
        Float? depth_size_similarity
        Int depth_breakend_window = 10000000

        # ExcludeIntervalsByIntervalOverlap
        File? exclude_intervals
        Float exclude_overlap_fraction = 0.5

        # GatkToSvtkVcf
        File? gatk_to_svtk_script
        Boolean svtk_set_pass = false

        RuntimeAttr? runtime_attr_sv_cluster
        RuntimeAttr? runtime_attr_exclude_intervals
        RuntimeAttr? runtime_attr_gatk_to_svtk_vcf
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    Array[String] contigs = read_lines(select_first([contig_subset_list, contig_list]))
    scatter (contig in contigs) {
        call SVCluster {
            input:
                vcf = depth_vcf,
                vcf_idx = depth_vcf_idx,
                prefix = "~{prefix}-~{contig}-depth_clustered",
                contig = contig,
                ploidy_table = ploidy_table,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                fast_mode = fast_mode,
                enable_cnv = enable_cnv,
                default_no_call = default_no_call,
                omit_members = omit_members,
                clustering_algorithm = clustering_algorithm,
                breakpoint_summary_strategy = breakpoint_summary_strategy,
                defrag_padding_fraction = defrag_padding_fraction,
                defrag_sample_overlap = defrag_sample_overlap,
                depth_sample_overlap = depth_sample_overlap,
                depth_interval_overlap = depth_interval_overlap,
                depth_size_similarity = depth_size_similarity,
                depth_breakend_window = depth_breakend_window,
                variant_prefix = "~{variant_prefix}_depth_~{contig}_",
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_sv_cluster
        }

        if (defined(exclude_intervals)) {
            call ExcludeIntervalsByIntervalOverlap {
                input:
                    vcf = SVCluster.clustered_vcf,
                    vcf_idx = SVCluster.clustered_vcf_idx,
                    overlap_fraction = exclude_overlap_fraction,
                    ref_fai = ref_fai,
                    prefix = "~{prefix}-~{contig}-depth-intervals_excluded",
                    intervals = select_first([exclude_intervals]),
                    docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_exclude_intervals
            }
        }

        call GatkToSvtkVcf {
            input:
                vcf = select_first([ExcludeIntervalsByIntervalOverlap.filtered_vcf, SVCluster.clustered_vcf]),
                vcf_idx = select_first([ExcludeIntervalsByIntervalOverlap.filtered_vcf_idx, SVCluster.clustered_vcf_idx]),
                prefix = "~{prefix}-~{contig}-depth-svtk_formatted",
                script = gatk_to_svtk_script,
                contig_list = contig_list,
                set_pass = svtk_set_pass,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_gatk_to_svtk_vcf
        }
    }

    RuntimeAttr default_attr_concat_vcfs = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(GatkToSvtkVcf.svtk_vcf, "GB") * 3) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr concat_attr = select_first([runtime_attr_concat_vcfs, default_attr_concat_vcfs])

    call Helpers.ConcatVcfs as ConcatVCFs {
        input:
            vcfs = GatkToSvtkVcf.svtk_vcf,
            vcf_idxs = GatkToSvtkVcf.svtk_vcf_idx,
            allow_overlaps = false,
            naive = true,
            no_version = true,
            no_address = true,
            prefix = "~{prefix}-depth",
            docker = sv_base_mini_docker,
            runtime_attr_override = concat_attr
    }

    output {
        File clustered_vcf = ConcatVCFs.concat_vcf
        File clustered_vcf_idx = ConcatVCFs.concat_vcf_idx
    }
}

task SVCluster {
    input {
        File vcf
        File vcf_idx
        String prefix
        String contig
        File ploidy_table
        File ref_fa
        File ref_fai
        File ref_dict
        Boolean? fast_mode
        Boolean? enable_cnv
        Boolean? default_no_call
        Boolean? omit_members
        String clustering_algorithm
        String? breakpoint_summary_strategy
        Float? defrag_padding_fraction
        Float? defrag_sample_overlap
        Float? depth_sample_overlap
        Float? depth_interval_overlap
        Float? depth_size_similarity
        Int? depth_breakend_window
        String? variant_prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcf: { localization_optional: true }
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    command <<<
        set -euo pipefail

        gatk --java-options '-Xmx~{command_mem_mb}m' SVCluster \
            --variant '~{vcf}' \
            --output '~{prefix}.vcf.gz' \
            --reference '~{ref_fa}' \
            --ploidy-table '~{ploidy_table}' \
            --intervals '~{contig}' \
            ~{true="--fast-mode" false="" fast_mode} \
            ~{true="--enable-cnv" false="" enable_cnv} \
            ~{true="--default-no-call" false="" default_no_call} \
            ~{true="--omit-members" false="" omit_members} \
            ~{"--variant-prefix '" + variant_prefix + "'"} \
            --algorithm '~{clustering_algorithm}' \
            ~{"--defrag-padding-fraction '" + defrag_padding_fraction + "'"} \
            ~{"--defrag-sample-overlap '" + defrag_sample_overlap + "'"} \
            ~{"--depth-sample-overlap '" + depth_sample_overlap + "'"} \
            ~{"--depth-interval-overlap '" + depth_interval_overlap + "'"} \
            ~{"--depth-size-similarity '" + depth_size_similarity + "'"} \
            ~{"--depth-breakend-window '" + depth_breakend_window + "'"} \
            ~{"--breakpoint-summary-strategy '" + breakpoint_summary_strategy + "'"}
    >>>

    output {
        File clustered_vcf = "~{prefix}.vcf.gz"
        File clustered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB") * 2 + size(ref_fa, "GB")) + 32,
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

task ExcludeIntervalsByIntervalOverlap {
    input {
        File vcf
        File vcf_idx
        Float overlap_fraction
        File intervals
        File ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cut -f1,2 '~{ref_fai}' > genome.file
        bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\n' '~{vcf}' > variants.bed
        bedtools coverage -sorted -g genome.file -f ~{overlap_fraction} -a variants.bed \
            -b '~{intervals}' \
            | awk -F"\t" '$6>0' \
            | cut -f4 > excluded_vids.list
        bcftools view --include '%ID!=@excluded_vids.list' --output-type z \
            --output '~{prefix}.vcf.gz' '~{vcf}'
        tabix '~{prefix}.vcf.gz'
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([vcf, intervals], "GB") * 3) + 32,
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
        noAddress: true
    }
}

task GatkToSvtkVcf {
    input {
        File vcf
        File vcf_idx
        File? script
        File contig_list
        Boolean set_pass
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python '~{default="/opt/sv-pipeline/scripts/format_gatk_vcf_for_svtk.py" script}' \
            --vcf '~{vcf}' \
            --out '~{prefix}.vcf.gz' \
            --source depth \
            --contigs '~{contig_list}' \
            --remove-formats CN \
            ~{if set_pass then "--set-pass" else ""}
        tabix '~{prefix}.vcf.gz'
    >>>

    output {
        File svtk_vcf = "~{prefix}.vcf.gz"
        File svtk_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) * 2 + 32,
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
        noAddress: true
    }
}
