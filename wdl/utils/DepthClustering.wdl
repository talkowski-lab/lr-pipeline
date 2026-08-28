version 1.0

import "Structs.wdl"

workflow DepthClustering {
    input {
        File depth_vcf
        String output_prefix
        String variant_prefix
        File pedigree

        File contig_list
        File? contig_subset_list

        File ref_fa
        File ref_fai
        File ref_dict

        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        # CreatePloidyTableFromPed
        String chr_x = "chrX"
        String chr_y = "chrY"

        # SVCluster
        Boolean fast_mode = true
        String clustering_algorithm = "SINGLE_LINKAGE"
        Boolean? enable_cnv
        Boolean? default_no_call
        Boolean? omit_members
        String? algorithm
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

        RuntimeAttr? runtime_attr_create_ploidy_table
        RuntimeAttr? runtime_attr_sv_cluster
        RuntimeAttr? runtime_attr_exclude_intervals
        RuntimeAttr? runtime_attr_gatk_to_svtk_vcf
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    call CreatePloidyTableFromPed {
        input:
            ped = pedigree,
            contig_list = contig_list,
            chr_x = chr_x,
            chr_y = chr_y,
            output_prefix = "~{output_prefix}",
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_create_ploidy_table
    }

    Array[String] contigs = read_lines(select_first([contig_subset_list, contig_list]))
    scatter (contig in contigs) {
        call SVCluster {
            input:
                vcf = depth_vcf,
                output_prefix = "~{contig}-depth_clustered",
                contig = contig,
                ploidy_table = CreatePloidyTableFromPed.ploidy_table,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                fast_mode = fast_mode,
                enable_cnv = enable_cnv,
                default_no_call = default_no_call,
                omit_members = omit_members,
                algorithm = algorithm,
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
                    overlap_fraction = exclude_overlap_fraction,
                    ref_fai = ref_fai,
                    output_prefix = "~{contig}-depth-intervals_excluded",
                    intervals = select_first([exclude_intervals]),
                    intervals_index = select_first([exclude_intervals]) + ".tbi",
                    docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_exclude_intervals
            }
        }

        call GatkToSvtkVcf {
            input:
                vcf = select_first([ExcludeIntervalsByIntervalOverlap.filtered_vcf, SVCluster.clustered_vcf]),
                output_prefix = "~{contig}-depth-svtk_formatted",
                script = gatk_to_svtk_script,
                contig_list = contig_list,
                set_pass = svtk_set_pass,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_gatk_to_svtk_vcf
        }
    }

    call ConcatVCFs {
        input:
            vcfs = GatkToSvtkVcf.svtk_vcf,
            vcf_idxs = GatkToSvtkVcf.svtk_vcf_index,
            output_prefix = "~{output_prefix}-depth",
            docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    output {
        File ploidy_table = CreatePloidyTableFromPed.ploidy_table
        File clustered_vcf = ConcatVCFs.concat_vcf
        File clustered_vcf_index = ConcatVCFs.concat_vcf_index
    }
}

task CreatePloidyTableFromPed {
    input {
        File ped
        File contig_list
        String chr_x
        String chr_y
        String output_prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python '/opt/sv-pipeline/scripts/ploidy_table_from_ped.py' \
            --ped ~{ped} \
            --out '~{output_prefix}.tsv' \
            --contigs '~{contig_list}' \
            --chr-x '~{chr_x}' \
            --chr-y '~{chr_y}'
    >>>

    output {
        File ploidy_table = "~{output_prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(ped, "GB") * 2) + 32,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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

task SVCluster {
    input {
        File vcf
        String output_prefix
        String contig
        File ploidy_table
        File ref_fa
        File ref_fai
        File ref_dict
        Boolean? fast_mode
        Boolean? enable_cnv
        Boolean? default_no_call
        Boolean? omit_members
        String? algorithm
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

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 800)

    command <<<
        set -euo pipefail

        gatk --java-options '-Xmx~{command_mem_mb}m' SVCluster \
            --variant '~{vcf}' \
            --output '~{output_prefix}.vcf.gz' \
            --reference '~{ref_fa}' \
            --ploidy-table '~{ploidy_table}' \
            --intervals '~{contig}' \
            ~{true="--fast-mode" false="" fast_mode} \
            ~{true="--enable-cnv" false="" enable_cnv} \
            ~{true="--default-no-call" false="" default_no_call} \
            ~{true="--omit-members" false="" omit_members} \
            ~{"--variant-prefix '" + variant_prefix + "'"} \
            ~{"--algorithm '" + algorithm + "'"} \
            ~{"--defrag-padding-fraction '" + defrag_padding_fraction + "'"} \
            ~{"--defrag-sample-overlap '" + defrag_sample_overlap + "'"} \
            ~{"--depth-sample-overlap '" + depth_sample_overlap + "'"} \
            ~{"--depth-interval-overlap '" + depth_interval_overlap + "'"} \
            ~{"--depth-size-similarity '" + depth_size_similarity + "'"} \
            ~{"--depth-breakend-window '" + depth_breakend_window + "'"} \
            ~{"--breakpoint-summary-strategy '" + breakpoint_summary_strategy + "'"}
    >>>

    output {
        File clustered_vcf = "~{output_prefix}.vcf.gz"
        File clustered_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB") * 2 + size(ref_fa, "GB")) + 32,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
        Float overlap_fraction
        File intervals
        File intervals_index
        File ref_fai
        String output_prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cut -f1,2 '~{ref_fai}' > genome.file
        bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\n' '~{vcf}' > variants.bed
        bedtools coverage -sorted -g genome.file -f ~{overlap_fraction} -a variants.bed
            -b '~{intervals}' \
            | awk -F"\t" '$6>0' \
            | cut -f4 > excluded_vids.list
        bcftools view --include '%ID!=@excluded_vids.list' --output-type z \
            --output '~{output_prefix}.vcf.gz' '~{vcf}'
        tabix '~{output_prefix}.vcf.gz'
    >>>

    output {
        File filtered_vcf = "~{output_prefix}.vcf.gz"
        File filtered_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) * 3 + 32,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
        File? script
        File contig_list
        Boolean set_pass
        String output_prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python '~{default="/opt/sv-pipeline/scripts/format_gatk_vcf_for_svtk.py" script}' \
            --vcf '~{vcf}' \
            --out '~{output_prefix}.vcf.gz' \
            --source depth \
            --contigs '~{contig_list}' \
            --remove-formats CN \
            ~{if set_pass then "--set-pass" else ""}
        tabix '~{output_prefix}.vcf.gz'
    >>>

    output {
        File svtk_vcf = "~{output_prefix}.vcf.gz"
        File svtk_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcf, "GB")) * 2 + 32,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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

task ConcatVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools concat --no-version --naive --output-type z --file-list '~{write_lines(vcfs)}' \
            --output '~{output_prefix}.vcf.gz'
        tabix '~{output_prefix}.vcf.gz'
    >>>

    output {
        File concat_vcf = "~{output_prefix}.vcf.gz"
        File concat_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(vcfs, "GB") * 3) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
