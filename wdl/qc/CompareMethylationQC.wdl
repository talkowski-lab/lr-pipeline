version 1.0

import "../utils/Structs.wdl"

# CompareMethylationQC — merge two MethylationQC.wdl runs (e.g. PacBio and ONT) into
# platform-labelled tables and side-by-side figures. The headline figure is
# compare_site_sd: the per-CpG SD of methylation across samples by platform, raw and
# after downsampling every sample to the same depth.
#
# Inputs are the stats_tar output of each run. Both runs must use the same reference,
# contigs and matched_coverage. Outputs mirror MethylationQC: one stats tarball, one
# figures tarball.
workflow CompareMethylationQC {
    input {
        File stats_tar_a
        File stats_tar_b
        String prefix

        String? label_a
        String? label_b
        Int min_samples_per_site = 10

        String methylation_docker = "quay.io/ymostovoy/lr-methylation:latest"

        RuntimeAttr? runtime_attr_compare
    }

    parameter_meta {
        stats_tar_a: "stats_tar output (<prefix>.stats_for_comparison.tar.gz) of the first MethylationQC run, e.g. PacBio."
        stats_tar_b: "stats_tar output of the second MethylationQC run, e.g. ONT."
        prefix: "Basename for all outputs."
        label_a: "Optional label for run A in column names and legends. Default: its platform (PacBio/ONT); if both runs share a platform, platform_betasource (e.g. PacBio_model vs PacBio_counts)."
        label_b: "Optional label for run B. Default as for label_a."
        min_samples_per_site: "Minimum samples contributing to a site in each cohort for the joint per-site comparison. Automatically capped at each cohort's size (minimum 2)."
        methylation_docker: "Image containing scripts/methylation/methylation_qc.py (dockerfiles/Dockerfile.Methylation)."
        runtime_attr_compare: "Optional RuntimeAttr override. Defaults: 2 CPU, 32 GiB."
    }

    call Compare {
        input:
            stats_tar_a = stats_tar_a,
            stats_tar_b = stats_tar_b,
            prefix = prefix,
            label_a = label_a,
            label_b = label_b,
            min_samples_per_site = min_samples_per_site,
            docker = methylation_docker,
            runtime_attr_override = runtime_attr_compare
    }

    output {
        File stats_tar = Compare.stats_tar
        File figures_tar = Compare.figures_tar
    }
}

task Compare {
    input {
        File stats_tar_a
        File stats_tar_b
        String prefix
        String? label_a
        String? label_b
        Int min_samples_per_site
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/methylation/methylation_qc.py compare \
            --stats ~{stats_tar_a} ~{stats_tar_b} \
            ~{"--label-a '" + label_a + "'"} \
            ~{"--label-b '" + label_b + "'"} \
            --min-samples-per-site ~{min_samples_per_site} \
            --prefix ~{prefix} \
            --outdir ~{prefix}.comparison_stats

        mkdir ~{prefix}.comparison_figures
        mv ~{prefix}.comparison_stats/*.png ~{prefix}.comparison_stats/*.pdf \
            ~{prefix}.comparison_figures/
        tar -czf ~{prefix}.comparison_stats.tar.gz ~{prefix}.comparison_stats
        tar -czf ~{prefix}.comparison_figures.tar.gz ~{prefix}.comparison_figures
    >>>

    output {
        File stats_tar = "~{prefix}.comparison_stats.tar.gz"
        File figures_tar = "~{prefix}.comparison_figures.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 32,
        disk_gb: 4 * ceil(size(stats_tar_a, "GB") + size(stats_tar_b, "GB")) + 20,
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
