version 1.0

import "../utils/Structs.wdl"

# CompareMethylationQC — merge two MethylationQC.wdl runs (e.g. PacBio and ONT) into
# platform-labelled tables and side-by-side figures. The headline output is
# compare_site_sd: the per-CpG SD of methylation across samples by platform, raw and
# after downsampling every sample to the same depth.
#
# Both runs must use the same reference, contigs and matched_coverage. Inputs are the
# outputs of each run, passed as [A, B] pairs (default labels PacBio, ONT).
workflow CompareMethylationQC {
    input {
        String label_a = "PacBio"
        String label_b = "ONT"
        String prefix

        File sample_qc_a
        File sample_qc_b
        File sample_hists_a
        File sample_hists_b
        File cohort_summary_a
        File cohort_summary_b
        File per_site_stats_a
        File per_site_stats_b
        File site_sd_hist_a
        File site_sd_hist_b
        File site_sd_by_mean_a
        File site_sd_by_mean_b
        File site_mean_sd_grid_a
        File site_mean_sd_grid_b
        File site_excess_var_hist_a
        File site_excess_var_hist_b
        File pca_matrix_a
        File pca_matrix_b

        Int min_samples_per_site = 10

        String methylation_docker = "quay.io/ymostovoy/lr-methylation:latest"

        RuntimeAttr? runtime_attr_compare
    }

    parameter_meta {
        label_a: "Platform label for the *_a inputs; used in column names and legends."
        label_b: "Platform label for the *_b inputs."
        prefix: "Basename for all outputs."
        sample_qc_a: "MethylationQC sample_qc_tsv output for platform A (same pattern for every *_a / *_b input: pass the identically named MethylationQC output of each run)."
        per_site_stats_a: "MethylationQC per_site_stats output for platform A. Joined with B on (chrom, start) for the same-site SD comparison."
        pca_matrix_a: "MethylationQC pca_matrix output for platform A. Shared sites are used for a joint PCA."
        min_samples_per_site: "Minimum samples contributing to a site in each cohort for the joint per-site comparison."
        methylation_docker: "Image containing scripts/methylation/methylation_qc.py (dockerfiles/Dockerfile.Methylation)."
        runtime_attr_compare: "Optional RuntimeAttr override. Defaults: 2 CPU, 32 GiB."
    }

    call Compare {
        input:
            label_a = label_a,
            label_b = label_b,
            prefix = prefix,
            sample_qc = [sample_qc_a, sample_qc_b],
            sample_hists = [sample_hists_a, sample_hists_b],
            cohort_summary = [cohort_summary_a, cohort_summary_b],
            per_site_stats = [per_site_stats_a, per_site_stats_b],
            site_sd_hist = [site_sd_hist_a, site_sd_hist_b],
            site_sd_by_mean = [site_sd_by_mean_a, site_sd_by_mean_b],
            site_mean_sd_grid = [site_mean_sd_grid_a, site_mean_sd_grid_b],
            site_excess_var_hist = [site_excess_var_hist_a, site_excess_var_hist_b],
            pca_matrix = [pca_matrix_a, pca_matrix_b],
            min_samples_per_site = min_samples_per_site,
            docker = methylation_docker,
            runtime_attr_override = runtime_attr_compare
    }

    output {
        Array[File] merged_tables = Compare.merged_tables
        File joint_site_stats = Compare.joint_site_stats
        File joint_sd_grid = Compare.joint_sd_grid
        File joint_delta_sd_summary = Compare.joint_delta_sd_summary
        File joint_pca_coords = Compare.joint_pca_coords
        File joint_pca_variance = Compare.joint_pca_variance
        Array[File] figures = Compare.figures
        File figures_tar = Compare.figures_tar
    }
}

task Compare {
    input {
        String label_a
        String label_b
        String prefix
        Array[File] sample_qc
        Array[File] sample_hists
        Array[File] cohort_summary
        Array[File] per_site_stats
        Array[File] site_sd_hist
        Array[File] site_sd_by_mean
        Array[File] site_mean_sd_grid
        Array[File] site_excess_var_hist
        Array[File] pca_matrix
        Int min_samples_per_site
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/methylation/methylation_qc.py compare \
            --label-a '~{label_a}' \
            --label-b '~{label_b}' \
            --sample-qc ~{sep=' ' sample_qc} \
            --sample-hists ~{sep=' ' sample_hists} \
            --cohort-summary ~{sep=' ' cohort_summary} \
            --per-site-stats ~{sep=' ' per_site_stats} \
            --site-sd-hist ~{sep=' ' site_sd_hist} \
            --site-sd-by-mean ~{sep=' ' site_sd_by_mean} \
            --site-mean-sd-grid ~{sep=' ' site_mean_sd_grid} \
            --site-excess-var-hist ~{sep=' ' site_excess_var_hist} \
            --pca-matrix ~{sep=' ' pca_matrix} \
            --min-samples-per-site ~{min_samples_per_site} \
            --prefix ~{prefix} \
            --outdir out

        mkdir figures
        mv out/*.png out/*.pdf figures/
        tar -czf ~{prefix}.compare_figures.tar.gz figures
    >>>

    output {
        Array[File] merged_tables = glob("out/~{prefix}.merged_*.tsv")
        File joint_site_stats = "out/~{prefix}.joint_site_stats.tsv.gz"
        File joint_sd_grid = "out/~{prefix}.joint_sd_grid.tsv"
        File joint_delta_sd_summary = "out/~{prefix}.joint_delta_sd_summary.tsv"
        File joint_pca_coords = "out/~{prefix}.joint.pca_coords.tsv"
        File joint_pca_variance = "out/~{prefix}.joint.pca_variance.tsv"
        Array[File] figures = glob("figures/*")
        File figures_tar = "~{prefix}.compare_figures.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 32,
        disk_gb: 3 * ceil(size(per_site_stats, "GB") + size(pca_matrix, "GB") + size(sample_hists, "GB")) + 20,
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
