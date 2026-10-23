version 1.0

import "Structs.wdl"

workflow CreateAnnotationPlots {
    meta {
        description: [
            "This utility draws the callset summary figures from the plotting Parquet tables that `SummarizeAnnotations` produces with `create_plotting` enabled, and returns them as one tarball of PNGs. Site counts and short-read overlap are drawn from the site table, sites per genome from the sample and allele tables, and trio inheritance from the de novo table, each only when its table is given.",
            "The allele-frequency and size distributions and the Hardy-Weinberg plots (overall, and by variant type, size and genomic region) are drawn from the per-variant list, which can hold hundreds of millions of rows. Each of these runs as its own task that reads only the columns its plot needs, so no single VM holds the whole table for every plot, and its default memory scales with the size of the list."
        ]
    }

    parameter_meta {
        plotting_sites_parquet: "Site-level plotting table (`plotting_sites_parquet` from `SummarizeAnnotations`)."
        plotting_samples_parquet: "Per-sample plotting table (`plotting_samples_parquet` from `SummarizeAnnotations`), for the sites-per-genome plot."
        plotting_alleles_parquet: "Per-allele plotting table (`plotting_alleles_parquet` from `SummarizeAnnotations`), adding allele counts to the sites-per-genome plot."
        plotting_denovo_parquet: "De novo transmission table (`plotting_denovo_parquet` from `SummarizeAnnotations`), for the trio inheritance plot. Leave unset when the cohort has no trios."
        plotting_variant_list_parquet: "Per-variant list (`plotting_variant_list_parquet` from `SummarizeAnnotations`), for the distribution and Hardy-Weinberg plots."
        random_seed: "Seed for every random subsample taken while plotting."
        hwe_subsample_count: "Maximum number of variants drawn per Hardy-Weinberg panel."
        plots_tarball: "Tarball of every plot drawn, one PNG per plot."
    }

    input {
        File plotting_sites_parquet
        File? plotting_samples_parquet
        File? plotting_alleles_parquet
        File? plotting_denovo_parquet
        File? plotting_variant_list_parquet
        String prefix

        Int random_seed = 42
        Int hwe_subsample_count = 500000

        String utils_docker

        RuntimeAttr? runtime_attr_create_site_plots
        RuntimeAttr? runtime_attr_create_genome_plots
        RuntimeAttr? runtime_attr_create_inheritance_plots
        RuntimeAttr? runtime_attr_create_variant_list_plots
        RuntimeAttr? runtime_attr_archive_plots
    }

    Array[String] variant_list_plots = ["distribution_af", "distribution_size", "hwe_overall", "hwe_type", "hwe_size", "hwe_region"]

    call CreatePlots as CreateSitePlots {
        input:
            plots = ["site_counts", "srwgs_overlap"],
            sites_parquet = plotting_sites_parquet,
            random_seed = random_seed,
            hwe_subsample_count = hwe_subsample_count,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_site_plots
    }

    if (defined(plotting_samples_parquet)) {
        call CreatePlots as CreateGenomePlots {
            input:
                plots = ["sites_per_genome"],
                samples_parquet = plotting_samples_parquet,
                alleles_parquet = plotting_alleles_parquet,
                random_seed = random_seed,
                hwe_subsample_count = hwe_subsample_count,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_genome_plots
        }
    }

    if (defined(plotting_denovo_parquet)) {
        call CreatePlots as CreateInheritancePlots {
            input:
                plots = ["inheritance"],
                denovo_parquet = plotting_denovo_parquet,
                random_seed = random_seed,
                hwe_subsample_count = hwe_subsample_count,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_inheritance_plots
        }
    }

    if (defined(plotting_variant_list_parquet)) {
        scatter (plot in variant_list_plots) {
            call CreatePlots as CreateVariantListPlots {
                input:
                    plots = [plot],
                    sites_parquet = plotting_sites_parquet,
                    variant_list_parquet = plotting_variant_list_parquet,
                    random_seed = random_seed,
                    hwe_subsample_count = hwe_subsample_count,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_create_variant_list_plots
            }
        }
    }

    call ArchivePlots {
        input:
            pngs = flatten([
                CreateSitePlots.pngs,
                select_first([CreateGenomePlots.pngs, []]),
                select_first([CreateInheritancePlots.pngs, []]),
                flatten(select_first([CreateVariantListPlots.pngs, []]))
            ]),
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_archive_plots
    }

    output {
        File plots_tarball = ArchivePlots.plots_tarball
    }
}

task CreatePlots {
    input {
        Array[String] plots
        File? sites_parquet
        File? samples_parquet
        File? alleles_parquet
        File? denovo_parquet
        File? variant_list_parquet
        Int random_seed
        Int hwe_subsample_count
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Float variant_list_size = size(variant_list_parquet, "GB")

    command <<<
        set -euo pipefail

        mkdir plots

        python3 /opt/scripts/qc/plot_annotation_summaries.py \
            --plots ~{sep=" " plots} \
            ~{"--sites_parquet " + sites_parquet} \
            ~{"--samples_parquet " + samples_parquet} \
            ~{"--alleles_parquet " + alleles_parquet} \
            ~{"--denovo_parquet " + denovo_parquet} \
            ~{"--variant_list_parquet " + variant_list_parquet} \
            --random-seed ~{random_seed} \
            --hwe-subsample-count ~{hwe_subsample_count} \
            --output_dir plots
    >>>

    output {
        Array[File] pngs = glob("plots/*.png")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4 + 12 * ceil(variant_list_size),
        disk_gb: 3 * ceil(variant_list_size + size([sites_parquet, samples_parquet, alleles_parquet, denovo_parquet], "GB")) + 10,
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

task ArchivePlots {
    input {
        Array[File] pngs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir "~{prefix}.plots"
        while read -r png; do
            cp "$png" "~{prefix}.plots/"
        done < ~{write_lines(pngs)}

        tar -czf "~{prefix}.plots.tar.gz" "~{prefix}.plots"
    >>>

    output {
        File plots_tarball = "~{prefix}.plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(pngs, "GB")) + 5,
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
