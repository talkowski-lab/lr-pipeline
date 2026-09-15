version 1.0

struct HumanReferenceBundle {
    File fasta
    File fai
    File dict

    File tandem_repeat_bed
    File PAR_bed

    File? size_balanced_scatter_interval_ids
    File? size_balanced_scatter_intervallists_locators

    File intervallists_autosomes
    File intervallists_allosomes

    File chromosome_ploidy_priors

    String mt_chr_name

    File? haplotype_map
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct SmallVarJobConfig {
    String? haploid_contigs

    Int dv_threads
    Int dv_memory
    Boolean use_gpu

    Boolean run_clair3
    Boolean phase_and_tag
    Boolean use_margin_for_tagging

    String? gcp_zones
}
