version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

# MethylationQC — cohort QC, summary statistics and figures for CpG methylation BEDs
# from a single platform: pb-cpg-tools combined BEDs (PacBio) or modkit bedMethyl (ONT).
# Run once per platform, then merge two runs with CompareMethylationQC.wdl.
#
# Every sample is projected onto a fixed reference CpG index, so per-site statistics are
# additive: samples are processed in a scatter, summed in batches (tree reduction), and
# the batch partial sums are combined once. No task localizes the whole cohort.
#
# Per-site variance across samples is reported raw, decomposed into expected binomial
# (read-sampling) variance plus excess, and after hypergeometric downsampling of every
# sample to matched_coverage reads, so platform comparisons are not confounded by depth.
# Figures are drawn only from the emitted summary tables, so they can be replotted.
workflow MethylationQC {
    input {
        Array[File] methylation_beds
        Array[String] sample_ids
        String platform
        String prefix
        String pacbio_beta_source = "model"

        File ref_fa
        File ref_fai
        File? cpg_islands_bed
        Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                                 "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
                                 "chr22", "chrX", "chrY"]

        Int min_coverage = 10
        Int matched_coverage = 10
        Int min_samples_per_site = 10
        Int max_neighbor_distance = 50
        Int pca_site_stride = 100
        Int samples_per_batch = 100
        Float min_ref_cpg_match_frac = 0.8
        Float outlier_mad = 4.0

        String methylation_docker = "quay.io/ymostovoy/lr-methylation:latest"

        RuntimeAttr? runtime_attr_build_index
        RuntimeAttr? runtime_attr_process_sample
        RuntimeAttr? runtime_attr_accumulate
        RuntimeAttr? runtime_attr_finalize
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_plot
    }

    parameter_meta {
        methylation_beds: "Per-sample combined (unphased) methylation BEDs, aligned by index with sample_ids. PacBio: pb-cpg-tools *.combined.bed(.gz). ONT: modkit pileup bedMethyl (stranded or --combine-strands; must be CpG-only, i.e. run with --cpg)."
        sample_ids: "Sample IDs, same order as methylation_beds. Must be unique."
        platform: "PacBio or ONT. Selects the parser; every file is checked against it and the run fails on a format mismatch. ONT uses 5mC (code m) only; 5hmC rows are dropped."
        pacbio_beta_source: "PacBio only (ignored for ONT). Which per-site methylation value is used as the raw beta: model = pb-cpg-tools mod_score, the model-based value normally used; counts = est_mod_count / cov, the fraction of reads called methylated, the same kind of value modkit reports for ONT. Counts makes raw, excess and matched metrics directly comparable to ONT; model may look less variable because the model smooths low-depth sites. Depth-matched metrics always use counts. Run once with each and different prefixes to measure the model's smoothing."
        prefix: "Basename for all outputs, e.g. the cohort and platform name."
        ref_fa: "Uncompressed reference FASTA the BEDs were generated against. Used to build the reference CpG index."
        ref_fai: "FASTA index for ref_fa."
        cpg_islands_bed: "Optional CpG-island BED (chrom, start, end; e.g. UCSC cpgIslandExt) with contig names matching ref_fa. Enables island/shore (<=2 kb)/shelf (2-4 kb)/open-sea stratification."
        contigs: "Reference contigs to index. Cohort variance summaries use autosomes only; chrX/chrY feed per-sample sex-check metrics."
        min_coverage: "Minimum read depth for a sample's site to count in raw per-site stats and per-sample beta metrics."
        matched_coverage: "Depth every sample is downsampled to (sites with fewer reads are dropped) for the coverage-matched metrics. Use the same value for both platforms."
        min_samples_per_site: "Minimum number of samples contributing to a site for it to enter the binned site summaries."
        max_neighbor_distance: "Maximum bp between adjacent CpGs for the within-sample neighbour-discordance metric."
        pca_site_stride: "Every Nth reference CpG is kept for the sample PCA matrix (100 gives ~290k sites on hg38)."
        samples_per_batch: "Samples per accumulation task. Larger batches mean fewer, longer tasks."
        min_ref_cpg_match_frac: "Fail a sample if fewer than this fraction of its records fall on reference CpGs, which catches a wrong reference build or non-CpG modkit output."
        outlier_mad: "Samples more than this many robust SDs (MAD-scaled) from the cohort median on any key metric are flagged in sample_qc_tsv."
        methylation_docker: "Image containing scripts/methylation/methylation_qc.py (dockerfiles/Dockerfile.Methylation)."
        runtime_attr_build_index: "Optional RuntimeAttr override. Defaults: 1 CPU, 8 GiB."
        runtime_attr_process_sample: "Optional RuntimeAttr override. Defaults: 2 CPU, 16 GiB, preemptible."
        runtime_attr_accumulate: "Optional RuntimeAttr override. Defaults: 2 CPU, 12 GiB."
        runtime_attr_finalize: "Optional RuntimeAttr override. Defaults: 4 CPU, 32 GiB. PCA memory grows with sample count; raise for >5000 samples."
        runtime_attr_concat: "Optional RuntimeAttr override for TSV concatenation."
        runtime_attr_plot: "Optional RuntimeAttr override. Defaults: 1 CPU, 8 GiB."
    }

    call BuildCpgIndex {
        input:
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            cpg_islands_bed = cpg_islands_bed,
            contigs = contigs,
            prefix = prefix,
            docker = methylation_docker,
            runtime_attr_override = runtime_attr_build_index
    }

    scatter (i in range(length(methylation_beds))) {
        call ProcessSample {
            input:
                bed = methylation_beds[i],
                sample_id = sample_ids[i],
                platform = platform,
                pacbio_beta_source = pacbio_beta_source,
                cpg_index = BuildCpgIndex.cpg_index,
                min_coverage = min_coverage,
                matched_coverage = matched_coverage,
                max_neighbor_distance = max_neighbor_distance,
                min_ref_cpg_match_frac = min_ref_cpg_match_frac,
                docker = methylation_docker,
                runtime_attr_override = runtime_attr_process_sample
        }
    }

    Int n_samples = length(methylation_beds)
    Int n_batches = (n_samples + samples_per_batch - 1) / samples_per_batch

    scatter (b in range(n_batches)) {
        scatter (j in range(samples_per_batch)) {
            Int k = b * samples_per_batch + j
            if (k < n_samples) {
                File batch_member = ProcessSample.sample_npz[k]
            }
        }

        call AccumulateBatch {
            input:
                sample_npzs = select_all(batch_member),
                cpg_index = BuildCpgIndex.cpg_index,
                min_coverage = min_coverage,
                pca_site_stride = pca_site_stride,
                prefix = "~{prefix}.batch~{b}",
                docker = methylation_docker,
                runtime_attr_override = runtime_attr_accumulate
        }
    }

    call Helpers.ConcatTsvs as ConcatSampleQc {
        input:
            tsvs = ProcessSample.qc_tsv,
            prefix = "~{prefix}.sample_qc.unflagged",
            preserve_header = true,
            docker = methylation_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call Helpers.ConcatTsvs as ConcatSampleHists {
        input:
            tsvs = ProcessSample.hists_tsv,
            prefix = "~{prefix}.sample_hists",
            preserve_header = true,
            docker = methylation_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call FinalizeSiteStats {
        input:
            partials = AccumulateBatch.partial_npz,
            cpg_index = BuildCpgIndex.cpg_index,
            sample_qc = ConcatSampleQc.concatenated_tsv,
            platform = platform,
            min_coverage = min_coverage,
            matched_coverage = matched_coverage,
            min_samples_per_site = min_samples_per_site,
            outlier_mad = outlier_mad,
            prefix = prefix,
            docker = methylation_docker,
            runtime_attr_override = runtime_attr_finalize
    }

    call PlotQC {
        input:
            sample_qc = FinalizeSiteStats.sample_qc_tsv,
            sample_hists = ConcatSampleHists.concatenated_tsv,
            site_sd_hist = FinalizeSiteStats.site_sd_hist,
            site_sd_by_mean = FinalizeSiteStats.site_sd_by_mean,
            site_mean_sd_grid = FinalizeSiteStats.site_mean_sd_grid,
            site_excess_var_hist = FinalizeSiteStats.site_excess_var_hist,
            pca_coords = FinalizeSiteStats.pca_coords,
            pca_variance = FinalizeSiteStats.pca_variance,
            prefix = prefix,
            docker = methylation_docker,
            runtime_attr_override = runtime_attr_plot
    }

    output {
        File cpg_index = BuildCpgIndex.cpg_index
        File sample_qc_tsv = FinalizeSiteStats.sample_qc_tsv
        File sample_hists_tsv = ConcatSampleHists.concatenated_tsv
        File cohort_summary = FinalizeSiteStats.cohort_summary
        File per_site_stats = FinalizeSiteStats.per_site_stats
        File per_site_stats_idx = FinalizeSiteStats.per_site_stats_idx
        File site_sd_hist = FinalizeSiteStats.site_sd_hist
        File site_sd_by_mean = FinalizeSiteStats.site_sd_by_mean
        File site_mean_sd_grid = FinalizeSiteStats.site_mean_sd_grid
        File site_excess_var_hist = FinalizeSiteStats.site_excess_var_hist
        File pca_matrix = FinalizeSiteStats.pca_matrix
        File pca_coords = FinalizeSiteStats.pca_coords
        File pca_variance = FinalizeSiteStats.pca_variance
        Array[File] figures = PlotQC.figures
        File figures_tar = PlotQC.figures_tar
    }
}

task BuildCpgIndex {
    input {
        File ref_fa
        File ref_fai
        File? cpg_islands_bed
        Array[String] contigs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/methylation/methylation_qc.py build-index \
            --ref-fa ~{ref_fa} \
            --ref-fai ~{ref_fai} \
            --contigs ~{sep=' ' contigs} \
            ~{"--cpg-islands-bed " + cpg_islands_bed} \
            --output ~{prefix}.cpg_index.npz
    >>>

    output {
        File cpg_index = "~{prefix}.cpg_index.npz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(ref_fa, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task ProcessSample {
    input {
        File bed
        String sample_id
        String platform
        String pacbio_beta_source
        File cpg_index
        Int min_coverage
        Int matched_coverage
        Int max_neighbor_distance
        Float min_ref_cpg_match_frac
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/methylation/methylation_qc.py sample \
            --bed ~{bed} \
            --sample-id '~{sample_id}' \
            --platform ~{platform} \
            --pacbio-beta-source ~{pacbio_beta_source} \
            --cpg-index ~{cpg_index} \
            --prefix '~{sample_id}' \
            --min-coverage ~{min_coverage} \
            --matched-coverage ~{matched_coverage} \
            --max-neighbor-distance ~{max_neighbor_distance} \
            --min-ref-cpg-match-frac ~{min_ref_cpg_match_frac}
    >>>

    output {
        File sample_npz = "~{sample_id}.methyl.npz"
        File qc_tsv = "~{sample_id}.qc.tsv"
        File hists_tsv = "~{sample_id}.hists.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(bed, "GB") + size(cpg_index, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 3,
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

task AccumulateBatch {
    input {
        Array[File] sample_npzs
        File cpg_index
        Int min_coverage
        Int pca_site_stride
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/methylation/methylation_qc.py accumulate \
            --cpg-index ~{cpg_index} \
            --min-coverage ~{min_coverage} \
            --pca-site-stride ~{pca_site_stride} \
            --output ~{prefix}.partial.npz \
            --samples ~{sep=' ' sample_npzs}
    >>>

    output {
        File partial_npz = "~{prefix}.partial.npz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 12,
        disk_gb: ceil(size(sample_npzs, "GB") + size(cpg_index, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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

task FinalizeSiteStats {
    input {
        Array[File] partials
        File cpg_index
        File sample_qc
        String platform
        Int min_coverage
        Int matched_coverage
        Int min_samples_per_site
        Float outlier_mad
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/methylation/methylation_qc.py finalize \
            --cpg-index ~{cpg_index} \
            --sample-qc ~{sample_qc} \
            --platform ~{platform} \
            --prefix ~{prefix} \
            --min-coverage ~{min_coverage} \
            --matched-coverage ~{matched_coverage} \
            --min-samples-per-site ~{min_samples_per_site} \
            --outlier-mad ~{outlier_mad} \
            --partials ~{sep=' ' partials}
    >>>

    output {
        File sample_qc_tsv = "~{prefix}.sample_qc.tsv"
        File cohort_summary = "~{prefix}.cohort_summary.tsv"
        File per_site_stats = "~{prefix}.per_site_stats.tsv.gz"
        File per_site_stats_idx = "~{prefix}.per_site_stats.tsv.gz.tbi"
        File site_sd_hist = "~{prefix}.site_sd_hist.tsv"
        File site_sd_by_mean = "~{prefix}.site_sd_by_mean.tsv"
        File site_mean_sd_grid = "~{prefix}.site_mean_sd_grid.tsv"
        File site_excess_var_hist = "~{prefix}.site_excess_var_hist.tsv"
        File pca_matrix = "~{prefix}.pca_matrix.npz"
        File pca_coords = "~{prefix}.pca_coords.tsv"
        File pca_variance = "~{prefix}.pca_variance.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 32,
        disk_gb: ceil(size(partials, "GB") + size(cpg_index, "GB")) + 50,
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

task PlotQC {
    input {
        File sample_qc
        File sample_hists
        File site_sd_hist
        File site_sd_by_mean
        File site_mean_sd_grid
        File site_excess_var_hist
        File pca_coords
        File pca_variance
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail
        python3 /opt/gnomad-lr/scripts/methylation/methylation_qc.py plot \
            --sample-qc ~{sample_qc} \
            --sample-hists ~{sample_hists} \
            --site-sd-hist ~{site_sd_hist} \
            --site-sd-by-mean ~{site_sd_by_mean} \
            --site-mean-sd-grid ~{site_mean_sd_grid} \
            --site-excess-var-hist ~{site_excess_var_hist} \
            --pca-coords ~{pca_coords} \
            --pca-variance ~{pca_variance} \
            --prefix ~{prefix} \
            --outdir figures
        tar -czf ~{prefix}.figures.tar.gz figures
    >>>

    output {
        Array[File] figures = glob("figures/*")
        File figures_tar = "~{prefix}.figures.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(sample_hists, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
