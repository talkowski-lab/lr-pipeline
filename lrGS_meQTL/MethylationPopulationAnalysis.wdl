version 1.0

import "../wdl/utils/Structs.wdl"

# Population-level analysis of long-read CpG methylation tables (per-sample
# and per-haplotype wide-format beds, one file per contig; values are percent
# methylation, `.` = missing). All logic lives in
# methylation_population_analysis.py (passed in as `analysis_script`).
#
# Part 1 (per-sample tables):
#   1. SampleContigStats (scattered per contig): per-sample call rate, mean,
#      median, fraction of low/intermediate/high CpGs, correlation with the
#      cohort per-CpG median, and a 0.1%-resolution value histogram.
#   2. SampleQC: contig x sample mean / median / call-rate tables; autosomal
#      genome-wide metrics; samples are excluded when any autosomal metric
#      (mean, median, fraction intermediate) is a robust-z outlier
#      (|z| >= sample_qc_z), the correlation with the cohort is low
#      (z <= -sample_qc_z), or call rate < sample_qc_min_call_rate. Sex is
#      inferred from chrY call rate.
#   3. VariableRegions (scattered per contig, QC-passing samples only; males
#      only on chrY): CpGs called in >= variable_min_call_rate of samples are
#      grouped into non-overlapping tiles of tile_cpgs consecutive CpGs; the
#      across-sample SD of each tile's mean is compared with tiles of similar
#      mean methylation (robust z within 5%-wide bins of the mean); tiles with
#      z >= variable_z are merged into regions. Within each region every
#      sample's CpG-level difference from the per-CpG cohort median is tested
#      (Wilcoxon signed-rank, BH across samples): q < variable_fdr and mean
#      shift >= +variable_min_delta = hyper (blue in plots), <=
#      -variable_min_delta = hypo (red in plots).
# Part 2 (per-haplotype tables):
#   4. HaplotypeContigStats (scattered per contig): per-sample hap1/hap2 call
#      counts and the noise (robust SD) of the per-tile hap1 - hap2 delta.
#   5. HaplotypeQC: excludes Part 1 exclusions, samples with few CpGs called
#      on both haplotypes, and samples with outlying delta noise or hap1/hap2
#      mean imbalance.
#   6. AlleleSpecificRegions (scattered per contig): per sample, tiles with
#      |hap1 - hap2| >= max(asm_min_delta, asm_sigma_k * sample noise) are
#      merged (same sign) into ASM regions, tested with a Wilcoxon signed-rank
#      test over CpGs (BH across all regions of the contig, q < asm_fdr), then
#      unioned across samples into population ASM regions.
#
# `methylation_beds` and `haplotype_methylation_beds` are per-contig files;
# the contig is the second dot-delimited field of the file name (e.g.
# hprc_methylated.chr22.combined.bed.gz -> chr22).

workflow MethylationPopulationAnalysis {
    input {
        Array[File] methylation_beds
        Array[File] haplotype_methylation_beds
        String prefix

        Float sample_qc_z = 5.0
        Float sample_qc_min_call_rate = 0.5
        Float variable_min_call_rate = 0.8
        Int tile_cpgs = 10
        Float variable_z = 4.0
        Float variable_min_delta = 20.0
        Float variable_fdr = 0.05
        Int variable_max_plots = 1000
        Float haplotype_qc_min_both_frac = 0.25
        Float asm_min_delta = 20.0
        Float asm_sigma_k = 4.0
        Float asm_fdr = 0.05
        Int asm_max_plots = 200

        File analysis_script
        File font_ttf

        String python_docker = "quay.io/jupyter/scipy-notebook@sha256:211986c06e20f36a2adb9342315f61cfa12df7124d3ad97f618b23dc7a4e97b9"

        RuntimeAttr? runtime_attr_sample_contig_stats
        RuntimeAttr? runtime_attr_sample_qc
        RuntimeAttr? runtime_attr_variable_regions
        RuntimeAttr? runtime_attr_haplotype_contig_stats
        RuntimeAttr? runtime_attr_haplotype_qc
        RuntimeAttr? runtime_attr_allele_specific_regions
        RuntimeAttr? runtime_attr_concat
    }

    scatter (bed in methylation_beds) {
        String sample_contig = sub(sub(basename(bed), "^[^.]*\\.", ""), "\\..*$", "")

        call SampleContigStats {
            input:
                bed = bed,
                contig = sample_contig,
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{sample_contig}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_sample_contig_stats
        }
    }

    call SampleQC {
        input:
            stats = SampleContigStats.stats,
            hists = SampleContigStats.hist,
            z = sample_qc_z,
            min_call_rate = sample_qc_min_call_rate,
            analysis_script = analysis_script,
            font_ttf = font_ttf,
            prefix = prefix,
            docker = python_docker,
            runtime_attr_override = runtime_attr_sample_qc
    }

    scatter (i in range(length(methylation_beds))) {
        call VariableRegions {
            input:
                bed = methylation_beds[i],
                contig = sample_contig[i],
                sample_qc = SampleQC.sample_qc,
                min_call_rate = variable_min_call_rate,
                tile_cpgs = tile_cpgs,
                z_var = variable_z,
                min_delta = variable_min_delta,
                fdr = variable_fdr,
                max_plots = variable_max_plots,
                analysis_script = analysis_script,
                font_ttf = font_ttf,
                prefix = "~{prefix}.~{sample_contig[i]}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_variable_regions
        }
    }

    call ConcatTables as ConcatVariableRegions {
        input:
            tables = VariableRegions.regions_bed,
            analysis_script = analysis_script,
            out_name = "~{prefix}.variable_regions.bed",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatVariableRegionCalls {
        input:
            tables = VariableRegions.sample_calls,
            analysis_script = analysis_script,
            out_name = "~{prefix}.variable_region_sample_calls.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    scatter (bed in haplotype_methylation_beds) {
        String haplotype_contig = sub(sub(basename(bed), "^[^.]*\\.", ""), "\\..*$", "")

        call HaplotypeContigStats {
            input:
                bed = bed,
                contig = haplotype_contig,
                tile_cpgs = tile_cpgs,
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{haplotype_contig}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_haplotype_contig_stats
        }
    }

    call HaplotypeQC {
        input:
            stats = HaplotypeContigStats.stats,
            sample_qc = SampleQC.sample_qc,
            z = sample_qc_z,
            min_both_frac = haplotype_qc_min_both_frac,
            analysis_script = analysis_script,
            font_ttf = font_ttf,
            prefix = prefix,
            docker = python_docker,
            runtime_attr_override = runtime_attr_haplotype_qc
    }

    scatter (i in range(length(haplotype_methylation_beds))) {
        call AlleleSpecificRegions {
            input:
                bed = haplotype_methylation_beds[i],
                contig = haplotype_contig[i],
                haplotype_qc = HaplotypeQC.haplotype_qc,
                tile_cpgs = tile_cpgs,
                min_delta = asm_min_delta,
                sigma_k = asm_sigma_k,
                fdr = asm_fdr,
                max_plots = asm_max_plots,
                analysis_script = analysis_script,
                font_ttf = font_ttf,
                prefix = "~{prefix}.~{haplotype_contig[i]}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_allele_specific_regions
        }
    }

    call ConcatTables as ConcatAsmPerSample {
        input:
            tables = AlleleSpecificRegions.asm_per_sample,
            analysis_script = analysis_script,
            out_name = "~{prefix}.asm_regions_per_sample.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatAsmPopulation {
        input:
            tables = AlleleSpecificRegions.asm_population,
            analysis_script = analysis_script,
            out_name = "~{prefix}.asm_regions_population.bed",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File mean_methylation_by_contig = SampleQC.mean_by_contig
        File median_methylation_by_contig = SampleQC.median_by_contig
        File call_rate_by_contig = SampleQC.call_rate_by_contig
        File mean_methylation_robust_z_by_contig = SampleQC.mean_robust_z_by_contig
        File per_contig_sample_stats = SampleQC.per_contig_stats
        File sample_qc = SampleQC.sample_qc
        File excluded_samples = SampleQC.excluded_samples
        File sample_qc_pdf = SampleQC.qc_pdf
        File variable_regions_bed = ConcatVariableRegions.merged
        File variable_region_sample_calls = ConcatVariableRegionCalls.merged
        Array[File] variable_region_tile_stats = VariableRegions.tile_stats
        Array[File] variable_region_pdfs = VariableRegions.regions_pdf
        File haplotype_sample_qc = HaplotypeQC.haplotype_qc
        File haplotype_excluded_samples = HaplotypeQC.excluded_samples
        File haplotype_n_both_by_contig = HaplotypeQC.n_both_by_contig
        File haplotype_sample_qc_pdf = HaplotypeQC.qc_pdf
        File asm_regions_per_sample = ConcatAsmPerSample.merged
        File asm_regions_population_bed = ConcatAsmPopulation.merged
        Array[File] asm_region_pdfs = AlleleSpecificRegions.asm_pdf
    }
}

task SampleContigStats {
    input {
        File bed
        String contig
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} contig-stats \
            --bed ~{bed} \
            --contig ~{contig} \
            --prefix ~{prefix}
    >>>

    output {
        File stats = "~{prefix}.sample_stats.tsv"
        File hist = "~{prefix}.sample_hist.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: ceil(12 * size(bed, "GB")) + 8,
        disk_gb: 2 * ceil(size(bed, "GB")) + 10,
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

task SampleQC {
    input {
        Array[File] stats
        Array[File] hists
        Float z
        Float min_call_rate
        File analysis_script
        File font_ttf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} sample-qc \
            --stats ~{sep=" " stats} \
            --hists ~{sep=" " hists} \
            --z ~{z} \
            --min-call-rate ~{min_call_rate} \
            --font-ttf ~{font_ttf} \
            --prefix ~{prefix}
    >>>

    output {
        File mean_by_contig = "~{prefix}.mean_methylation_by_contig.tsv"
        File median_by_contig = "~{prefix}.median_methylation_by_contig.tsv"
        File call_rate_by_contig = "~{prefix}.call_rate_by_contig.tsv"
        File mean_robust_z_by_contig = "~{prefix}.mean_methylation_robust_z_by_contig.tsv"
        File per_contig_stats = "~{prefix}.per_contig_sample_stats.tsv"
        File sample_qc = "~{prefix}.sample_qc.tsv"
        File excluded_samples = "~{prefix}.excluded_samples.txt"
        File qc_pdf = "~{prefix}.sample_qc.pdf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(stats, "GB") + size(hists, "GB")) + 10,
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

task VariableRegions {
    input {
        File bed
        String contig
        File sample_qc
        Float min_call_rate
        Int tile_cpgs
        Float z_var
        Float min_delta
        Float fdr
        Int max_plots
        File analysis_script
        File font_ttf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} variable-regions \
            --bed ~{bed} \
            --contig ~{contig} \
            --sample-qc ~{sample_qc} \
            --min-call-rate ~{min_call_rate} \
            --tile-cpgs ~{tile_cpgs} \
            --z-var ~{z_var} \
            --min-delta ~{min_delta} \
            --fdr ~{fdr} \
            --max-plots ~{max_plots} \
            --font-ttf ~{font_ttf} \
            --prefix ~{prefix}
    >>>

    output {
        File regions_bed = "~{prefix}.variable_regions.bed"
        File sample_calls = "~{prefix}.variable_region_sample_calls.tsv.gz"
        File tile_stats = "~{prefix}.tile_stats.tsv.gz"
        File regions_pdf = "~{prefix}.variable_regions.pdf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: ceil(16 * size(bed, "GB")) + 8,
        disk_gb: 3 * ceil(size(bed, "GB")) + 20,
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

task HaplotypeContigStats {
    input {
        File bed
        String contig
        Int tile_cpgs
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} haplotype-stats \
            --bed ~{bed} \
            --contig ~{contig} \
            --tile-cpgs ~{tile_cpgs} \
            --prefix ~{prefix}
    >>>

    output {
        File stats = "~{prefix}.haplotype_stats.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: ceil(10 * size(bed, "GB")) + 8,
        disk_gb: 2 * ceil(size(bed, "GB")) + 10,
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

task HaplotypeQC {
    input {
        Array[File] stats
        File sample_qc
        Float z
        Float min_both_frac
        File analysis_script
        File font_ttf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} haplotype-qc \
            --stats ~{sep=" " stats} \
            --sample-qc ~{sample_qc} \
            --z ~{z} \
            --min-both-frac ~{min_both_frac} \
            --font-ttf ~{font_ttf} \
            --prefix ~{prefix}
    >>>

    output {
        File haplotype_qc = "~{prefix}.haplotype_sample_qc.tsv"
        File excluded_samples = "~{prefix}.haplotype_excluded_samples.txt"
        File n_both_by_contig = "~{prefix}.n_both_haplotypes_by_contig.tsv"
        File qc_pdf = "~{prefix}.haplotype_sample_qc.pdf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task AlleleSpecificRegions {
    input {
        File bed
        String contig
        File haplotype_qc
        Int tile_cpgs
        Float min_delta
        Float sigma_k
        Float fdr
        Int max_plots
        File analysis_script
        File font_ttf
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} allele-specific \
            --bed ~{bed} \
            --contig ~{contig} \
            --haplotype-qc ~{haplotype_qc} \
            --tile-cpgs ~{tile_cpgs} \
            --min-delta ~{min_delta} \
            --sigma-k ~{sigma_k} \
            --fdr ~{fdr} \
            --max-plots ~{max_plots} \
            --font-ttf ~{font_ttf} \
            --prefix ~{prefix}
    >>>

    output {
        File asm_per_sample = "~{prefix}.asm_regions_per_sample.tsv.gz"
        File asm_population = "~{prefix}.asm_regions_population.bed"
        File asm_pdf = "~{prefix}.asm_regions.pdf"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: ceil(10 * size(bed, "GB")) + 8,
        disk_gb: 2 * ceil(size(bed, "GB")) + 20,
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

task ConcatTables {
    input {
        Array[File] tables
        File analysis_script
        String out_name
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} concat \
            --inputs ~{sep=" " tables} \
            --out ~{out_name}
    >>>

    output {
        File merged = "~{out_name}"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(tables, "GB")) + 10,
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
