version 1.0

import "../wdl/utils/Structs.wdl"

# Per-haplotype methylation summaries within target regions (e.g. ABC enhancers).
#
# `haplotype_methylation_beds` are per-contig wide tables (#chrom start end
# <sample>_hap1 ... <sample>_hap2 ...; percent methylation, `.` = missing);
# the contig is the second dot-delimited field of the file name (e.g.
# hprc_methylated.chr22.haplotype.bed.gz -> chr22). `regions_bed` is any
# BED-like file (optional `#` header, optional gzip); all of its columns are
# carried through. All logic lives in region_haplotype_methylation.py
# (passed in as `analysis_script`). If `excluded_samples` is given (one sample
# per line, optional tab-separated reason), those samples are removed before
# any statistics, so they are absent from every output and from sample QC.
#
#   1. RegionStats (scattered per contig): every region with >= min_cpg CpG
#      sites in the methylation table within [start, end) is summarised per
#      haplotype: number of called CpGs, mean, median and SD (SD needs >= 2
#      called CpGs) across the region's CpG sites. Also writes per-haplotype
#      CpG-site statistics for sample QC.
#   2. SampleQC: autosomal per-haplotype call rate, mean, median, fraction of
#      intermediate CpGs (20-80%) and correlation with the per-CpG cohort
#      median; a sample is excluded if either haplotype is a robust-z outlier
#      (|z| >= qc_z; low side only for the correlation), has a call rate below
#      qc_min_call_rate_frac x the cohort median, or if its autosomal
#      hap1 - hap2 mean is a robust-z outlier.
#   3. RegionOutliers (scattered per contig): MethBat background-comparison
#      rule on the region means of QC-passing samples. A sample is examined in
#      a region only if both haplotypes are called at >= outlier_min_hap_call_frac
#      of the region's CpG sites; otherwise it is left out of the background
#      and of the calls. Per region, combined =
#      mean of the hap1/hap2 means and abs_hap_delta = |hap2 - hap1|; the
#      background is their mean and SD across QC-passing samples (needs >=
#      outlier_min_samples). HyperMethylated / HypoMethylated: z >= outlier_min_z
#      (<= -) and delta from the background mean >= outlier_min_delta (<= -),
#      with the sample itself >= outlier_methylated_min (<= outlier_unmethylated_max);
#      HyperASM / HypoASM: the same on abs_hap_delta, HyperASM also needing
#      abs_hap_delta >= outlier_asm_min_abs_delta.
#   4. ConcatTables: per-contig outputs are concatenated genome-wide.
#
# Outputs share the key `region_index` (0-based row of the region in
# `regions_bed`, header excluded):
#   regions   original region columns + region_n_cpg + region_n_haplotypes_called
#   n_called / mean / median / sd   region x haplotype matrices for all samples (NA = not computable)
#   region_outliers   per region background statistics, counts and sample lists per outlier label
#   outlier_calls     one row per outlier sample x region with its values, deltas and z-scores

workflow RegionHaplotypeMethylation {
    input {
        Array[File] haplotype_methylation_beds
        File regions_bed
        String prefix

        File? excluded_samples

        Int min_cpg = 1
        Float qc_z = 5.0
        Float qc_min_call_rate_frac = 0.5
        Float outlier_min_z = 3.0
        Float outlier_min_delta = 20.0
        Float outlier_methylated_min = 80.0
        Float outlier_unmethylated_max = 20.0
        Float outlier_asm_min_abs_delta = 50.0
        Int outlier_min_samples = 10
        Float outlier_min_hap_call_frac = 0.8

        File analysis_script

        String python_docker

        RuntimeAttr? runtime_attr_region_stats
        RuntimeAttr? runtime_attr_sample_qc
        RuntimeAttr? runtime_attr_region_outliers
        RuntimeAttr? runtime_attr_concat
    }

    scatter (bed in haplotype_methylation_beds) {
        String contig = sub(sub(basename(bed), "^[^.]*\\.", ""), "\\..*$", "")

        call RegionStats {
            input:
                bed = bed,
                regions_bed = regions_bed,
                contig = contig,
                min_cpg = min_cpg,
                excluded_samples = excluded_samples,
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{contig}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_region_stats
        }
    }

    call SampleQC {
        input:
            site_qc_stats = RegionStats.site_qc_stats,
            site_hists = RegionStats.site_hist,
            qc_z = qc_z,
            qc_min_call_rate_frac = qc_min_call_rate_frac,
            analysis_script = analysis_script,
            prefix = prefix,
            docker = python_docker,
            runtime_attr_override = runtime_attr_sample_qc
    }

    scatter (i in range(length(RegionStats.mean))) {
        call RegionOutliers {
            input:
                mean_matrix = RegionStats.mean[i],
                n_called_matrix = RegionStats.n_called[i],
                regions_table = RegionStats.regions[i],
                excluded_samples = SampleQC.excluded_samples,
                min_z = outlier_min_z,
                min_delta = outlier_min_delta,
                methylated_min = outlier_methylated_min,
                unmethylated_max = outlier_unmethylated_max,
                asm_min_abs_delta = outlier_asm_min_abs_delta,
                min_samples = outlier_min_samples,
                min_hap_call_frac = outlier_min_hap_call_frac,
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{contig[i]}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_region_outliers
        }
    }

    call ConcatTables as ConcatRegions {
        input:
            tables = RegionStats.regions,
            analysis_script = analysis_script,
            out_name = "~{prefix}.regions.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatNCalled {
        input:
            tables = RegionStats.n_called,
            analysis_script = analysis_script,
            out_name = "~{prefix}.n_called.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatMean {
        input:
            tables = RegionStats.mean,
            analysis_script = analysis_script,
            out_name = "~{prefix}.mean.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatMedian {
        input:
            tables = RegionStats.median,
            analysis_script = analysis_script,
            out_name = "~{prefix}.median.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatSd {
        input:
            tables = RegionStats.sd,
            analysis_script = analysis_script,
            out_name = "~{prefix}.sd.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatRegionOutliers {
        input:
            tables = RegionOutliers.region_outliers,
            analysis_script = analysis_script,
            out_name = "~{prefix}.region_outliers.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call ConcatTables as ConcatOutlierCalls {
        input:
            tables = RegionOutliers.outlier_calls,
            analysis_script = analysis_script,
            out_name = "~{prefix}.outlier_calls.tsv.gz",
            docker = python_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File sample_qc = SampleQC.sample_qc
        File qc_excluded_samples = SampleQC.excluded_samples
        File region_outliers = ConcatRegionOutliers.merged
        File outlier_calls = ConcatOutlierCalls.merged
        File region_table = ConcatRegions.merged
        File n_called_matrix = ConcatNCalled.merged
        File mean_matrix = ConcatMean.merged
        File median_matrix = ConcatMedian.merged
        File sd_matrix = ConcatSd.merged
    }
}

task RegionStats {
    input {
        File bed
        File regions_bed
        String contig
        Int min_cpg
        File? excluded_samples
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} contig \
            --bed ~{bed} \
            --regions ~{regions_bed} \
            --contig ~{contig} \
            --min-cpg ~{min_cpg} \
            ~{"--excluded-samples " + excluded_samples} \
            --prefix ~{prefix}
    >>>

    output {
        File regions = "~{prefix}.regions.tsv.gz"
        File n_called = "~{prefix}.n_called.tsv.gz"
        File mean = "~{prefix}.mean.tsv.gz"
        File median = "~{prefix}.median.tsv.gz"
        File sd = "~{prefix}.sd.tsv.gz"
        File site_qc_stats = "~{prefix}.site_qc_stats.tsv.gz"
        File site_hist = "~{prefix}.site_hist.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: ceil(14 * size(bed, "GB") + 4 * size(regions_bed, "GB")) + 8,
        disk_gb: 4 * ceil(size(bed, "GB") + size(regions_bed, "GB")) + 20,
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
        Array[File] site_qc_stats
        Array[File] site_hists
        Float qc_z
        Float qc_min_call_rate_frac
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} sample-qc \
            --stats ~{sep=" " site_qc_stats} \
            --hists ~{sep=" " site_hists} \
            --z ~{qc_z} \
            --min-call-rate-frac ~{qc_min_call_rate_frac} \
            --prefix ~{prefix}
    >>>

    output {
        File sample_qc = "~{prefix}.sample_qc.tsv"
        File excluded_samples = "~{prefix}.excluded_samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(site_qc_stats, "GB") + size(site_hists, "GB")) + 10,
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

task RegionOutliers {
    input {
        File mean_matrix
        File n_called_matrix
        File regions_table
        File excluded_samples
        Float min_z
        Float min_delta
        Float methylated_min
        Float unmethylated_max
        Float asm_min_abs_delta
        Int min_samples
        Float min_hap_call_frac
        File analysis_script
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python ~{analysis_script} outliers \
            --mean-matrix ~{mean_matrix} \
            --n-called-matrix ~{n_called_matrix} \
            --regions-table ~{regions_table} \
            --excluded-samples ~{excluded_samples} \
            --min-z ~{min_z} \
            --min-delta ~{min_delta} \
            --methylated-min ~{methylated_min} \
            --unmethylated-max ~{unmethylated_max} \
            --asm-min-abs-delta ~{asm_min_abs_delta} \
            --min-samples ~{min_samples} \
            --min-hap-call-frac ~{min_hap_call_frac} \
            --prefix ~{prefix}
    >>>

    output {
        File region_outliers = "~{prefix}.region_outliers.tsv.gz"
        File outlier_calls = "~{prefix}.outlier_calls.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(40 * size([mean_matrix, n_called_matrix], "GB")) + 4,
        disk_gb: 3 * ceil(size([mean_matrix, n_called_matrix, regions_table], "GB")) + 10,
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
