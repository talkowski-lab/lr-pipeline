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
# (passed in as `analysis_script`).
#
#   1. RegionStats (scattered per contig): every region with >= min_cpg CpG
#      sites in the methylation table within [start, end) is summarised per
#      haplotype: number of called CpGs, mean, median and SD (SD needs >= 2
#      called CpGs) across the region's CpG sites.
#   2. ConcatTables: per-contig outputs are concatenated genome-wide.
#
# Outputs share the key `region_index` (0-based row of the region in
# `regions_bed`, header excluded):
#   regions   original region columns + region_n_cpg + region_n_haplotypes_called
#   n_called / mean / median / sd   region x haplotype matrices (NA = not computable)

workflow RegionHaplotypeMethylation {
    input {
        Array[File] haplotype_methylation_beds
        File regions_bed
        String prefix

        Int min_cpg = 1

        File analysis_script

        String python_docker

        RuntimeAttr? runtime_attr_region_stats
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
                analysis_script = analysis_script,
                prefix = "~{prefix}.~{contig}",
                docker = python_docker,
                runtime_attr_override = runtime_attr_region_stats
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

    output {
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
            --prefix ~{prefix}
    >>>

    output {
        File regions = "~{prefix}.regions.tsv.gz"
        File n_called = "~{prefix}.n_called.tsv.gz"
        File mean = "~{prefix}.mean.tsv.gz"
        File median = "~{prefix}.median.tsv.gz"
        File sd = "~{prefix}.sd.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: ceil(10 * size(bed, "GB") + 4 * size(regions_bed, "GB")) + 8,
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
