version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CreateSampleReadCounts {
    input {
        Array[File] mosdepth_bed_files
        Array[File] mosdepth_bed_idx
        File ref_dict
        Array[String] contigs
        String prefix

        Int bin_size
        String sample_id

        String utils_docker

        RuntimeAttr? runtime_attr_bin
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_sort
    }

    scatter (i in range(length(contigs))) {
        call BinMosDepthCounts {
            input:
                mosdepth_bed = mosdepth_bed_files[i],
                mosdepth_bed_idx = mosdepth_bed_idx[i],
                contig = contigs[i],
                bin_size = bin_size,
                prefix = "~{prefix}.~{contigs[i]}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_bin
        }
    }

    call MergeBinnedCounts {
        input:
            binned_counts = BinMosDepthCounts.binned_counts,
            binned_counts_idx = BinMosDepthCounts.binned_counts_idx,
            ref_dict = ref_dict,
            sample_id = sample_id,
            prefix = "~{prefix}.counts.unsorted",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge
    }

    call Helpers.SortReadCounts {
        input:
            read_counts = MergeBinnedCounts.merged_counts,
            prefix = "~{prefix}.counts",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_sort
    }

    output {
        File binned_read_counts = SortReadCounts.sorted_read_counts
    }
}

task BinMosDepthCounts {
    input {
        File mosdepth_bed
        File mosdepth_bed_idx
        String contig
        Int bin_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip
import sys
from statistics import median

bed_file = "~{mosdepth_bed}"
contig = "~{contig}"
bin_size = ~{bin_size}
output_file = "~{prefix}.tsv"

bins = {}

with gzip.open(bed_file, 'rt') as f:
    for line in f:
        parts = line.strip().split('\t')
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        coverage = int(float(parts[3]))
        
        start_bin = (start // bin_size) * bin_size
        end_bin = ((end - 1) // bin_size) * bin_size
        
        for bin_start in range(start_bin, end_bin + 1, bin_size):
            bin_end = bin_start + bin_size
            overlap_start = max(start, bin_start)
            overlap_end = min(end, bin_end)
            overlap_length = overlap_end - overlap_start
            
            if bin_start not in bins:
                bins[bin_start] = []
            bins[bin_start].extend([coverage] * overlap_length)

with open(output_file, 'w') as out:
    for bin_start in sorted(bins.keys()):
        if len(bins[bin_start]) == bin_size:
            bin_end = bin_start + bin_size
            median_coverage = int(median(bins[bin_start]))
            out.write(f"{contig}\t{bin_start + 1}\t{bin_end}\t{median_coverage}\n")
CODE

        bgzip ~{prefix}.tsv
        tabix -s 1 -b 2 -e 3 ~{prefix}.tsv.gz
    >>>

    output {
        File binned_counts = "~{prefix}.tsv.gz"
        File binned_counts_idx = "~{prefix}.tsv.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(mosdepth_bed, "GB")) + 10,
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

task MergeBinnedCounts {
    input {
        Array[File] binned_counts
        Array[File] binned_counts_idx
        File ref_dict
        String sample_id
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        grep "^@" ~{ref_dict} > ~{prefix}.tsv
        echo -e "@RG\tID:GATKCopyNumber\tSM:~{sample_id}" >> ~{prefix}.tsv
        echo -e "CONTIG\tSTART\tEND\tCOUNT" >> ~{prefix}.tsv
        zcat ~{sep=' ' binned_counts} >> ~{prefix}.tsv

        bgzip ~{prefix}.tsv
    >>>

    output {
        File merged_counts = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(binned_counts, "GB")) + 10,
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
