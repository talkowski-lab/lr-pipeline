version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow CreateCohortCoverageSummary {
    input {
        Array[File] mosdepth_bed_files
        Array[File] mosdepth_bed_idx
        File ref_fai
        Array[String] contigs
        String prefix

        Int window_size
        Int bin_size
        Array[Int] thresholds

        String utils_docker

        RuntimeAttr? runtime_attr_make_windows
        RuntimeAttr? runtime_attr_compute_coverage
        RuntimeAttr? runtime_attr_merge_coverages
    }

    call Helpers.MakeWindows {
        input:
            ref_fai = ref_fai,
            contigs = contigs,
            window_size = window_size,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_make_windows
    }

    scatter (region in MakeWindows.regions) {
        call ComputeBinnedCoverage {
            input:
                mosdepth_beds = mosdepth_bed_files,
                mosdepth_bed_idx = mosdepth_bed_idx,
                region = region,
                bin_size = bin_size,
                thresholds = thresholds,
                prefix = "~{prefix}." + sub(region, ":", "_"),
                docker = utils_docker,
                runtime_attr_override = runtime_attr_compute_coverage
        }
    }

    call ConcatenateCoverages {
        input:
            formatted_coverages = ComputeBinnedCoverage.binned_coverage,
            thresholds = thresholds,
            prefix = "~{prefix}.coverage",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_coverages
    }

    output {
        File binned_coverage_tsv = ConcatenateCoverages.merged_coverage
    }
}

task ComputeBinnedCoverage {
    input {
        Array[File] mosdepth_beds
        Array[File] mosdepth_bed_idx
        String region
        Int bin_size
        Array[Int] thresholds
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import gzip
import subprocess
import sys

mosdepth_files = "~{sep=',' mosdepth_beds}".split(',')
bin_size = ~{bin_size}
thresholds = [~{sep=',' thresholds}]
region = "~{region}"
output_file = "~{prefix}.tsv"
num_samples = len(mosdepth_files)

# Parse boundaries
if ":" in region:
    contig, span = region.split(":")
    reg_start, reg_end = map(int, span.split("-"))
else:
    contig = region
    reg_start = 0
    reg_end = float('inf')

iterators = []
current_rows = []

for bed_file in mosdepth_files:
    proc = subprocess.Popen(
        ['tabix', bed_file, region], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.DEVNULL,
        text=True
    )
    iterators.append(proc.stdout)
    
    # Fast-forward to first valid overlapping interval
    line = proc.stdout.readline()
    while line:
        parts = line.strip().split('\t')

        # Clamp coordinates to boundaries
        s = max(int(parts[1]), reg_start)
        e = min(int(parts[2]), reg_end)
        if s < e:
            current_rows.append((s, e, float(parts[3])))
            break
        line = proc.stdout.readline()
    else:
        current_rows.append(None)

with gzip.open(output_file + '.gz', 'wt') as out:
    active_starts = [r[0] for r in current_rows if r]
    current_pos = max((min(active_starts) // bin_size) * bin_size, reg_start)

    while True:
        # Advance any iterator that has fallen behind current_pos
        for i in range(num_samples):
            while current_rows[i] is not None and current_rows[i][1] <= current_pos:
                line = iterators[i].readline()
                while line:
                    parts = line.strip().split('\t')
                    s = max(int(parts[1]), reg_start)
                    e = min(int(parts[2]), reg_end)
                    if s < e:
                        current_rows[i] = (s, e, float(parts[3]))
                        break
                    line = iterators[i].readline()
                else:
                    current_rows[i] = None
        
        if any(r is None for r in current_rows):
            break
            
        # Find the end of this block of constant coverage across all samples
        min_end = min((r[1] for r in current_rows))
        
        if min_end > current_pos:
            # Calculate stats once for the whole constant block
            depths = [r[2] for r in current_rows]
            total_dp = sum(depths)
            mean_cov = round(total_dp / num_samples, 6)
            depths_sorted = sorted(depths)
            median_cov = int(depths_sorted[num_samples // 2])
            over_fracs = [round(sum(1 for d in depths if d > t) / num_samples, 6) for t in thresholds]
            stats_str = f"{mean_cov:.6f}\t{median_cov}\t{int(total_dp)}\t" + "\t".join(f"{f:.6f}" for f in over_fracs)
            
            # Output bins matching the bin_size
            while current_pos < min_end:
                next_bin_boundary = ((current_pos // bin_size) + 1) * bin_size
                actual_end = min(next_bin_boundary, min_end)
                out.write(f"{contig}\t{current_pos + 1}\t{actual_end}\t{stats_str}\n")
                current_pos = actual_end

for it in iterators:
    it.close()
CODE
    >>>

    output {
        File binned_coverage = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(mosdepth_beds, "GB")) + 20,
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

task ConcatenateCoverages {
    input {
        Array[File] formatted_coverages
        Array[Int] thresholds
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        zcat ~{sep=' ' formatted_coverages} | bgzip > ~{prefix}.tsv.gz
    >>>

    output {
        File merged_coverage = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 3 * ceil(size(formatted_coverages, "GB")) + 10,
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
