version 1.0

import "../utils/Structs.wdl"

workflow CreateDepthProfile {
    input {
        Array[String] sample_ids
        Array[File] mosdepth_bed_files
        Array[File] mosdepth_bed_idx
        String contig
        Int window_start
        Int window_end
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_extract
        RuntimeAttr? runtime_attr_build_matrix
    }

    scatter (i in range(length(sample_ids))) {
        call ExtractSampleDepth {
            input:
                mosdepth_bed = mosdepth_bed_files[i],
                mosdepth_bed_idx = mosdepth_bed_idx[i],
                sample_id = sample_ids[i],
                contig = contig,
                window_start = window_start,
                window_end = window_end,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_extract
        }
    }

    call BuildDepthMatrix {
        input:
            per_sample_depths = ExtractSampleDepth.sample_depth,
            sample_ids = sample_ids,
            contig = contig,
            window_start = window_start,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_build_matrix
    }

    output {
        File region_depth_profile = BuildDepthMatrix.depth_profile
    }
}

task ExtractSampleDepth {
    input {
        File mosdepth_bed
        File mosdepth_bed_idx
        String sample_id
        String contig
        Int window_start
        Int window_end
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import subprocess

mosdepth_bed = "~{mosdepth_bed}"
contig = "~{contig}"
window_start = ~{window_start}
window_end = ~{window_end}
output_file = "~{sample_id}.depth.txt"

proc = subprocess.Popen(
    ['tabix', mosdepth_bed, f'{contig}:{window_start + 1}-{window_end}'],
    stdout=subprocess.PIPE, text=True
)

pos = window_start
with open(output_file, 'w') as out:
    for line in proc.stdout:
        parts = line.strip().split('\t')
        s = int(parts[1])
        e = int(parts[2])
        cov = parts[3].strip()

        s_clamped = max(s, window_start)
        e_clamped = min(e, window_end)
        if s_clamped >= e_clamped:
            continue

        while pos < s_clamped:
            out.write("0\n")
            pos += 1

        while pos < e_clamped:
            out.write(f"{cov}\n")
            pos += 1

    while pos < window_end:
        out.write("0\n")
        pos += 1

proc.wait()
CODE
    >>>

    output {
        File sample_depth = "~{sample_id}.depth.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(mosdepth_bed, "GB")) + 10,
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

task BuildDepthMatrix {
    input {
        Array[File] per_sample_depths
        Array[String] sample_ids
        String contig
        Int window_start
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
sample_ids = "~{sep=',' sample_ids}".split(',')
depth_files = "~{sep=',' per_sample_depths}".split(',')
contig = "~{contig}"
window_start = ~{window_start}
output_file = "~{prefix}.depth_profile.bed"

sample_to_file = {}
for sid, fpath in zip(sample_ids, depth_files):
    sample_to_file[sid] = fpath

handles = [open(sample_to_file[sid], 'r') for sid in sample_ids]

with open(output_file, 'w') as out:
    out.write("#CHROM\tSTART\tEND\t" + "\t".join(sample_ids) + "\n")

    pos = window_start
    while True:
        depths = [h.readline() for h in handles]
        if not depths[0]:
            break
        cols = "\t".join(d.strip() for d in depths)
        out.write(f"{contig}\t{pos}\t{pos + 1}\t{cols}\n")
        pos += 1

for h in handles:
    h.close()
CODE

        bgzip ~{prefix}.depth_profile.bed
    >>>

    output {
        File depth_profile = "~{prefix}.depth_profile.bed.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(per_sample_depths, "GB")) + 10,
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
