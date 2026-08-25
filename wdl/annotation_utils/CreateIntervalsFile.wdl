version 1.0

import "../utils/Structs.wdl"

workflow CreateIntervalsFile {
    input {
        File ref_fai
        Array[String] contigs
        String prefix

        Int bin_size

        String utils_docker

        RuntimeAttr? runtime_attr_create_intervals
    }

    File contigs_file = write_lines(contigs)

    call CreateIntervals {
        input:
            ref_fai = ref_fai,
            contigs_file = contigs_file,
            prefix = prefix,
            bin_size = bin_size,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_intervals
    }

    output {
        File intervals = CreateIntervals.intervals
    }
}

task CreateIntervals {
    input {
        File ref_fai
        File contigs_file
        String prefix
        Int bin_size
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
from pathlib import Path

bin_size = ~{bin_size}
if bin_size <= 0:
    raise ValueError("bin_size must be greater than zero")

contigs = Path("~{contigs_file}").read_text().splitlines()
lengths = {}
with open("~{ref_fai}") as fai:
    for line in fai:
        fields = line.rstrip().split("\t")
        lengths[fields[0]] = int(fields[1])

missing_contigs = [contig for contig in contigs if contig not in lengths]
if missing_contigs:
    raise ValueError(f"Contigs absent from reference index: {', '.join(missing_contigs)}")

with open("~{prefix}.intervals", "w") as intervals:
    for contig in contigs:
        for end in range(bin_size, lengths[contig] + 1, bin_size):
            intervals.write(f"{contig}:{end - bin_size + 1}-{end}\n")
CODE
    >>>

    output {
        File intervals = "~{prefix}.intervals"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: ceil(size(ref_fai, "GB")) + 10,
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
