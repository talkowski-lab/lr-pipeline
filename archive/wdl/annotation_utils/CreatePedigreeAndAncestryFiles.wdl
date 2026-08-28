version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CreatePedigreeAndAncestryFiles {
    input {
        Array[String] sample_ids
        Array[String] sexes
        String prefix = "cohort"

        String utils_docker

        RuntimeAttr? runtime_attr_override
    }

    call CreatePedFile {
        input:
            sample_ids = sample_ids,
            sexes = sexes,
            prefix = "~{prefix}.ped",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_override
    }

    call CreateAncestryFile {
        input:
            sample_ids = sample_ids,
            prefix = "~{prefix}.ancestry",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File ped = CreatePedFile.ped_file
        File ancestry = CreateAncestryFile.ancestry_file
    }
}

task CreatePedFile {
    input {
        Array[String] sample_ids
        Array[String] sexes
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import sys

sample_ids = ["~{sep='", "' sample_ids}"]
sexes = ["~{sep='", "' sexes}"]

if len(sample_ids) != len(sexes):
    print("Error: sample_ids and sexes must have the same length", file=sys.stderr)
    sys.exit(1)

with open("~{prefix}.ped", "w") as f:
    for sample_id, sex in zip(sample_ids, sexes):
        if sex == "M":
            sex_code = "1"
        elif sex == "F":
            sex_code = "2"
        else:
            sex_code = "0"
        f.write(f"{sample_id}\t{sample_id}\t0\t0\t{sex_code}\t0\n")
CODE
    >>>

    output {
        File ped_file = "~{prefix}.ped"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task CreateAncestryFile {
    input {
        Array[String] sample_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
sample_ids = ["~{sep='", "' sample_ids}"]

with open("~{prefix}.ancestry.tsv", "w") as f:
    for sample_id in sample_ids:
        f.write(f"{sample_id}\tafr\n")
CODE
    >>>

    output {
        File ancestry_file = "~{prefix}.ancestry.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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
