version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow CreateCohortMetadata {
    input {
        File ped_file
        File ancestry_file
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_override
    }

    call CreateMetadata {
        input:
            ped_file = ped_file,
            ancestry_file = ancestry_file,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File metadata = CreateMetadata.metadata_file
    }
}

task CreateMetadata {
    input {
        File ped_file
        File ancestry_file
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
sex_map = {}
with open("~{ped_file}") as f:
    for line in f:
        parts = line.strip().split('\t')
        sample_id = parts[1]
        sex_code = parts[4]
        if sex_code == "1":
            sex_map[sample_id] = "male"
        elif sex_code == "2":
            sex_map[sample_id] = "female"
        else:
            sex_map[sample_id] = "unknown"

pop_map = {}
with open("~{ancestry_file}") as f:
    for line in f:
        parts = line.strip().split('\t')
        pop_map[parts[0]] = parts[1]

with open("~{prefix}.metadata.tsv", "w") as f:
    f.write("SampleId\tSex\tPopulation\n")
    for sample_id in sex_map:
        sex = sex_map[sample_id]
        population = pop_map.get(sample_id, "unknown")
        f.write(f"{sample_id}\t{sex}\t{population}\n")
CODE
    >>>

    output {
        File metadata_file = "~{prefix}.metadata.tsv"
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
