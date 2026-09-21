version 1.0

import "../utils/Structs.wdl"

workflow CreateCohortAncestryFileAoUPhase2 {
    input {
        File ancestry_predictions
        Array[String] sample_ids
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_create_ancestry_file
    }

    call CreateAncestryFile {
        input:
            ancestry_predictions = ancestry_predictions,
            sample_ids = sample_ids,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_ancestry_file
    }

    output {
        File ancestry = CreateAncestryFile.ancestry_file
        File missing_samples = CreateAncestryFile.missing_samples_file
    }
}

task CreateAncestryFile {
    input {
        File ancestry_predictions
        Array[String] sample_ids
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import csv

with open("~{write_lines(sample_ids)}") as f:
    sample_ids = [line.strip() for line in f if line.strip()]

# Map each research ID to its predicted ancestry, ignoring the probability and PCA columns
pop_map = {}
with open("~{ancestry_predictions}") as f:
    for row in csv.DictReader(f, delimiter='\t'):
        pop_map[row["research_id"]] = row["ancestry_pred"]

absent = [sample_id for sample_id in sample_ids if sample_id not in pop_map]
if absent:
    raise ValueError(f"Samples absent from the ancestry predictions: {', '.join(absent)}")

with open("~{prefix}.ancestry.tsv", "w") as f:
    for sample_id in sample_ids:
        f.write(f"{sample_id}\t{pop_map[sample_id]}\n")

# Report predicted samples that the cohort does not include
with open("~{prefix}.missing_samples.txt", "w") as f:
    for sample_id in sorted(set(pop_map) - set(sample_ids)):
        f.write(f"{sample_id}\n")
CODE
    >>>

    output {
        File ancestry_file = "~{prefix}.ancestry.tsv"
        File missing_samples_file = "~{prefix}.missing_samples.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(ancestry_predictions, "GB")) + 10,
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
