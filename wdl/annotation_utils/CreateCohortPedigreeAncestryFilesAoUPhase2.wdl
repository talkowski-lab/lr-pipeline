version 1.0

import "../utils/Structs.wdl"

workflow CreateCohortPedigreeAncestryFilesAoUPhase2 {
    input {
        File ancestry_predictions
        Array[String] sample_ids
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_create_pedigree_ancestry_files
    }

    call CreatePedigreeAncestryFiles {
        input:
            ancestry_predictions = ancestry_predictions,
            sample_ids = sample_ids,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_pedigree_ancestry_files
    }

    output {
        File ped = CreatePedigreeAncestryFiles.ped_file
        File ancestry = CreatePedigreeAncestryFiles.ancestry_file
        File missing_samples = CreatePedigreeAncestryFiles.missing_samples_file
        File samples_missing_ancestry = CreatePedigreeAncestryFiles.samples_missing_ancestry_file
    }
}

task CreatePedigreeAncestryFiles {
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

# Rename the All of Us European label to the gnomAD label that compute_AFs.py accepts
POP_LABELS = {"eur": "nfe"}

with open("~{write_lines(sample_ids)}") as f:
    sample_ids = [line.strip() for line in f if line.strip()]

# Map each research ID to its predicted ancestry, ignoring the probability and PCA columns
pop_map = {}
with open("~{ancestry_predictions}") as f:
    for row in csv.DictReader(f, delimiter='\t'):
        ancestry_pred = row["ancestry_pred"]
        pop_map[row["research_id"]] = POP_LABELS.get(ancestry_pred, ancestry_pred)

# Write one singleton PED row per sample, with unknown parents, sex and phenotype
with open("~{prefix}.ped", "w") as f:
    for sample_id in sample_ids:
        f.write(f"{sample_id}\t{sample_id}\t0\t0\t0\t0\n")

# Label unpredicted samples '.' so compute_AFs.py counts them globally but in no population
with open("~{prefix}.ancestry.tsv", "w") as f:
    for sample_id in sample_ids:
        f.write(f"{sample_id}\t{pop_map.get(sample_id, '.')}\n")

# Report predicted samples that the cohort does not include
with open("~{prefix}.missing_samples.txt", "w") as f:
    for sample_id in sorted(set(pop_map) - set(sample_ids)):
        f.write(f"{sample_id}\n")

# Report cohort samples that the predictions do not cover
with open("~{prefix}.samples_missing_ancestry.txt", "w") as f:
    for sample_id in sample_ids:
        if sample_id not in pop_map:
            f.write(f"{sample_id}\n")
CODE
    >>>

    output {
        File ped_file = "~{prefix}.ped"
        File ancestry_file = "~{prefix}.ancestry.tsv"
        File missing_samples_file = "~{prefix}.missing_samples.txt"
        File samples_missing_ancestry_file = "~{prefix}.samples_missing_ancestry.txt"
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
