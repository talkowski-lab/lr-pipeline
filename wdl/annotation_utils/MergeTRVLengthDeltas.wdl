version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeTRVLengthDeltas {
    meta {
        description: [
            "This utility merges the per-contig outputs of `ExtractTRVLengthDeltas` into callset-wide tables. It concatenates the per-allele tables and sums the per-length-delta count tables within each length delta."
        ]
    }

    parameter_meta {
        allele_deltas_tsvs: "Per-contig `ExtractTRVLengthDeltas.allele_deltas_tsv` outputs."
        length_deltas_tsvs: "Per-contig `ExtractTRVLengthDeltas.length_deltas_tsv` outputs."
        merged_allele_deltas_tsv: "Gzipped per-allele table concatenated across all inputs, with the same columns as the inputs."
        merged_length_deltas_tsv: "Gzipped count table with columns `length_delta`, `n_alleles` and `allele_count`, summed across all inputs."
    }

    input {
        Array[File] allele_deltas_tsvs
        Array[File] length_deltas_tsvs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_sum
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = allele_deltas_tsvs,
            sort_output = false,
            preserve_header = true,
            compressed_tsvs = true,
            compressed_output = true,
            prefix = "~{prefix}.trv_allele_deltas",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call SumLengthDeltas {
        input:
            tsvs = length_deltas_tsvs,
            prefix = "~{prefix}.trv_length_deltas",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_sum
    }

    output {
        File merged_allele_deltas_tsv = ConcatTsvs.concatenated_tsv
        File merged_length_deltas_tsv = SumLengthDeltas.summed_tsv
    }
}

task SumLengthDeltas {
    input {
        Array[File] tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv
import gzip
from collections import defaultdict

INPUTS = ["~{sep='", "' tsvs}"]
OUTPUT = "~{prefix}.tsv.gz"

n_alleles = defaultdict(int)
allele_count = defaultdict(int)
for path in INPUTS:
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            length_delta = int(row["length_delta"])
            n_alleles[length_delta] += int(row["n_alleles"])
            allele_count[length_delta] += int(row["allele_count"])

with gzip.open(OUTPUT, "wt", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["length_delta", "n_alleles", "allele_count"])
    for length_delta in sorted(n_alleles):
        writer.writerow([length_delta, n_alleles[length_delta], allele_count[length_delta]])
PYCODE
    >>>

    output {
        File summed_tsv = "~{prefix}.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(tsvs, "GB")) + 10,
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
