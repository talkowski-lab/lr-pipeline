version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow SubsetTsvToColumns {
    input {
        File annotations_tsv
        File annotations_header
        Array[String] subset_columns
        String prefix

        Array[Array[String]]? subset_column_values

        String utils_docker

        RuntimeAttr? runtime_attr_subset
    }

    call SubsetAndFilterTsv {
        input:
            tsv_file = annotations_tsv,
            header_file = annotations_header,
            columns_to_keep = subset_columns,
            column_filter_values = subset_column_values,
            prefix = "~{prefix}.subset",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset
    }

    output {
        File subset_tsv = SubsetAndFilterTsv.out_tsv
    }
}

task SubsetAndFilterTsv {
    input {
        File tsv_file
        File header_file
        Array[String] columns_to_keep
        Array[Array[String]]? column_filter_values
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import json
import csv
import sys

tsv_path = "~{tsv_file}"
header_path = "~{header_file}"
out_path = "~{prefix}.tsv"

with open("~{write_json(columns_to_keep)}") as f:
    cols_requested = json.load(f)

filter_json_path = "~{write_json(select_first([column_filter_values, []]))}"
with open(filter_json_path) as f:
    filter_vals_list = json.load(f)

col_name_to_idx = {}
with open(header_path, 'r') as f:
    for i, line in enumerate(f):
        name = line.strip()
        if name:
            col_name_to_idx[name] = i + 5

output_idx = [0, 1, 2, 3, 4]
active_filters = {}
for i, col_name in enumerate(cols_requested):
    if col_name not in col_name_to_idx:
        sys.stderr.write(f"Error: Column '{col_name}' not found in header file.\n")
        sys.exit(1)
    tsv_idx = col_name_to_idx[col_name]
    output_idx.append(tsv_idx)

    if i < len(filter_vals_list):
        allowed_values = filter_vals_list[i]
        if allowed_values:
            active_filters[tsv_idx] = set(allowed_values)

with open(tsv_path, 'r') as fin, open(out_path, 'w', newline='') as fout:
    reader = csv.reader(fin, delimiter='\t')
    writer = csv.writer(fout, delimiter='\t')

    for row in reader:
        if not row:
            continue

        keep_row = True
        for filter_col_idx, allowed_set in active_filters.items():
            if filter_col_idx >= len(row):
                keep_row = False
                break
            if row[filter_col_idx] not in allowed_set:
                keep_row = False
                break

        if keep_row:
            out_row = []
            for idx in output_idx:
                if idx < len(row):
                    out_row.append(row[idx])
                else:
                    out_row.append(".")

            writer.writerow(out_row)
CODE
    >>>

    output {
        File out_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv_file, "GB")) + 5,
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
