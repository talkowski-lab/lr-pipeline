version 1.0

import "../utils/Structs.wdl"

workflow MergeSummarizedCallsetCounts {
    meta {
        description: [
            "This utility sums the site and per-sample count tables that `SummarizeMergedCallsets` writes for separate contigs into one site table and one per-sample table. Rows are matched on their label columns and every count column is summed, so a row that appears for only some contigs keeps the counts of those contigs."
        ]
    }

    parameter_meta {
        site_counts_tsvs: "Site count tables from `SummarizeMergedCallsets`, such as one per contig."
        sample_counts_tsvs: "Per-sample count tables from `SummarizeMergedCallsets`, from the same runs as `site_counts_tsvs`."
        site_counts_tsv: "Number of matched sites and of sites unique to each callset, summed across the inputs."
        sample_counts_tsv: "Per-sample counts of matched and unique sites, summed across the inputs."
    }

    input {
        Array[File] site_counts_tsvs
        Array[File] sample_counts_tsvs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_sum
    }

    call SumCallsetCountTables {
        input:
            site_counts_tsvs = site_counts_tsvs,
            sample_counts_tsvs = sample_counts_tsvs,
            prefix = prefix,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_sum
    }

    output {
        File site_counts_tsv = SumCallsetCountTables.site_counts_tsv
        File sample_counts_tsv = SumCallsetCountTables.sample_counts_tsv
    }
}

task SumCallsetCountTables {
    input {
        Array[File] site_counts_tsvs
        Array[File] sample_counts_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import pandas as pd

SITE_TSVS = "~{sep=',' site_counts_tsvs}".split(",")
SAMPLE_TSVS = "~{sep=',' sample_counts_tsvs}".split(",")
SITE_OUTPUT = "~{prefix}.site_counts.tsv"
SAMPLE_OUTPUT = "~{prefix}.sample_counts.tsv"
SITE_KEYS = ["concordance", "region", "status"]
SAMPLE_KEYS = ["sample", "callset", "concordance", "region", "status"]
SAMPLE_ORDER = ["callset", "sample", "concordance", "region", "status"]
FIXED_ORDERS = {
    "concordance": ["All", "dbSNP Missing", "gnomAD Missing", "dbSNP/gnomAD Missing"],
    "region": ["All", "US", "RM", "SD", "SR"],
}


# Sum the tables on their label columns, ordering concordance and region as SummarizeMergedCallsets does
def sum_tables(paths, keys, order, output):
    table = pd.concat(
        [pd.read_csv(p, sep="\t", dtype={k: str for k in keys}, keep_default_na=False) for p in paths],
        ignore_index=True,
    )
    for key in keys:
        table[key] = pd.Categorical(table[key], FIXED_ORDERS.get(key, list(dict.fromkeys(table[key]))))
    summed = table.groupby(keys, observed=True).sum().reset_index().sort_values(order)
    summed.to_csv(output, sep="\t", index=False)


sum_tables(SITE_TSVS, SITE_KEYS, SITE_KEYS, SITE_OUTPUT)
sum_tables(SAMPLE_TSVS, SAMPLE_KEYS, SAMPLE_ORDER, SAMPLE_OUTPUT)
CODE
    >>>

    output {
        File site_counts_tsv = "~{prefix}.site_counts.tsv"
        File sample_counts_tsv = "~{prefix}.sample_counts.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 2 * ceil(size(site_counts_tsvs, "GB") + size(sample_counts_tsvs, "GB")) + 10,
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
