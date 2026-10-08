version 1.0

import "../utils/Structs.wdl"

workflow MergeVariantBaseTables {
    meta {
        description: [
            "This utility merges the per-contig site-level and sample-level tables written by `SummarizeVariantBases` into one site-level and one sample-level table, typically covering the whole genome.",
            "Every row's `contig_bases` and every `<name>_bases` column are summed across the input tables, and each `<name>_proportion` column is recomputed from those sums, so the proportions are of the combined length of the contigs merged. Every input table must have the same columns in the same order. Summing per-sample means across contigs gives the mean per genome only when every contig was summarized over the same samples."
        ]
    }

    parameter_meta {
        site_bases_tsvs: "Per-contig `site_bases_tsv` outputs of `SummarizeVariantBases`."
        sample_bases_tsvs: "Per-contig `sample_bases_tsv` outputs of `SummarizeVariantBases`."
        site_bases_tsv: "Site-level table with the bases summed across `site_bases_tsvs`."
        sample_bases_tsv: "Sample-level table with the bases summed across `sample_bases_tsvs`."
    }

    input {
        Array[File] site_bases_tsvs
        Array[File] sample_bases_tsvs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_sum_site_bases
        RuntimeAttr? runtime_attr_sum_sample_bases
    }

    call SumVariantBaseTables as SumSiteBases {
        input:
            bases_tsvs = site_bases_tsvs,
            prefix = "~{prefix}.site_bases",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_sum_site_bases
    }

    call SumVariantBaseTables as SumSampleBases {
        input:
            bases_tsvs = sample_bases_tsvs,
            prefix = "~{prefix}.sample_bases",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_sum_sample_bases
    }

    output {
        File site_bases_tsv = SumSiteBases.merged_tsv
        File sample_bases_tsv = SumSampleBases.merged_tsv
    }
}

task SumVariantBaseTables {
    input {
        Array[File] bases_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import sys
from collections import defaultdict

PATHS = "~{sep='\t' bases_tsvs}".split("\t")

header = None
categories = []
sums = defaultdict(float)

for path in PATHS:
    with open(path) as handle:
        file_header = handle.readline().rstrip("\n").split("\t")
        if header is None:
            header = file_header
        elif file_header != header:
            sys.exit(f"{path} has columns {file_header}, expected {header}")
        for line in handle:
            fields = dict(zip(header, line.rstrip("\n").split("\t")))
            category = fields["category"]
            if category not in categories:
                categories.append(category)
            for column in header:
                if column == "contig_bases" or column.endswith("_bases"):
                    sums[(category, column)] += float(fields[column])

# Recompute each proportion from the summed bases over the summed contig length
with open("~{prefix}.tsv", "w") as out:
    out.write("\t".join(header) + "\n")
    for category in categories:
        contig_bases = sums[(category, "contig_bases")]
        values = [category]
        for column in header[1:]:
            if column == "contig_bases":
                values.append(str(int(contig_bases)))
            elif column.endswith("_proportion"):
                bases = sums[(category, column[:-len("_proportion")] + "_bases")]
                values.append(f"{bases / contig_bases:.6f}")
            else:
                values.append(f"{sums[(category, column)]:.6f}")
        out.write("\t".join(values) + "\n")
CODE
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(bases_tsvs, "GB")) + 10,
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
