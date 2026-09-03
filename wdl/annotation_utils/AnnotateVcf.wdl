version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateVcf {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[String] contigs
        String prefix

        Int? records_per_shard
        Array[String]? strip_info_fields

        Array[Boolean] sort_tsvs = []
        Array[String] subset_vcf_strings = []
        Array[String] awk_tsv_conditions = []
        Array[Array[Int]] subset_tsv_columns = []

        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_strip
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1
    Boolean sort_requested = length(sort_tsvs) > 0

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }
        }

        if (!single_contig || sort_requested) {
            scatter (i in range(length(annotations_tsvs))) {
                call Helpers.SubsetTsvToContig {
                    input:
                        tsv = annotations_tsvs[i],
                        contig = contig,
                        sort_output = if sort_requested then sort_tsvs[i] else false,
                        prefix = "~{prefix}.~{contig}.tsv~{i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_tsv
                }
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])
        Array[File] contig_tsvs = select_first([SubsetTsvToContig.subset_tsv, annotations_tsvs])

        if (defined(strip_info_fields)) {
            call Helpers.StripInfoFields {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    info_fields = select_first([strip_info_fields]),
                    prefix = "~{prefix}.~{contig}.stripped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_strip
            }
        }

        File stripped_vcf = select_first([StripInfoFields.stripped_vcf, contig_vcf])
        File stripped_vcf_idx = select_first([StripInfoFields.stripped_vcf_idx, contig_vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = stripped_vcf,
                    vcf_idx = stripped_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [stripped_vcf]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [stripped_vcf_idx]])

        scatter (s in range(length(shard_vcfs))) {
            call AnnotateSequentially {
                input:
                    vcf = shard_vcfs[s],
                    vcf_idx = shard_vcf_idxs[s],
                    annotations_tsvs = contig_tsvs,
                    subset_vcf_strings = subset_vcf_strings,
                    awk_tsv_conditions = awk_tsv_conditions,
                    subset_tsv_columns = subset_tsv_columns,
                    info_names = info_names,
                    info_descriptions = info_descriptions,
                    info_types = info_types,
                    info_numbers = info_numbers,
                    prefix = "~{prefix}.~{contig}.annotated.shard_~{s}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatShards {
                input:
                    vcfs = AnnotateSequentially.annotated_vcf,
                    vcf_idxs = AnnotateSequentially.annotated_vcf_idx,
                    allow_overlaps = false,
                    naive = true,
                    prefix = "~{prefix}.~{contig}.annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File contig_annotated_vcf = select_first([ConcatShards.concat_vcf, AnnotateSequentially.annotated_vcf[0]])
        File contig_annotated_vcf_idx = select_first([ConcatShards.concat_vcf_idx, AnnotateSequentially.annotated_vcf_idx[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = contig_annotated_vcf,
                vcf_idxs = contig_annotated_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotated_vcf = select_first([ConcatVcfs.concat_vcf, contig_annotated_vcf[0]])
        File annotated_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, contig_annotated_vcf_idx[0]])
    }
}

task AnnotateSequentially {
    input {
        File vcf
        File vcf_idx
        Array[File] annotations_tsvs
        Array[String] subset_vcf_strings = []
        Array[String] awk_tsv_conditions = []
        Array[Array[Int]] subset_tsv_columns = []
        Array[Array[String]] info_names
        Array[Array[String]] info_descriptions
        Array[Array[String]] info_types
        Array[Array[String]] info_numbers
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<EOF
import sys

info_names = [line.strip().split('\t') for line in open('~{write_tsv(info_names)}')]
info_descriptions = [line.strip().split('\t') for line in open('~{write_tsv(info_descriptions)}')]
info_types = [line.strip().split('\t') for line in open('~{write_tsv(info_types)}')]
info_numbers = [line.strip().split('\t') for line in open('~{write_tsv(info_numbers)}')]
subset_tsv_columns = [line.strip().split('\t') for line in open('~{write_tsv(subset_tsv_columns)}')]
subset_tsv_columns = [cols if cols != [''] else [] for cols in subset_tsv_columns]

if len(info_names) != len(info_descriptions) or len(info_names) != len(info_types) or len(info_names) != len(info_numbers):
    sys.stderr.write("Error: All info arrays must have the same length.\n")
    sys.exit(1)

for i, (names, descs, types, numbers) in enumerate(zip(info_names, info_descriptions, info_types, info_numbers)):
    if len(names) != len(descs) or len(names) != len(types) or len(names) != len(numbers):
        sys.stderr.write(f"Error: info arrays at index {i} must have the same length.\n")
        sys.exit(1)
    
    with open(f"header_{i}.txt", "w") as f:
        for name, desc, type_val, number in zip(names, descs, types, numbers):
            f.write(f'##INFO=<ID={name},Number={number},Type={type_val},Description="{desc}">\n')
    
    column_spec = ','.join(['CHROM', 'POS', 'REF', 'ALT', '~ID'] + [f'INFO/{name}' for name in names])
    with open(f"columns_{i}.txt", "w") as f:
        f.write(column_spec)

    cols = subset_tsv_columns[i] if i < len(subset_tsv_columns) else []
    with open(f"colsel_{i}.txt", "w") as f:
        if cols:
            fields = ','.join(['$1', '$2', '$3', '$4', '$5'] + [f'${c}' for c in cols])
            f.write('BEGIN{OFS="\\t"}{print ' + fields + '}')

with open("num_tsvs.txt", "w") as f:
    f.write(str(len(info_names)))
EOF

        current_vcf="~{vcf}"
        SUBSET_FILE="~{write_lines(subset_vcf_strings)}"
        AWK_COND_FILE="~{write_lines(awk_tsv_conditions)}"     
        i=0
        for tsv_file in ~{sep=' ' annotations_tsvs}; do
            AWK_ARG=""
            if [ -s "$AWK_COND_FILE" ]; then
                AWK_ARG=$(sed -n "$((i + 1))p" "$AWK_COND_FILE")
            fi

            COLSEL_ARG=$(cat "colsel_${i}.txt")

            if [ -n "$AWK_ARG" ] && [ -n "$COLSEL_ARG" ]; then
                awk -F'\t' "$AWK_ARG" "$tsv_file" | awk -F'\t' "$COLSEL_ARG" | bgzip -c > "annotations_${i}.tsv.gz"
            elif [ -n "$AWK_ARG" ]; then
                awk -F'\t' "$AWK_ARG" "$tsv_file" | bgzip -c > "annotations_${i}.tsv.gz"
            elif [ -n "$COLSEL_ARG" ]; then
                awk -F'\t' "$COLSEL_ARG" "$tsv_file" | bgzip -c > "annotations_${i}.tsv.gz"
            else
                bgzip -c "$tsv_file" > "annotations_${i}.tsv.gz"
            fi
            tabix -s1 -b2 -e2 "annotations_${i}.tsv.gz"
            
            COLUMN_SPEC=$(cat "columns_${i}.txt")

            SUBSET_ARG=""
            if [ -s "$SUBSET_FILE" ]; then
                SUBSET_ARG=$(sed -n "$((i + 1))p" "$SUBSET_FILE")
            fi
            
            bcftools annotate \
                -a "annotations_${i}.tsv.gz" \
                -h "header_${i}.txt" \
                -c "$COLUMN_SPEC" \
                $SUBSET_ARG \
                -Oz -o "temp_${i}.vcf.gz" \
                "$current_vcf"
            current_vcf="temp_${i}.vcf.gz"
            
            i=$((i + 1))
        done
        
        mv "$current_vcf" ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * length(annotations_tsvs) * ceil(size(vcf, "GB") + size(annotations_tsvs, "GB")) + 10,
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
