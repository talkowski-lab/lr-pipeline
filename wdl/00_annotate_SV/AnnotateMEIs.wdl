version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateMEIs {
    input {
        File annotations_tsv_l1meaid
        File annotations_tsv_palmer
        File annotations_tsv_svan
        File annotations_header_svan
        Array[String] contigs
        String prefix

        String utils_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetTsvToContig as SubsetL1MEAID {
                input:
                    tsv = annotations_tsv_l1meaid,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.l1meaid",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }

            call Helpers.SubsetTsvToContig as SubsetPALMER {
                input:
                    tsv = annotations_tsv_palmer,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.palmer",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }

            call Helpers.SubsetTsvToContig as SubsetSVAN {
                input:
                    tsv = annotations_tsv_svan,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.svan",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }
        }

        File contig_l1meaid = select_first([SubsetL1MEAID.subset_tsv, annotations_tsv_l1meaid])
        File contig_palmer = select_first([SubsetPALMER.subset_tsv, annotations_tsv_palmer])
        File contig_svan = select_first([SubsetSVAN.subset_tsv, annotations_tsv_svan])

        call FilterMEIs {
            input:
                tsv_l1meaid = contig_l1meaid,
                tsv_palmer = contig_palmer,
                tsv_svan = contig_svan,
                header_svan = annotations_header_svan,
                prefix = "~{prefix}.~{contig}.filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter
        }
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = FilterMEIs.filtered_tsv,
                sort_output = false,
                prefix = "~{prefix}.mei_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File annotations_tsv_meis = select_first([ConcatTsvs.concatenated_tsv, FilterMEIs.filtered_tsv[0]])
    }
}

task FilterMEIs {
    input {
        File tsv_l1meaid
        File tsv_palmer
        File tsv_svan
        File header_svan
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
        import csv

        l1_path = "~{tsv_l1meaid}"
        palmer_path = "~{tsv_palmer}"
        svan_path = "~{tsv_svan}"
        header_path = "~{header_svan}"
        out_path = "~{prefix}.filtered.tsv"

        fam_n_idx = -1
        with open(header_path, 'r') as f:
            for i, line in enumerate(f):
                if line.strip() == "FAM_N":
                    fam_n_idx = i + 5
                    break

        palmer_set = set()
        with open(palmer_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if row:
                    key = tuple(row[:6])
                    palmer_set.add(key)

        svan_set = set()
        with open(svan_path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if row:
                    key = tuple(row[:5]) + (row[fam_n_idx].upper(),)
                    svan_set.add(key)

        with open(l1_path, 'r') as fin, open(out_path, 'w') as fout:
            reader = csv.reader(fin, delimiter='\t')
            writer = csv.writer(fout, delimiter='\t')
            for row in reader:
                if not row:
                    continue
                me_type = row[5]
                if me_type == "LINE":
                    writer.writerow(row)
                else:
                    key_6 = tuple(row[:6])
                    if key_6 in palmer_set or key_6 in svan_set:
                        writer.writerow(row)
        CODE
    >>>

    output {
        File filtered_tsv = "~{prefix}.filtered.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size([tsv_l1meaid, tsv_palmer, tsv_svan], "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
