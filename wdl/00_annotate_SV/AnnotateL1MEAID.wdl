version 1.0

import "../utils/Helpers.wdl"
import "../tools/RepeatMasker.wdl"
import "../utils/Structs.wdl"

workflow AnnotateL1MEAID {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Int min_length

        String utils_docker
        String repeatmasker_docker
        String l1meaid_docker
        String intact_mei_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_ins_to_fa
        RuntimeAttr? runtime_attr_repeat_masker
        RuntimeAttr? runtime_attr_limeaid
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_concat_shards_intactmei
        RuntimeAttr? runtime_attr_concat_shards_annotations
        RuntimeAttr? runtime_attr_concat_contigs
        RuntimeAttr? runtime_attr_concat_intactmei
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        call Helpers.SubsetVcfByArgs {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                include_args = "abs(INFO/allele_length) >= ~{min_length} && INFO/allele_type = \"ins\"",
                extra_args = if single_contig then "" else "--regions " + contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetVcfByArgs.subset_vcf,
                    vcf_idx = SubsetVcfByArgs.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [SubsetVcfByArgs.subset_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [SubsetVcfByArgs.subset_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call RepeatMasker.RepeatMasker {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    prefix = "~{prefix}.~{contig}.shard_~{i}.rm",
                    utils_docker = utils_docker,
                    repeatmasker_docker = repeatmasker_docker,
                    runtime_attr_ins_to_fa = runtime_attr_ins_to_fa,
                    runtime_attr_repeat_masker = runtime_attr_repeat_masker
            }

            call L1MEAID {
                input:
                    rm_fa = RepeatMasker.rm_fa,
                    rm_out = RepeatMasker.rm_out,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.l1meaid",
                    docker = l1meaid_docker,
                    runtime_attr_override = runtime_attr_limeaid
            }

            call IntactMEI {
                input:
                    l1meaid_output = L1MEAID.l1meaid_output,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.intactmei",
                    docker = intact_mei_docker,
                    runtime_attr_override = runtime_attr_filter
            }

            call GenerateAnnotationTable {
                input:
                    vcf = vcfs_to_process[i],
                    filtered_tsv = IntactMEI.filtered_output,
                    prefix = "~{prefix}.~{contig}.shard_~{i}.intactmei_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatAnnotationShards {
                input:
                    tsvs = GenerateAnnotationTable.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.intactmei_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards_annotations
            }
        }

        File final_annotations_tsv = select_first([ConcatAnnotationShards.concatenated_tsv, GenerateAnnotationTable.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotations {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.intactmei_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_contigs
        }
    }

    output {
        File annotations_tsv_l1meaid = select_first([MergeAnnotations.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task L1MEAID {
    input {
        File rm_fa
        File rm_out
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/src/L1ME-AID/limeaid.py \
            -i ~{rm_fa} \
            -r ~{rm_out} \
            -o ~{prefix}.txt
    >>>

    output {
        File l1meaid_output = "~{prefix}.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(rm_fa, "GB")) + 5,
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

task IntactMEI {
    input {
        File l1meaid_output
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        perl /opt/src/utility/limeaid.filter.pl \
            ~{l1meaid_output} \
            > ~{prefix}.tsv
    >>>

    output {
        File filtered_output = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(l1meaid_output, "GB")) + 5,
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

task GenerateAnnotationTable {
    input {
        File vcf
        File filtered_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ~{vcf} > vcf_lookup.tsv

        python3 <<EOF
import sys

vcf_lookup_file = "vcf_lookup.tsv"
input_tsv = "~{filtered_tsv}"
output_anno = "~{prefix}.tsv"

vcf_lookup = {}
with open(vcf_lookup_file, 'r') as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) >= 5:
            chrom, pos, ref, alt, var_id = fields[0], fields[1], fields[2], fields[3], fields[4]
            vcf_lookup[(chrom, pos, ref, alt)] = var_id

with open(input_tsv, 'r') as f_in, open(output_anno, 'w') as f_out:
    for line in f_in:
        parts = line.strip().split('\t')
        if len(parts) < 12:
            continue

        classification = parts[8]
        structure = parts[10]
        me_type = None
        if classification == "SINE/Alu" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "ALU"
        elif classification == "Retroposon/SVA" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "SVA"
        elif classification == "LINE/L1" and (structure == "INTACT" or structure == "INTACT_3end"):
            me_type = "LINE"

        if me_type:
            full_id = parts[0]
            full_id_parts = full_id.split(';')
            chrom = full_id_parts[0].split(':')[0]
            pos = full_id_parts[0].split(':')[1]
            ref = full_id_parts[1].split('_')[0]
            alt = parts[1]
            subfam = parts[4]
            key = (chrom, pos, ref, alt)
            if key in vcf_lookup:
                f_out.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{vcf_lookup[key]}\t{me_type}\t{subfam}\n")
EOF
    >>>

    output {
        File annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(filtered_tsv, "GB")) + 5,
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
