version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow VcfToBed {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Int? min_length
        String length_field = "allele_length"
        Boolean convert_to_ins_del = false
        Array[Array[String]]? switch_info_fields
        Array[String] info_columns = ["ALL"]
        Boolean include_samples = true
        Boolean include_filters = true
        Boolean split_bnd = false
        Boolean split_cpx = false
        Boolean output_gz = false
        Boolean output_bed = false

        String gatk_sv_lr_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_length
        RuntimeAttr? runtime_attr_modify_info
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_rename
    }

    scatter (i in range(length(contigs))) {
        if (defined(min_length)) {
            call Helpers.SubsetVcfByLength {
                input:
                    vcf = vcfs[i],
                    vcf_idx = vcf_idxs[i],
                    length_field = length_field,
                    min_length = min_length,
                    prefix = "~{prefix}.~{contigs[i]}.length_filtered",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_length
            }
        }

        File length_vcf = select_first([SubsetVcfByLength.subset_vcf, vcfs[i]])
        File length_vcf_idx = select_first([SubsetVcfByLength.subset_vcf_idx, vcf_idxs[i]])

        if (convert_to_ins_del || defined(switch_info_fields)) {
            call ModifyInfoFields {
                input:
                    vcf = length_vcf,
                    vcf_idx = length_vcf_idx,
                    convert_to_ins_del = convert_to_ins_del,
                    switch_info_fields = select_first([switch_info_fields, []]),
                    prefix = "~{prefix}.~{contigs[i]}.modified",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_modify_info
            }
        }

        File processed_vcf = select_first([ModifyInfoFields.modified_vcf, length_vcf])
        File processed_vcf_idx = select_first([ModifyInfoFields.modified_vcf_idx, length_vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = processed_vcf,
                    vcf_idx = processed_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contigs[i]}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [processed_vcf]])

        scatter (j in range(length(shard_vcfs))) {
            call SvtkVcfToBed {
                input:
                    vcf = shard_vcfs[j],
                    info_columns = info_columns,
                    include_samples = include_samples,
                    include_filters = include_filters,
                    split_bnd = split_bnd,
                    split_cpx = split_cpx,
                    prefix = "~{prefix}.~{contigs[i]}.shard_~{j}",
                    docker = gatk_sv_lr_docker,
                    runtime_attr_override = runtime_attr_vcf2bed
            }
        }
    }

    call Helpers.ConcatTsvs as ConcatBedShards {
        input:
            tsvs = flatten(SvtkVcfToBed.bed),
            sort_output = false,
            preserve_header = true,
            compressed_output = false,
            prefix = "~{prefix}.vcf2bed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    call RenameOutput {
        input:
            bed = ConcatBedShards.concatenated_tsv,
            output_gz = output_gz,
            output_bed = output_bed,
            prefix = "~{prefix}.vcf2bed",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_rename
    }

    output {
        File bed = RenameOutput.renamed_bed
    }
}

task SvtkVcfToBed {
    input {
        File vcf
        Array[String] info_columns
        Boolean include_samples
        Boolean include_filters
        Boolean split_bnd
        Boolean split_cpx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        svtk vcf2bed \
            ~{sep=" " prefix("--info ", info_columns)} \
            ~{if include_samples then "" else "--no-samples"} \
            ~{if include_filters then "--include-filters" else ""} \
            ~{if split_bnd then "--split-bnd" else ""} \
            ~{if split_cpx then "--split-cpx" else ""} \
            ~{vcf} \
            ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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

task RenameOutput {
    input {
        File bed
        Boolean output_gz
        Boolean output_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    String base_ext = if output_bed then "bed" else "tsv"

    command <<<
        set -euo pipefail

        mv ~{bed} ~{prefix}.~{base_ext}
        
        if [ "~{output_gz}" == "true" ]; then
            gzip -1 ~{prefix}.~{base_ext}
        fi
    >>>

    output {
        File renamed_bed = if output_gz then "~{prefix}.~{base_ext}.gz" else "~{prefix}.~{base_ext}"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: 2 * ceil(size(bed, "GB")) + 5,
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

task ModifyInfoFields {
    input {
        File vcf
        File vcf_idx
        Boolean convert_to_ins_del
        Array[Array[String]] switch_info_fields
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    File switch_tsv = write_tsv(switch_info_fields)

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

convert = ~{true="True" false="False" convert_to_ins_del}

pairs = []
with open("~{switch_tsv}") as f:
    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue
        src, dst = line.split("\t")
        pairs.append((src, dst))

vcf_in = pysam.VariantFile("~{vcf}")
header = vcf_in.header

if convert and 'SVTYPE' not in header.info:
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type inferred from variant ID (INS/DEL)">')
for src, dst in pairs:
    if dst not in header.info:
        src_hdr = header.info[src]
        header.add_line(f'##INFO=<ID={dst},Number={src_hdr.number},Type={src_hdr.type},Description="Copied from {src}">')

vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=header)
for record in vcf_in:
    if convert:
        vid = record.id or ""
        if 'INS' in vid:
            record.info['SVTYPE'] = 'INS'
        elif 'DEL' in vid:
            record.info['SVTYPE'] = 'DEL'
    for src, dst in pairs:
        if src in record.info:
            record.info[dst] = record.info[src]
    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File modified_vcf = "~{prefix}.vcf.gz"
        File modified_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
