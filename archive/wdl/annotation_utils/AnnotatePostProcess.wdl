version 1.0

import "../utils/Helpers.wdl"

workflow AnnotatePostProcess {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        Array[String] include_contigs
        String prefix

        Int? records_per_shard

        Boolean filter_singletons = false
        File seqrepo_tar

        String utils_docker
        String vrs_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_post_process
        RuntimeAttr? runtime_attr_annotate_vrs
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.post_process",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call PostProcessVcf {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    include_contigs = include_contigs,
                    filter_singletons = filter_singletons,
                    prefix = "~{prefix}.~{contig}.post_process.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_post_process
            }

            call AnnotateVcfWithVRS {
                input:
                    vcf = PostProcessVcf.post_processed_vcf,
                    vcf_idx = PostProcessVcf.post_processed_vcf_idx,
                    seqrepo_tar = seqrepo_tar,
                    prefix = "~{prefix}.~{contig}.post_process.vrs.shard_~{i}",
                    docker = vrs_docker,
                    runtime_attr_override = runtime_attr_annotate_vrs
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatVcfs as ConcatShards {
                input:
                    vcfs = AnnotateVcfWithVRS.annotated_vcf,
                    vcf_idxs = AnnotateVcfWithVRS.annotated_vcf_idx,
                    allow_overlaps = false,
                    naive = true,
                    prefix = "~{prefix}.~{contig}.post_processed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_post_processed_vcf = select_first([ConcatShards.concat_vcf, AnnotateVcfWithVRS.annotated_vcf[0]])
        File final_post_processed_vcf_idx = select_first([ConcatShards.concat_vcf_idx, AnnotateVcfWithVRS.annotated_vcf_idx[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = final_post_processed_vcf,
                vcf_idxs = final_post_processed_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.post_processed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File post_processed_vcf = select_first([ConcatVcfs.concat_vcf, final_post_processed_vcf[0]])
        File post_processed_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, final_post_processed_vcf_idx[0]])
    }
}

task PostProcessVcf {
    input {
        File vcf
        File vcf_idx
        Array[String] include_contigs
        Boolean filter_singletons
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import sys
import pysam


filter_singletons = ~{true="True" false="False" filter_singletons}


def get_scalar(value):
    if isinstance(value, (list, tuple)):
        for item in value:
            if item is not None:
                return item
        return None
    return value


def prune_meis(record):
    allele_type = get_scalar(record.info.get("allele_type"))
    allele_length = get_scalar(record.info.get("allele_length"))
    if allele_type is None or allele_length is None:
        return

    length = abs(int(allele_length))
    if allele_type in {"alu_ins", "alu_del"} and (length < 250 or length > 350):
        record.info["allele_type"] = "ins" if "ins" in allele_type else "del"
        if "SUB_FAMILY" in record.info:
            del record.info["SUB_FAMILY"]
    elif allele_type in {"sva_ins", "sva_del"} and (length < 1000 or length > 4000):
        record.info["allele_type"] = "ins" if "ins" in allele_type else "del"
        if "SUB_FAMILY" in record.info:
            del record.info["SUB_FAMILY"]


def shortest_motif_length(record):
    motifs = record.info.get("MOTIFS")
    if motifs is None:
        return None

    motif_values = []
    if isinstance(motifs, (list, tuple)):
        raw_values = motifs
    else:
        raw_values = [motifs]

    for value in raw_values:
        if value is None:
            continue
        motif_values.extend(part for part in str(value).split(",") if part and part != ".")

    if not motif_values:
        return None
    return min(len(motif) for motif in motif_values)


def has_single_read_support(record):
    if not filter_singletons:
        return False

    if 'AC' in record.info and len(record.alts) == 1 and record.info['AC'][0] > 2:
        return False

    alt_depths = []
    for sample_data in record.samples.values():
        ad = sample_data.get("AD")
        if ad is None or len(ad) < 2:
            continue

        alt_depth = 0
        has_alt_depth = False
        for value in ad[1:]:
            if value is not None:
                alt_depth += value
                has_alt_depth = True

        if has_alt_depth and alt_depth > 0:
            alt_depths.append(alt_depth)

    return len(alt_depths) == 1 and alt_depths[0] == 1


def flush_buffer(buf, out_vcf):
    def custom_sort_key(rec):
        # 1. Sort by abs(INFO/allele_length)
        al_val = get_scalar(rec.info.get("allele_length"))
        abs_al = abs(int(al_val)) if al_val is not None else 0

        # 2. Natural sort by Variant ID (e.g., _2 before _10)
        id_val = rec.id if rec.id else ""
        parts = id_val.rsplit('_', 1)
        if len(parts) == 2 and parts[1].isdigit():
            id_sort = (parts[0], int(parts[1]))
        else:
            id_sort = (id_val, 0)

        return (abs_al, id_sort)

    buf.sort(key=custom_sort_key)
    for r in buf:
        out_vcf.write(r)


# Update header
vcf_in = pysam.VariantFile("~{vcf}")
header = vcf_in.header.copy()
if "HOMOPOLYMER_TRV" not in header.info:
    header.info.add("HOMOPOLYMER_TRV", 0, "Flag", "Tandem repeat call where the shortest motif has length 1.")
if filter_singletons and "SINGLE_READ_SUPPORT" not in header.filters:
    header.filters.add("SINGLE_READ_SUPPORT", None, None, "Variant supported by a single read in a single sample.")

# Set up variables
vcf_out = pysam.VariantFile("~{prefix}.processed.vcf.gz", "wz", header=header)

buffer = []
current_chrom = None
current_pos = None

for record in vcf_in:
    record.translate(vcf_out.header)

    # Revise MEIs outside expected size range to indel
    prune_meis(record)

    # Flag homopolymer TR variants
    if get_scalar(record.info.get("allele_type")) == "trv" and shortest_motif_length(record) == 1:
        record.info["HOMOPOLYMER_TRV"] = True

    # Flag variants with single read support
    if has_single_read_support(record):
        record.filter.add("SINGLE_READ_SUPPORT")

    # Buffer records to sort them if they fall on the exact same coordinate
    if record.chrom != current_chrom or record.pos != current_pos:
        if buffer:
            flush_buffer(buffer, vcf_out)
        buffer = [record]
        current_chrom = record.chrom
        current_pos = record.pos
    else:
        buffer.append(record)

if buffer:
    flush_buffer(buffer, vcf_out)

vcf_in.close()
vcf_out.close()
CODE

        bcftools view -h ~{prefix}.processed.vcf.gz \
            | awk '/^##fileformat=|^##contig=|^##FILTER=|^##INFO=|^##FORMAT=|^#CHROM/ { print }' \
            > clean_header.txt

        python3 <<CODE
import json
import re

with open("~{write_json(include_contigs)}", "r") as handle:
    include_contigs = json.load(handle)

fileformat_line = None
chrom_line = None
contig_lines = {}
filter_lines = []
info_lines = []
format_lines = []

with open("clean_header.txt", "r") as handle:
    for line in handle:
        if line.startswith("##fileformat="):
            fileformat_line = line
        elif line.startswith("##contig="):
            match = re.search(r"ID=([^,>]+)", line)
            if match:
                contig_lines[match.group(1)] = line
        elif line.startswith("##FILTER="):
            filter_lines.append(line)
        elif line.startswith("##INFO="):
            info_lines.append(line)
        elif line.startswith("##FORMAT="):
            format_lines.append(line)
        elif line.startswith("#CHROM"):
            chrom_line = line

with open("clean_header.grouped.txt", "w") as handle:
    if fileformat_line is not None:
        handle.write(fileformat_line)
    for contig_name in include_contigs:
        if contig_name in contig_lines:
            handle.write(contig_lines[contig_name])
    for line in filter_lines:
        handle.write(line)
    for line in info_lines:
        handle.write(line)
    for line in format_lines:
        handle.write(line)
    if chrom_line is not None:
        handle.write(chrom_line)
CODE

        bcftools reheader \
            -h clean_header.grouped.txt \
            -o ~{prefix}.vcf.gz \
            ~{prefix}.processed.vcf.gz

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File post_processed_vcf = "~{prefix}.vcf.gz"
        File post_processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(size(vcf, "GB")) + 5,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 20,
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

task AnnotateVcfWithVRS {
    input {
        File vcf
        File vcf_idx
        File seqrepo_tar
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir seqrepo_data
        tar -xzf ~{seqrepo_tar} -C seqrepo_data --strip-components 1

        SEQREPO_PATH=$(find $(pwd)/seqrepo_data -name "aliases.sqlite3" | xargs dirname)

        vrs-annotate vcf \
            --dataproxy-uri="seqrepo+file://${SEQREPO_PATH}" \
            --vcf-out ~{prefix}.vcf.gz \
            --vrs-attributes \
            ~{vcf}

        rm -rf seqrepo_data

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(seqrepo_tar, "GB")) + 20,
        boot_disk_gb: 50,
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
