version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateGQMetrics {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Array[String] gq_fields
        Array[Array[Int]] gq_bins
        Array[String] gq_variant_filters
        Array[Boolean] gq_larger_field
        Boolean ab_annotation
        Array[Float] ab_bins

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_generate_tsv
        RuntimeAttr? runtime_attr_generate_ab_tsv
        RuntimeAttr? runtime_attr_merge_aligned_tsv
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_vcf
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
                    runtime_attr_override = runtime_attr_subset_vcf
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
                    prefix = "~{prefix}.~{contig}.gq_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (shard_i in range(length(vcfs_to_process))) {
            scatter (field_i in range(length(gq_fields))) {
                call GenerateGQAnnotationTsv {
                    input:
                        vcf = vcfs_to_process[shard_i],
                        vcf_idx = vcf_idxs_to_process[shard_i],
                        gq_field = gq_fields[field_i],
                        gq_bins = gq_bins[field_i],
                        gq_variant_filter = gq_variant_filters[field_i],
                        gq_larger_field = gq_larger_field[field_i],
                        prefix = "~{prefix}.~{contig}.~{gq_fields[field_i]}.shard_~{shard_i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_generate_tsv
                }
            }

            if (ab_annotation) {
                call GenerateABAnnotationTsv {
                    input:
                        vcf = vcfs_to_process[shard_i],
                        vcf_idx = vcf_idxs_to_process[shard_i],
                        ab_bins = ab_bins,
                        prefix = "~{prefix}.~{contig}.ab_hist_alt.shard_~{shard_i}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_generate_ab_tsv
                }
            }

            call Helpers.MergeAlignedTsvs as MergeShardAnnotations {
                input:
                    tsvs = flatten([GenerateGQAnnotationTsv.annotation_tsv, select_all([GenerateABAnnotationTsv.annotation_tsv])]),
                    prefix = "~{prefix}.~{contig}.gq_annotations.shard_~{shard_i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_merge_aligned_tsv
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = MergeShardAnnotations.merged_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.gq_annotations",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, MergeShardAnnotations.merged_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.gq_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    output {
        File annotations_tsv_gq = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task GenerateGQAnnotationTsv {
    input {
        File vcf
        File vcf_idx
        String gq_field
        Array[Int] gq_bins
        String gq_variant_filter
        Boolean gq_larger_field
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        FILTER='~{gq_variant_filter}'
        INPUT_VCF="~{vcf}"
        if [ -n "$FILTER" ] && [ "$FILTER" != "." ] && [ "$FILTER" != "None" ]; then
            bcftools view -i "$FILTER" -Oz -o filtered.vcf.gz "~{vcf}"
            tabix -p vcf filtered.vcf.gz
            INPUT_VCF="filtered.vcf.gz"
        fi

        python3 <<CODE
import json
import pysam

field = "~{gq_field}"
with open("~{write_json(gq_bins)}") as f:
    bins = json.load(f)
larger_flag = ~{true="True" false="False" gq_larger_field}
prefix = "~{prefix}"

bin_edges_str = "|".join(map(str, bins))
col_names = [f"{field}_hist_all_bin_freq", f"{field}_hist_alt_bin_freq"]
if larger_flag:
    col_names += [f"{field}_hist_all_n_larger", f"{field}_hist_alt_n_larger"]
n_intervals = len(bins) - 1

def get_bin_index(val):
    for j in range(n_intervals):
        if val < bins[j + 1]:
            return j
    return None

vcf_in = pysam.VariantFile("$INPUT_VCF")

with open(f"{prefix}.tsv", "w") as out:
    out.write("\t".join(["#CHROM", "POS", "REF", "ALT", "ID"] + col_names) + "\n")
    for record in vcf_in:
        counts_all = [0] * n_intervals
        counts_alt = [0] * n_intervals
        n_larger_all = 0
        n_larger_alt = 0

        for sample_data in record.samples.values():
            val = sample_data.get(field)
            if val is None:
                continue
            if isinstance(val, tuple):
                if not val or val[0] is None:
                    continue
                val = val[0]

            gt = sample_data.get("GT")
            alt_allele_count = 0
            if gt is not None:
                alt_allele_count = sum(1 for a in gt if a is not None and a > 0)
            is_alt = alt_allele_count > 0

            val_f = float(val)
            if val_f < bins[0]:
                continue
            elif val_f >= bins[-1]:
                n_larger_all += 1
                if is_alt:
                    n_larger_alt += 1
            else:
                b_idx = get_bin_index(val_f)
                if b_idx is not None:
                    counts_all[b_idx] += 1
                    if is_alt:
                        counts_alt[b_idx] += 1

        row = [
            record.chrom, str(record.pos), record.ref,
            ",".join(record.alts) if record.alts else ".",
            record.id if record.id else ".",
            "|".join(map(str, counts_all)),
            "|".join(map(str, counts_alt)),
        ]
        if larger_flag:
            row += [str(n_larger_all), str(n_larger_alt)]

        out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File annotation_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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

task GenerateABAnnotationTsv {
    input {
        File vcf
        File vcf_idx
        Array[Float] ab_bins
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import json
import pysam

with open("~{write_json(ab_bins)}") as f:
    bins = json.load(f)
prefix = "~{prefix}"

bin_edges_str = "|".join(map(str, bins))
n_intervals = len(bins) - 1

def get_bin_index(val):
    for j in range(n_intervals):
        if val < bins[j + 1]:
            return j
    return None

vcf_in = pysam.VariantFile("~{vcf}")

with open(f"{prefix}.tsv", "w") as out:
    out.write("\t".join(["#CHROM", "POS", "REF", "ALT", "ID", "ab_hist_alt_bin_freq"]) + "\n")
    for record in vcf_in:
        counts_alt = [0] * n_intervals

        for sample_data in record.samples.values():
            gt = sample_data.get("GT")
            if gt is None or None in gt or len(gt) != 2:
                continue
            if sorted(gt) != [0, 1]:
                continue

            ad = sample_data.get("AD")
            if ad is None or len(ad) < 2 or ad[0] is None or ad[1] is None:
                continue
            total = ad[0] + ad[1]
            if total == 0:
                continue

            ab = ad[0] / float(total)
            if ab < bins[0] or ab >= bins[-1]:
                continue
            b_idx = get_bin_index(ab)
            if b_idx is not None:
                counts_alt[b_idx] += 1

        row = [
            record.chrom, str(record.pos), record.ref,
            ",".join(record.alts) if record.alts else ".",
            record.id if record.id else ".",
            "|".join(map(str, counts_alt)),
        ]
        out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File annotation_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
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
