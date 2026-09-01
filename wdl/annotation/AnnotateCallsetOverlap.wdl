version 1.0

import "../utils/BedtoolsClosestSV.wdl"
import "../utils/ExactMatch.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "../utils/TruvariMatch.wdl"

workflow AnnotateCallsetOverlap {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        Array[String] contigs
        String prefix

        Int min_sv_length_truvari_vcf
        Int min_sv_length_truvari_truth_vcf
        Int min_sv_length_bedtools_closest_vcf
        Int min_sv_length_bedtools_closest_truth_vcf

        Int? shard_bin_size_exact_match

        String type_field_vcf = "allele_type"
        String length_field_vcf = "allele_length"
        String source_tag_truth_snv_indel_vcf = "SNV_indel"
        String source_tag_truth_sv_vcf = "SV"

        String? args_string_vcf
        String? args_string_truth_snv_indel_vcf
        String? args_string_truth_sv_vcf
        String? rename_id_string_vcf
        String? rename_id_string_truth_snv_indel_vcf
        String? rename_id_string_truth_sv_vcf
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_truth_snv_indel_vcf
        Boolean? rename_id_strip_chr_truth_sv_vcf

        File ref_fa
        File ref_fai

        String gatk_sv_lr_docker
        String utils_docker

        RuntimeAttr? runtime_attr_strip_genotypes
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_subset_sv_truth
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_rename_sv_truth
        RuntimeAttr? runtime_attr_create_exact_shards
        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_concat_exact_annotations
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_concat_truvari_eval
        RuntimeAttr? runtime_attr_concat_truvari_truth
        RuntimeAttr? runtime_attr_truvari_run_truvari_09
        RuntimeAttr? runtime_attr_truvari_run_truvari_07
        RuntimeAttr? runtime_attr_truvari_run_truvari_05
        RuntimeAttr? runtime_attr_truvari_concat_matched
        RuntimeAttr? runtime_attr_truvari_concat_matched_truth
        RuntimeAttr? runtime_attr_append_truvari_annotations
        RuntimeAttr? runtime_attr_bedtools_subset_vcf
        RuntimeAttr? runtime_attr_bedtools_subset_truth
        RuntimeAttr? runtime_attr_bedtools_convert_to_symbolic
        RuntimeAttr? runtime_attr_bedtools_split_vcf
        RuntimeAttr? runtime_attr_bedtools_split_truth
        RuntimeAttr? runtime_attr_bedtools_compare
        RuntimeAttr? runtime_attr_bedtools_calculate
        RuntimeAttr? runtime_attr_bedtools_merge_comparisons
        RuntimeAttr? runtime_attr_append_bedtools_annotations
        RuntimeAttr? runtime_attr_build_annotation_tsv
        RuntimeAttr? runtime_attr_merge_annotation_tsvs
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig || defined(args_string_vcf)) {
            call Helpers.SubsetVcfByArgs as SubsetEval {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    include_args = args_string_vcf,
                    extra_args = if single_contig then "" else "--regions ~{contig}",
                    prefix = "~{prefix}.~{contig}.eval",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }
        }

        if (!single_contig || defined(args_string_truth_snv_indel_vcf)) {
            call Helpers.SubsetVcfByArgs as SubsetTruth {
                input:
                    vcf = truth_snv_indel_vcf,
                    vcf_idx = truth_snv_indel_vcf_idx,
                    include_args = args_string_truth_snv_indel_vcf,
                    extra_args = if single_contig then "" else "--regions ~{contig}",
                    prefix = "~{prefix}.~{contig}.truth",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_truth
            }
        }

        if (!single_contig || defined(args_string_truth_sv_vcf)) {
            call Helpers.SubsetVcfByArgs as SubsetSVTruth {
                input:
                    vcf = truth_sv_vcf,
                    vcf_idx = truth_sv_vcf_idx,
                    include_args = args_string_truth_sv_vcf,
                    extra_args = if single_contig then "" else "--regions ~{contig}",
                    prefix = "~{prefix}.~{contig}.sv_truth",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_sv_truth
            }
        }

        File vcf_subsetted = select_first([SubsetEval.subset_vcf, vcf])
        File vcf_subsetted_idx = select_first([SubsetEval.subset_vcf_idx, vcf_idx])
        File truth_snv_indel_vcf_subsetted = select_first([SubsetTruth.subset_vcf, truth_snv_indel_vcf])
        File truth_snv_indel_vcf_subsetted_idx = select_first([SubsetTruth.subset_vcf_idx, truth_snv_indel_vcf_idx])
        File truth_sv_vcf_subsetted = select_first([SubsetSVTruth.subset_vcf, truth_sv_vcf])
        File truth_sv_vcf_subsetted_idx = select_first([SubsetSVTruth.subset_vcf_idx, truth_sv_vcf_idx])

        if (defined(rename_id_string_truth_sv_vcf)) {
            call Helpers.RenameVariantIds as RenameSVTruthIds {
                input:
                    vcf = truth_sv_vcf_subsetted,
                    vcf_idx = truth_sv_vcf_subsetted_idx,
                    prefix = "~{prefix}.~{contig}.sv_truth.renamed",
                    id_format = select_first([rename_id_string_truth_sv_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_truth_sv_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_sv_truth
            }
        }

        File truth_sv_vcf_final = select_first([RenameSVTruthIds.renamed_vcf, truth_sv_vcf_subsetted])
        File truth_sv_vcf_final_idx = select_first([RenameSVTruthIds.renamed_vcf_idx, truth_sv_vcf_subsetted_idx])

        call ExactMatch.ExactMatch {
            input:
                vcf = vcf_subsetted,
                vcf_idx = vcf_subsetted_idx,
                truth_snv_indel_vcf = truth_snv_indel_vcf_subsetted,
                truth_snv_indel_vcf_idx = truth_snv_indel_vcf_subsetted_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                shard_bin_size_exact_match = shard_bin_size_exact_match,
                min_sv_length_truvari_vcf = min_sv_length_truvari_vcf,
                min_sv_length_truvari_truth_vcf = min_sv_length_truvari_truth_vcf,
                length_field_vcf = length_field_vcf,
                source_tag_truth_snv_indel_vcf = source_tag_truth_snv_indel_vcf,
                rename_id_string_vcf = rename_id_string_vcf,
                rename_id_string_truth_snv_indel_vcf = rename_id_string_truth_snv_indel_vcf,
                rename_id_strip_chr_vcf = rename_id_strip_chr_vcf,
                rename_id_strip_chr_truth_snv_indel_vcf = rename_id_strip_chr_truth_snv_indel_vcf,
                utils_docker = utils_docker,
                runtime_attr_strip_genotypes = runtime_attr_strip_genotypes,
                runtime_attr_create_exact_shards = runtime_attr_create_exact_shards,
                runtime_attr_subset_exact_vcf = runtime_attr_subset_exact_vcf,
                runtime_attr_subset_exact_truth = runtime_attr_subset_exact_truth,
                runtime_attr_rename_vcf = runtime_attr_rename_vcf,
                runtime_attr_rename_truth = runtime_attr_rename_truth,
                runtime_attr_exact_match = runtime_attr_exact_match,
                runtime_attr_concat_exact_annotations = runtime_attr_concat_exact_annotations,
                runtime_attr_append_exact_annotations = runtime_attr_append_exact_annotations,
                runtime_attr_truvari_subset_vcf = runtime_attr_truvari_subset_vcf,
                runtime_attr_truvari_subset_truth = runtime_attr_truvari_subset_truth,
                runtime_attr_concat_truvari_eval = runtime_attr_concat_truvari_eval,
                runtime_attr_concat_truvari_truth = runtime_attr_concat_truvari_truth
        }

        call TruvariMatch.TruvariMatch {
            input:
                vcf = ExactMatch.truvari_eval_vcf,
                vcf_idx = ExactMatch.truvari_eval_vcf_idx,
                truth_snv_indel_vcf = ExactMatch.truvari_truth_vcf,
                truth_snv_indel_vcf_idx = ExactMatch.truvari_truth_vcf_idx,
                prefix = "~{prefix}.~{contig}.truvari",
                source_tag = source_tag_truth_snv_indel_vcf,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                utils_docker = utils_docker,
                runtime_attr_run_truvari_09 = runtime_attr_truvari_run_truvari_09,
                runtime_attr_run_truvari_07 = runtime_attr_truvari_run_truvari_07,
                runtime_attr_run_truvari_05 = runtime_attr_truvari_run_truvari_05,
                runtime_attr_concat_matched = runtime_attr_truvari_concat_matched,
                runtime_attr_concat_matched_truth = runtime_attr_truvari_concat_matched_truth
        }

        call Helpers.AppendAnnotationsFromVcf as AppendTruvariAnnotations {
            input:
                annotation_tsv = TruvariMatch.annotation_tsv,
                truth_vcf = TruvariMatch.matched_truth_vcf,
                truth_vcf_idx = TruvariMatch.matched_truth_vcf_idx,
                is_sv_truth = false,
                prefix = "~{prefix}.~{contig}.truvari_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_append_truvari_annotations
        }

        call BedtoolsClosestSV.BedtoolsClosestSV {
            input:
                vcf = TruvariMatch.unmatched_vcf,
                vcf_idx = TruvariMatch.unmatched_vcf_idx,
                truth_sv_vcf = truth_sv_vcf_final,
                truth_sv_vcf_idx = truth_sv_vcf_final_idx,
                prefix = "~{prefix}.~{contig}.bedtools_closest",
                min_sv_length = min_sv_length_bedtools_closest_vcf,
                min_sv_length_truth = min_sv_length_bedtools_closest_truth_vcf,
                type_field = type_field_vcf,
                length_field = length_field_vcf,
                source_tag = source_tag_truth_sv_vcf,
                gatk_sv_lr_docker = gatk_sv_lr_docker,
                utils_docker = utils_docker,
                runtime_attr_subset_vcf = runtime_attr_bedtools_subset_vcf,
                runtime_attr_subset_truth = runtime_attr_bedtools_subset_truth,
                runtime_attr_convert_to_symbolic = runtime_attr_bedtools_convert_to_symbolic,
                runtime_attr_split_vcf = runtime_attr_bedtools_split_vcf,
                runtime_attr_split_truth = runtime_attr_bedtools_split_truth,
                runtime_attr_compare = runtime_attr_bedtools_compare,
                runtime_attr_calculate = runtime_attr_bedtools_calculate,
                runtime_attr_merge_comparisons = runtime_attr_bedtools_merge_comparisons
        }

        call Helpers.AppendAnnotationsFromVcf as AppendBedtoolsAnnotations {
            input:
                annotation_tsv = BedtoolsClosestSV.annotation_tsv,
                truth_vcf = truth_sv_vcf_final,
                truth_vcf_idx = truth_sv_vcf_final_idx,
                is_sv_truth = true,
                prefix = "~{prefix}.~{contig}.bedtools_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_append_bedtools_annotations
        }

        Array[File] extended_annotation_tsvs = [
            ExactMatch.annotated_tsv,
            AppendTruvariAnnotations.annotated_tsv,
            AppendBedtoolsAnnotations.annotated_tsv
        ]

        call BuildBenchmarkAnnotationTsv {
            input:
                tsvs = extended_annotation_tsvs,
                prefix = "~{prefix}.~{contig}.benchmark_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_build_annotation_tsv
        }
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotationTsvs {
            input:
                tsvs = BuildBenchmarkAnnotationTsv.merged_tsv,
                sort_output = false,
                preserve_header = false,
                prefix = "~{prefix}.benchmark_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_annotation_tsvs
        }
    }

    output {
        File annotations_tsv_benchmark = select_first([MergeAnnotationTsvs.concatenated_tsv, BuildBenchmarkAnnotationTsv.merged_tsv[0]])
        File annotations_header_benchmark = BuildBenchmarkAnnotationTsv.merged_header[0]
    }
}

task BuildBenchmarkAnnotationTsv {
    input {
        Array[File] tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import re

input_files = "~{sep=',' tsvs}".split(',')
prefix = "~{prefix}"

fixed_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'ID']
static_extra = ['match_type', 'truth_ID', 'source_tag', 'filter']
genotype_cols = ['N_HOMREF', 'N_HET', 'N_HOMALT']
skip_cols = set(static_extra + genotype_cols)

# Collect AC_/AF_/AN_ field names from all TSV headers
all_ac, all_af, all_an = set(), set(), set()
for f in input_files:
    with open(f) as fh:
        header = fh.readline().strip().split('\t')
    for col in header[5:]:
        if col in skip_cols:
            continue
        if col == 'AC' or col.startswith('AC_'):
            all_ac.add(col)
        elif col == 'AF' or col.startswith('AF_'):
            all_af.add(col)
        elif col == 'AN' or col.startswith('AN_'):
            all_an.add(col)

dyn_cols = sorted(all_ac) + sorted(all_af) + sorted(all_an)
all_extra = static_extra + dyn_cols + genotype_cols
master_header = fixed_cols + all_extra

with open(f"{prefix}.tsv", 'w') as fout:
    for f in input_files:
        with open(f) as fh:
            file_cols = fh.readline().strip().split('\t')
            col_map = {name: i for i, name in enumerate(file_cols)}
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                row = []
                for col in master_header:
                    if col in col_map:
                        try:
                            row.append(parts[col_map[col]])
                        except IndexError:
                            row.append('.')
                    else:
                        row.append('.')
                fout.write('\t'.join(row) + '\n')

with open(f"{prefix}.header.txt", 'w') as hout:
    for col in all_extra:
        hout.write(col + '\n')

EOF
    >>>

    output {
        File merged_tsv = "~{prefix}.tsv"
        File merged_header = "~{prefix}.header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
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
