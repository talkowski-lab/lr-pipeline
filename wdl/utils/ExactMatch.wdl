version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExactMatch {
    meta {
        description: [
            "This sub-workflow performs the first callset-comparison round, matching records to a truth callset on exact position and allele. Both callsets are optionally renamed to a common ID scheme, sharded, matched, and the annotations concatenated. Records left unmatched are emitted alongside the renamed truth callset, both unfiltered by length, for the caller to subset before `TruvariMatch`."
        ]
    }

    parameter_meta {
        vcf: "Callset being compared."
        vcf_idx: "Index for vcf."
        truth_snv_indel_vcf: "Truth callset."
        truth_snv_indel_vcf_idx: "Index for truth_snv_indel_vcf."
        contig: "Contig being processed."
        shard_bin_size_exact_match: "Shard size for the matching step."
        source_tag_truth_snv_indel_vcf: "Tag identifying the truth callset in the annotations."
        rename_id_string_vcf: "ID rename templates."
        rename_id_string_truth_snv_indel_vcf: "ID rename templates."
        rename_id_strip_chr_vcf: "Strip the `chr` prefix while renaming."
        rename_id_strip_chr_truth_snv_indel_vcf: "Strip the `chr` prefix while renaming."
        annotated_tsv: "Exact-match annotations."
        unmatched_vcf: "Callset records left unmatched, not yet filtered by length."
        unmatched_vcf_idx: "Index for unmatched_vcf."
        truth_vcf: "Truth callset after optional ID renaming, not yet filtered by length."
        truth_vcf_idx: "Index for truth_vcf."
    }

    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String contig
        String prefix

        Int? shard_bin_size_exact_match

        String source_tag_truth_snv_indel_vcf

        String? rename_id_string_vcf
        String? rename_id_string_truth_snv_indel_vcf
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_truth_snv_indel_vcf

        String utils_docker

        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_create_exact_shards
        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_concat_exact_annotations
        RuntimeAttr? runtime_attr_concat_exact_unmatched
    }

    if (defined(rename_id_string_vcf)) {
        call Helpers.RenameVariantIds as RenameEvalIds {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                prefix = "~{prefix}.eval.renamed",
                id_format = select_first([rename_id_string_vcf]),
                strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_vcf
        }
    }

    if (defined(rename_id_string_truth_snv_indel_vcf)) {
        call Helpers.RenameVariantIds as RenameTruthIds {
            input:
                vcf = truth_snv_indel_vcf,
                vcf_idx = truth_snv_indel_vcf_idx,
                prefix = "~{prefix}.truth.renamed",
                id_format = select_first([rename_id_string_truth_snv_indel_vcf]),
                strip_chr = select_first([rename_id_strip_chr_truth_snv_indel_vcf, false]),
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_truth
        }
    }

    File eval_vcf_final = select_first([RenameEvalIds.renamed_vcf, vcf])
    File eval_vcf_final_idx = select_first([RenameEvalIds.renamed_vcf_idx, vcf_idx])
    File truth_vcf_final = select_first([RenameTruthIds.renamed_vcf, truth_snv_indel_vcf])
    File truth_vcf_final_idx = select_first([RenameTruthIds.renamed_vcf_idx, truth_snv_indel_vcf_idx])

    if (defined(shard_bin_size_exact_match)) {
        call Helpers.CreateContigShards as CreateExactShards {
            input:
                vcfs = [eval_vcf_final, truth_vcf_final],
                vcf_idxs = [eval_vcf_final_idx, truth_vcf_final_idx],
                contig = contig,
                shard_bin_size = select_first([shard_bin_size_exact_match]),
                prefix = "~{prefix}.exact_shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_exact_shards
        }

        scatter (k in range(length(CreateExactShards.shard_regions))) {
            call Helpers.SubsetVcfToRegionStreaming as SubsetExactEval {
                input:
                    vcf = eval_vcf_final,
                    vcf_idx = eval_vcf_final_idx,
                    region = CreateExactShards.shard_regions[k],
                    prefix = "~{prefix}.exact_eval_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_exact_vcf
            }

            call Helpers.SubsetVcfToRegionStreaming as SubsetExactTruth {
                input:
                    vcf = truth_vcf_final,
                    vcf_idx = truth_vcf_final_idx,
                    region = CreateExactShards.shard_regions[k],
                    prefix = "~{prefix}.exact_truth_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_exact_truth
            }

            call Helpers.ExactMatch as ExactMatchShard {
                input:
                    vcf = SubsetExactEval.subset_vcf,
                    vcf_idx = SubsetExactEval.subset_vcf_idx,
                    truth_snv_indel_vcf = SubsetExactTruth.subset_vcf,
                    truth_snv_indel_vcf_idx = SubsetExactTruth.subset_vcf_idx,
                    source_tag = source_tag_truth_snv_indel_vcf,
                    prefix = "~{prefix}.exact_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_exact_match
            }

            call Helpers.AppendAnnotationsFromVcf as AppendExactAnnotationsShard {
                input:
                    annotation_tsv = ExactMatchShard.annotation_tsv,
                    truth_vcf = ExactMatchShard.matched_truth_vcf,
                    truth_vcf_idx = ExactMatchShard.matched_truth_vcf_idx,
                    is_sv_truth = false,
                    prefix = "~{prefix}.exact_annotated_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_append_exact_annotations
            }
        }

        call Helpers.ConcatTsvs as ConcatExactAnnotations {
            input:
                tsvs = AppendExactAnnotationsShard.annotated_tsv,
                sort_output = true,
                preserve_header = true,
                prefix = "~{prefix}.exact_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_exact_annotations
        }

        call Helpers.ConcatVcfs as ConcatExactUnmatched {
            input:
                vcfs = ExactMatchShard.unmatched_vcf,
                vcf_idxs = ExactMatchShard.unmatched_vcf_idx,
                allow_overlaps = false,
                naive = false,
                prefix = "~{prefix}.exact_unmatched",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_exact_unmatched
        }
    }

    if (!defined(shard_bin_size_exact_match)) {
        call Helpers.ExactMatch as ExactMatchFull {
            input:
                vcf = eval_vcf_final,
                vcf_idx = eval_vcf_final_idx,
                truth_snv_indel_vcf = truth_vcf_final,
                truth_snv_indel_vcf_idx = truth_vcf_final_idx,
                source_tag = source_tag_truth_snv_indel_vcf,
                prefix = "~{prefix}.exact",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_exact_match
        }

        call Helpers.AppendAnnotationsFromVcf as AppendExactAnnotationsFull {
            input:
                annotation_tsv = ExactMatchFull.annotation_tsv,
                truth_vcf = ExactMatchFull.matched_truth_vcf,
                truth_vcf_idx = ExactMatchFull.matched_truth_vcf_idx,
                is_sv_truth = false,
                prefix = "~{prefix}.exact_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_append_exact_annotations
        }
    }

    File annotated_tsv_final = select_first([ConcatExactAnnotations.concatenated_tsv, AppendExactAnnotationsFull.annotated_tsv])
    File unmatched_vcf_final = select_first([ConcatExactUnmatched.concat_vcf, ExactMatchFull.unmatched_vcf])
    File unmatched_vcf_final_idx = select_first([ConcatExactUnmatched.concat_vcf_idx, ExactMatchFull.unmatched_vcf_idx])

    output {
        File annotated_tsv = annotated_tsv_final
        File unmatched_vcf = unmatched_vcf_final
        File unmatched_vcf_idx = unmatched_vcf_final_idx
        File truth_vcf = truth_vcf_final
        File truth_vcf_idx = truth_vcf_final_idx
    }
}
