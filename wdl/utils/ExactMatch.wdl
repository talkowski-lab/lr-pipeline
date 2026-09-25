version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow ExactMatch {
    meta {
        description: [
            "This sub-workflow performs the first callset-comparison round, matching records to a truth callset on exact position and allele. The contig is cut into fixed-width regions and each is processed independently: both VCFs are streamed down to the region straight from the bucket, with genotypes dropped and any include expression applied on the way, then optionally renamed to a common ID scheme and matched. Each region's unmatched callset records and its truth records are then subset to the caller's minimum Truvari lengths, so no step ever passes over a whole contig. The per-region annotations and VCFs are concatenated into the form `TruvariMatch` expects."
        ]
    }

    parameter_meta {
        vcf: "Callset being compared. Read region by region, so it may hold any set of contigs."
        vcf_idx: "Index for vcf."
        truth_snv_indel_vcf: "Truth callset. Read region by region, so it may hold any set of contigs."
        truth_snv_indel_vcf_idx: "Index for truth_snv_indel_vcf."
        contig: "Contig being processed."
        shard_bin_size_exact_match: "Width in base pairs of the contig regions the matching step is sharded into."
        min_sv_length_truvari_vcf: "Minimum length for an unmatched callset record to be emitted for Truvari, measured by `length_field_vcf`."
        min_sv_length_truvari_truth_snv_indel_vcf: "Minimum length for a truth record to be emitted for Truvari, measured by `length_field_truth_snv_indel_vcf`."
        length_field_vcf: "INFO field in the callset holding allele length."
        length_field_truth_snv_indel_vcf: "Length used to filter the truth callset, either an INFO field or `ILEN`, the bcftools built-in indel length computed from REF and ALT."
        source_tag_truth_snv_indel_vcf: "Tag identifying the truth callset in the annotations."
        args_string_vcf: "`bcftools view` include expression applied to the callset as each region is read."
        args_string_truth_snv_indel_vcf: "`bcftools view` include expression applied to the truth callset as each region is read."
        rename_id_string_vcf: "ID rename templates."
        rename_id_string_truth_snv_indel_vcf: "ID rename templates."
        rename_id_strip_chr_vcf: "Strip the `chr` prefix while renaming."
        rename_id_strip_chr_truth_snv_indel_vcf: "Strip the `chr` prefix while renaming."
        annotated_tsv: "Exact-match annotations."
        truvari_eval_vcf: "Unmatched callset records at or above `min_sv_length_truvari_vcf`, passed to `TruvariMatch`."
        truvari_eval_vcf_idx: "Index for truvari_eval_vcf."
        truvari_truth_vcf: "Renamed truth records at or above `min_sv_length_truvari_truth_snv_indel_vcf`, passed to `TruvariMatch`."
        truvari_truth_vcf_idx: "Index for truvari_truth_vcf."
    }

    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String contig
        String prefix

        Int shard_bin_size_exact_match = 5000000

        Int min_sv_length_truvari_vcf
        Int min_sv_length_truvari_truth_snv_indel_vcf
        String length_field_vcf
        String length_field_truth_snv_indel_vcf
        String source_tag_truth_snv_indel_vcf

        String? args_string_vcf
        String? args_string_truth_snv_indel_vcf
        String? rename_id_string_vcf
        String? rename_id_string_truth_snv_indel_vcf
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_truth_snv_indel_vcf

        String utils_docker

        RuntimeAttr? runtime_attr_create_exact_shards
        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_concat_exact_annotations
        RuntimeAttr? runtime_attr_concat_exact_unmatched
        RuntimeAttr? runtime_attr_concat_exact_truth
    }

    call Helpers.CreateContigShards as CreateExactShards {
        input:
            vcfs = [vcf, truth_snv_indel_vcf],
            vcf_idxs = [vcf_idx, truth_snv_indel_vcf_idx],
            contig = contig,
            shard_bin_size = shard_bin_size_exact_match,
            prefix = "~{prefix}.exact_shards",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_create_exact_shards
    }

    scatter (k in range(length(CreateExactShards.shard_regions))) {
        # Read each region straight out of the bucket without genotypes, applying any include expression on the way
        call Helpers.SubsetVcfToRegionStreaming as SubsetExactEval {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                region = CreateExactShards.shard_regions[k],
                include_args = args_string_vcf,
                drop_genotypes = true,
                prefix = "~{prefix}.exact_eval_~{k}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_exact_vcf
        }

        call Helpers.SubsetVcfToRegionStreaming as SubsetExactTruth {
            input:
                vcf = truth_snv_indel_vcf,
                vcf_idx = truth_snv_indel_vcf_idx,
                region = CreateExactShards.shard_regions[k],
                include_args = args_string_truth_snv_indel_vcf,
                drop_genotypes = true,
                prefix = "~{prefix}.exact_truth_~{k}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_exact_truth
        }

        if (defined(rename_id_string_vcf)) {
            call Helpers.RenameVariantIds as RenameEvalIds {
                input:
                    vcf = SubsetExactEval.subset_vcf,
                    vcf_idx = SubsetExactEval.subset_vcf_idx,
                    prefix = "~{prefix}.exact_eval_~{k}.renamed",
                    id_format = select_first([rename_id_string_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_vcf
            }
        }

        if (defined(rename_id_string_truth_snv_indel_vcf)) {
            call Helpers.RenameVariantIds as RenameTruthIds {
                input:
                    vcf = SubsetExactTruth.subset_vcf,
                    vcf_idx = SubsetExactTruth.subset_vcf_idx,
                    prefix = "~{prefix}.exact_truth_~{k}.renamed",
                    id_format = select_first([rename_id_string_truth_snv_indel_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_truth_snv_indel_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_truth
            }
        }

        File shard_eval_vcf = select_first([RenameEvalIds.renamed_vcf, SubsetExactEval.subset_vcf])
        File shard_eval_vcf_idx = select_first([RenameEvalIds.renamed_vcf_idx, SubsetExactEval.subset_vcf_idx])
        File shard_truth_vcf = select_first([RenameTruthIds.renamed_vcf, SubsetExactTruth.subset_vcf])
        File shard_truth_vcf_idx = select_first([RenameTruthIds.renamed_vcf_idx, SubsetExactTruth.subset_vcf_idx])

        call Helpers.ExactMatch as ExactMatchShard {
            input:
                vcf = shard_eval_vcf,
                vcf_idx = shard_eval_vcf_idx,
                truth_snv_indel_vcf = shard_truth_vcf,
                truth_snv_indel_vcf_idx = shard_truth_vcf_idx,
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

        # Admit only records long enough for Truvari, measuring each VCF by its own length field
        call Helpers.SubsetVcfByLength as SubsetTruvariEval {
            input:
                vcf = ExactMatchShard.unmatched_vcf,
                vcf_idx = ExactMatchShard.unmatched_vcf_idx,
                length_field = length_field_vcf,
                min_length = min_sv_length_truvari_vcf,
                prefix = "~{prefix}.truvari_eval_~{k}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_truvari_subset_vcf
        }

        call Helpers.SubsetVcfByLength as SubsetTruvariTruth {
            input:
                vcf = shard_truth_vcf,
                vcf_idx = shard_truth_vcf_idx,
                length_field = length_field_truth_snv_indel_vcf,
                min_length = min_sv_length_truvari_truth_snv_indel_vcf,
                prefix = "~{prefix}.truvari_truth_~{k}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_truvari_subset_truth
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
            vcfs = SubsetTruvariEval.subset_vcf,
            vcf_idxs = SubsetTruvariEval.subset_vcf_idx,
            allow_overlaps = false,
            naive = false,
            prefix = "~{prefix}.truvari_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_exact_unmatched
    }

    call Helpers.ConcatVcfs as ConcatExactTruth {
        input:
            vcfs = SubsetTruvariTruth.subset_vcf,
            vcf_idxs = SubsetTruvariTruth.subset_vcf_idx,
            allow_overlaps = false,
            naive = false,
            prefix = "~{prefix}.truvari_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_exact_truth
    }

    output {
        File annotated_tsv = ConcatExactAnnotations.concatenated_tsv
        File truvari_eval_vcf = ConcatExactUnmatched.concat_vcf
        File truvari_eval_vcf_idx = ConcatExactUnmatched.concat_vcf_idx
        File truvari_truth_vcf = ConcatExactTruth.concat_vcf
        File truvari_truth_vcf_idx = ConcatExactTruth.concat_vcf_idx
    }
}
