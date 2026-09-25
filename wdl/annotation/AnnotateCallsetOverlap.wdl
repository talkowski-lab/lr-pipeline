version 1.0

import "../utils/BedtoolsClosestSV.wdl"
import "../utils/ExactMatch.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "../utils/TruvariMatch.wdl"

workflow AnnotateCallsetOverlap {
    meta {
        description: [
            "This workflow ingests a callset VCF and two truth VCFs - one of SNVs & indels and one of SVs - and finds matching variants across them, annotating each matched callset variant with the truth callset's AC/AF/AN and genotype-count fields. This enables benchmarking annotations against an existing cohort (e.g. gnomAD) and surfacing variants that are outliers relative to it.",
            "The workflow undergoes multiple rounds of variant matching in order to determine matched pairs: (1) Exact match across CHROM, POS, REF and ALT. (2) Truvari match with overlap percentages of 90%, 70% and 50%. (3) Matching based on `bedtools closest`, finetuned for SVs. Here the callset and truth variants are split by type and converted to a symbolic representation, after which separate `bedtools closest` passes are run - one tuned for deletions and duplications via reciprocal positional overlap, and one tuned for insertions via breakpoint proximity - so that each callset variant is paired with the nearest same-type truth variant above the per-callset minimum SV-length thresholds.",
            "Note: When converting to symbolic representation, only canonical DUPs (allele_type = `DUP` exactly) are treated as DUP; other DUP subtypes (e.g., `dup_interspersed`, `inv_dup`) are treated as insertions.",
            "Note: Callset DUPs are compared twice, because the two matching rules need different coordinates. Against truth DUPs they are compared by reciprocal overlap, over their `ORIGIN` interval when `move_dup_to_origin` is true and over their own span from POS otherwise; against truth insertions they are always collapsed to a point at their VCF position and compared by breakpoint proximity and length ratio. Setting `move_dup_to_origin` to false is what lets the workflow run on a callset with no `INFO/ORIGIN`.",
            "Note: The SV truth VCF is expected to be symbolic already. Set `convert_symbolic_truth_sv_vcf` when it instead carries sequence alleles, in which case its allele type and length are read from `type_field_truth_sv_vcf` and `length_field_truth_sv_vcf`, which need not match the fields used for the callset. Its canonical DUPs then follow the same `move_dup_to_origin` positioning as the callset.",
            "Every minimum-length filter is applied in this workflow, measuring each VCF by its own `length_field_*` input, and the SV truth VCF is filtered before renaming and conversion. With `convert_symbolic_truth_sv_vcf` a canonical DUP is therefore measured by `length_field_truth_sv_vcf` rather than by the `ORIGIN` span that conversion writes into `SVLEN`.",
            "The workflow runs on one contig. Each VCF is either streamed down to that contig straight from the bucket or, when its `subset_contig_*` input is false, taken as already containing only that contig. Its optional `args_string_*` expression is then applied in a separate pass that also drops genotypes, which dominate localization but go unused here.",
            "Both the exact-match and Truvari rounds can be sharded within a contig. Truvari shard boundaries are snapped forward to the next gap wider than the `min_shard_gap_truvari_match` input of `TruvariMatch`, which keeps results identical to an unsharded run because Truvari only groups records into a new comparison chunk once the next record clears the running end by more than its chunk size. Fixed-width bins alone would split colocated record pairs and silently lose matches."
        ]
    }

    parameter_meta {
        vcf: "Callset VCF being annotated."
        vcf_idx: "Index for `vcf`."
        truth_snv_indel_vcf: "Truth VCF containing SNVs & indels to match against."
        truth_snv_indel_vcf_idx: "Index for `truth_snv_indel_vcf`."
        truth_sv_vcf: "Truth VCF containing SVs to match against."
        truth_sv_vcf_idx: "Index for `truth_sv_vcf`."
        contig: "Contig to evaluate."
        min_sv_length_truvari_vcf: "Minimum length for a callset variant to enter the Truvari matching round, measured by `length_field_vcf`."
        min_sv_length_truvari_truth_snv_indel_vcf: "Minimum length for a SNV & indel truth variant to enter the Truvari matching round, measured by `length_field_truth_snv_indel_vcf`."
        min_sv_length_bedtools_closest_vcf: "Minimum length for a callset variant to enter the `bedtools closest` matching round, measured by `length_field_vcf`."
        min_sv_length_bedtools_closest_truth_vcf: "Minimum length for an SV truth variant to enter the `bedtools closest` matching round, measured by `length_field_truth_sv_vcf` before any renaming or conversion."
        shard_bin_size_exact_match: "If set, shards the exact-match round into contig regions of roughly this many base pairs, run in parallel."
        shard_bin_size_truvari_match: "If set, shards the Truvari round into contig regions of at least this many base pairs, run in parallel. Each region is extended to the next safe gap, so a value of 1000000 or more is recommended."
        subset_contig_vcf: "Whether to stream `vcf` down to `contig` first. When false it must already contain only that contig."
        subset_contig_truth_snv_indel_vcf: "Whether to stream `truth_snv_indel_vcf` down to `contig` first. When false it must already contain only that contig."
        subset_contig_truth_sv_vcf: "Whether to stream `truth_sv_vcf` down to `contig` first. When false it must already contain only that contig."
        convert_symbolic_truth_sv_vcf: "Whether the SV truth VCF represents alleles as sequence rather than symbolically. When true it is converted to a symbolic representation first, reading the `type_field_truth_sv_vcf` and `length_field_truth_sv_vcf` INFO fields."
        move_dup_to_origin: "Whether canonical DUPs are repositioned onto their `INFO/ORIGIN` interval before the DUP-vs-DUP reciprocal-overlap comparison. When false each DUP instead spans its own coordinates, from POS over its allele length, and `INFO/ORIGIN` is not required."
        type_field_vcf: "INFO field in the callset VCF giving each variant's allele type."
        type_field_truth_sv_vcf: "INFO field in the SV truth VCF giving each variant's allele type, read only when `convert_symbolic_truth_sv_vcf` is true. A symbolic truth VCF carries `SVTYPE`; a sequence-allele one needs its own field named here, e.g. `allele_type`."
        length_field_vcf: "INFO field in the callset VCF giving each variant's allele length, read by both callset length filters and by the symbolic conversion."
        length_field_truth_snv_indel_vcf: "Length used to filter the SNV & indel truth VCF, either an INFO field or `ILEN`, the bcftools built-in indel length computed from REF and ALT for a truth VCF such as gnomAD that carries no length field."
        length_field_truth_sv_vcf: "INFO field in the SV truth VCF giving each variant's allele length, used to apply `min_sv_length_bedtools_closest_truth_vcf` and, when `convert_symbolic_truth_sv_vcf` is true, read by the conversion. A symbolic truth VCF carries `SVLEN`; a sequence-allele one needs its own field named here, e.g. `allele_length`."
        source_tag_truth_snv_indel_vcf: "Label used to tag matches against the SNV & indel truth VCF."
        source_tag_truth_sv_vcf: "Label used to tag matches against the SV truth VCF."
        args_string_vcf: "`bcftools view` include expression used to pre-subset the callset VCF."
        args_string_truth_snv_indel_vcf: "`bcftools view` include expression used to pre-subset the SNV & indel truth VCF."
        args_string_truth_sv_vcf: "`bcftools view` include expression used to pre-subset the SV truth VCF."
        rename_id_string_vcf: "Expression used to rename variant IDs in the callset VCF prior to matching."
        rename_id_string_truth_snv_indel_vcf: "Expression used to rename variant IDs in the SNV & indel truth VCF prior to matching."
        rename_id_string_truth_sv_vcf: "Expression used to rename variant IDs in the SV truth VCF prior to matching."
        rename_id_strip_chr_vcf: "Whether to strip the `chr` prefix when renaming callset variant IDs."
        rename_id_strip_chr_truth_snv_indel_vcf: "Whether to strip the `chr` prefix when renaming SNV & indel truth variant IDs."
        rename_id_strip_chr_truth_sv_vcf: "Whether to strip the `chr` prefix when renaming SV truth variant IDs."
        ref_fa: "From references. Only needed when either VCF represents alleles symbolically, since Truvari uses it solely to resolve those alleles to sequence."
        ref_fai: "From references."
        annotations_tsv_benchmark: "TSV mapping callset variants to their matched truth variants, match type, and the truth callset's AC/AF/AN and genotype-count fields."
        annotations_header_benchmark: "Header listing the extra annotation columns present in `annotations_tsv_benchmark`."
    }

    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        String contig
        String prefix

        Int min_sv_length_truvari_vcf
        Int min_sv_length_truvari_truth_snv_indel_vcf
        Int min_sv_length_bedtools_closest_vcf
        Int min_sv_length_bedtools_closest_truth_vcf

        Int? shard_bin_size_exact_match
        Int? shard_bin_size_truvari_match

        Boolean subset_contig_vcf = true
        Boolean subset_contig_truth_snv_indel_vcf = true
        Boolean subset_contig_truth_sv_vcf = true
        Boolean convert_symbolic_truth_sv_vcf = false
        Boolean move_dup_to_origin = true

        String type_field_vcf = "allele_type"
        String type_field_truth_sv_vcf = "SVTYPE"
        String length_field_vcf = "allele_length"
        String length_field_truth_snv_indel_vcf = "ILEN"
        String length_field_truth_sv_vcf = "SVLEN"
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

        File? ref_fa
        File? ref_fai

        String gatk_sv_lr_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_contig_vcf
        RuntimeAttr? runtime_attr_subset_contig_truth
        RuntimeAttr? runtime_attr_subset_contig_sv_truth
        RuntimeAttr? runtime_attr_subset_args_vcf
        RuntimeAttr? runtime_attr_subset_args_truth
        RuntimeAttr? runtime_attr_subset_args_sv_truth
        RuntimeAttr? runtime_attr_rename_sv_truth
        RuntimeAttr? runtime_attr_convert_sv_truth
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_create_exact_shards
        RuntimeAttr? runtime_attr_subset_exact_vcf
        RuntimeAttr? runtime_attr_subset_exact_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_concat_exact_annotations
        RuntimeAttr? runtime_attr_concat_exact_unmatched
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_truvari_create_shards
        RuntimeAttr? runtime_attr_truvari_subset_region_vcf
        RuntimeAttr? runtime_attr_truvari_subset_region_truth
        RuntimeAttr? runtime_attr_truvari_run_truvari_09
        RuntimeAttr? runtime_attr_truvari_run_truvari_07
        RuntimeAttr? runtime_attr_truvari_run_truvari_05
        RuntimeAttr? runtime_attr_truvari_concat_matched
        RuntimeAttr? runtime_attr_truvari_concat_matched_truth
        RuntimeAttr? runtime_attr_truvari_concat_unmatched
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
    }

    # Stream each VCF down to the contig straight out of the bucket, unless it already holds only that contig
    if (subset_contig_vcf) {
        call Helpers.SubsetVcfToRegionStreaming as SubsetContigEval {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                region = contig,
                drop_genotypes = true,
                prefix = "~{prefix}.~{contig}.eval.contig",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_vcf
        }
    }

    if (subset_contig_truth_snv_indel_vcf) {
        call Helpers.SubsetVcfToRegionStreaming as SubsetContigTruth {
            input:
                vcf = truth_snv_indel_vcf,
                vcf_idx = truth_snv_indel_vcf_idx,
                region = contig,
                drop_genotypes = true,
                prefix = "~{prefix}.~{contig}.truth.contig",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_truth
        }
    }

    if (subset_contig_truth_sv_vcf) {
        call Helpers.SubsetVcfToRegionStreaming as SubsetContigSVTruth {
            input:
                vcf = truth_sv_vcf,
                vcf_idx = truth_sv_vcf_idx,
                region = contig,
                drop_genotypes = true,
                prefix = "~{prefix}.~{contig}.sv_truth.contig",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_contig_sv_truth
        }
    }

    # Apply each optional include expression in its own pass, dropping genotypes so a VCF that skipped the contig subset arrives sites-only
    call Helpers.SubsetVcfByArgs as SubsetEval {
        input:
            vcf = select_first([SubsetContigEval.subset_vcf, vcf]),
            vcf_idx = select_first([SubsetContigEval.subset_vcf_idx, vcf_idx]),
            include_args = args_string_vcf,
            extra_args = "-G",
            prefix = "~{prefix}.~{contig}.eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_args_vcf
    }

    call Helpers.SubsetVcfByArgs as SubsetTruth {
        input:
            vcf = select_first([SubsetContigTruth.subset_vcf, truth_snv_indel_vcf]),
            vcf_idx = select_first([SubsetContigTruth.subset_vcf_idx, truth_snv_indel_vcf_idx]),
            include_args = args_string_truth_snv_indel_vcf,
            extra_args = "-G",
            prefix = "~{prefix}.~{contig}.truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_args_truth
    }

    call Helpers.SubsetVcfByArgs as SubsetSVTruth {
        input:
            vcf = select_first([SubsetContigSVTruth.subset_vcf, truth_sv_vcf]),
            vcf_idx = select_first([SubsetContigSVTruth.subset_vcf_idx, truth_sv_vcf_idx]),
            include_args = args_string_truth_sv_vcf,
            extra_args = "-G",
            prefix = "~{prefix}.~{contig}.sv_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_args_sv_truth
    }

    # Apply the SV truth length filter before renaming and conversion, so every later step sees only records long enough
    call Helpers.SubsetVcfByLength as SubsetBedtoolsTruth {
        input:
            vcf = SubsetSVTruth.subset_vcf,
            vcf_idx = SubsetSVTruth.subset_vcf_idx,
            length_field = length_field_truth_sv_vcf,
            min_length = min_sv_length_bedtools_closest_truth_vcf,
            prefix = "~{prefix}.~{contig}.sv_truth.subset",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_bedtools_subset_truth
    }

    if (defined(rename_id_string_truth_sv_vcf)) {
        call Helpers.RenameVariantIds as RenameSVTruthIds {
            input:
                vcf = SubsetBedtoolsTruth.subset_vcf,
                vcf_idx = SubsetBedtoolsTruth.subset_vcf_idx,
                prefix = "~{prefix}.~{contig}.sv_truth.renamed",
                id_format = select_first([rename_id_string_truth_sv_vcf]),
                strip_chr = select_first([rename_id_strip_chr_truth_sv_vcf, false]),
                docker = utils_docker,
                runtime_attr_override = runtime_attr_rename_sv_truth
        }
    }

    File truth_sv_vcf_renamed = select_first([RenameSVTruthIds.renamed_vcf, SubsetBedtoolsTruth.subset_vcf])
    File truth_sv_vcf_renamed_idx = select_first([RenameSVTruthIds.renamed_vcf_idx, SubsetBedtoolsTruth.subset_vcf_idx])

    # Give a sequence-allele SV truth VCF the symbolic ALTs, SVTYPE, SVLEN and END the bedtools closest round reads
    if (convert_symbolic_truth_sv_vcf) {
        call Helpers.ConvertToSymbolic as ConvertSVTruth {
            input:
                vcf = truth_sv_vcf_renamed,
                vcf_idx = truth_sv_vcf_renamed_idx,
                move_dup_to_origin = move_dup_to_origin,
                type_field = type_field_truth_sv_vcf,
                length_field = length_field_truth_sv_vcf,
                prefix = "~{prefix}.~{contig}.sv_truth.symbolic",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert_sv_truth
        }
    }

    File truth_sv_vcf_final = select_first([ConvertSVTruth.processed_vcf, truth_sv_vcf_renamed])
    File truth_sv_vcf_final_idx = select_first([ConvertSVTruth.processed_vcf_idx, truth_sv_vcf_renamed_idx])

    call ExactMatch.ExactMatch {
        input:
            vcf = SubsetEval.subset_vcf,
            vcf_idx = SubsetEval.subset_vcf_idx,
            truth_snv_indel_vcf = SubsetTruth.subset_vcf,
            truth_snv_indel_vcf_idx = SubsetTruth.subset_vcf_idx,
            contig = contig,
            prefix = "~{prefix}.~{contig}",
            shard_bin_size_exact_match = shard_bin_size_exact_match,
            source_tag_truth_snv_indel_vcf = source_tag_truth_snv_indel_vcf,
            rename_id_string_vcf = rename_id_string_vcf,
            rename_id_string_truth_snv_indel_vcf = rename_id_string_truth_snv_indel_vcf,
            rename_id_strip_chr_vcf = rename_id_strip_chr_vcf,
            rename_id_strip_chr_truth_snv_indel_vcf = rename_id_strip_chr_truth_snv_indel_vcf,
            utils_docker = utils_docker,
            runtime_attr_rename_vcf = runtime_attr_rename_vcf,
            runtime_attr_rename_truth = runtime_attr_rename_truth,
            runtime_attr_create_exact_shards = runtime_attr_create_exact_shards,
            runtime_attr_subset_exact_vcf = runtime_attr_subset_exact_vcf,
            runtime_attr_subset_exact_truth = runtime_attr_subset_exact_truth,
            runtime_attr_exact_match = runtime_attr_exact_match,
            runtime_attr_append_exact_annotations = runtime_attr_append_exact_annotations,
            runtime_attr_concat_exact_annotations = runtime_attr_concat_exact_annotations,
            runtime_attr_concat_exact_unmatched = runtime_attr_concat_exact_unmatched
    }

    # Admit only records long enough for Truvari, measuring each VCF by its own length field
    call Helpers.SubsetVcfByLength as SubsetTruvariEval {
        input:
            vcf = ExactMatch.unmatched_vcf,
            vcf_idx = ExactMatch.unmatched_vcf_idx,
            length_field = length_field_vcf,
            min_length = min_sv_length_truvari_vcf,
            prefix = "~{prefix}.~{contig}.truvari_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_truvari_subset_vcf
    }

    call Helpers.SubsetVcfByLength as SubsetTruvariTruth {
        input:
            vcf = ExactMatch.truth_vcf,
            vcf_idx = ExactMatch.truth_vcf_idx,
            length_field = length_field_truth_snv_indel_vcf,
            min_length = min_sv_length_truvari_truth_snv_indel_vcf,
            prefix = "~{prefix}.~{contig}.truvari_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_truvari_subset_truth
    }

    call TruvariMatch.TruvariMatch {
        input:
            vcf = SubsetTruvariEval.subset_vcf,
            vcf_idx = SubsetTruvariEval.subset_vcf_idx,
            truth_snv_indel_vcf = SubsetTruvariTruth.subset_vcf,
            truth_snv_indel_vcf_idx = SubsetTruvariTruth.subset_vcf_idx,
            contig = contig,
            prefix = "~{prefix}.~{contig}.truvari",
            source_tag = source_tag_truth_snv_indel_vcf,
            shard_bin_size_truvari_match = shard_bin_size_truvari_match,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            utils_docker = utils_docker,
            runtime_attr_create_truvari_shards = runtime_attr_truvari_create_shards,
            runtime_attr_subset_truvari_vcf = runtime_attr_truvari_subset_region_vcf,
            runtime_attr_subset_truvari_truth = runtime_attr_truvari_subset_region_truth,
            runtime_attr_run_truvari_09 = runtime_attr_truvari_run_truvari_09,
            runtime_attr_run_truvari_07 = runtime_attr_truvari_run_truvari_07,
            runtime_attr_run_truvari_05 = runtime_attr_truvari_run_truvari_05,
            runtime_attr_concat_matched = runtime_attr_truvari_concat_matched,
            runtime_attr_concat_matched_truth = runtime_attr_truvari_concat_matched_truth,
            runtime_attr_concat_unmatched = runtime_attr_truvari_concat_unmatched
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

    # Admit only records long enough for the bedtools closest round; the SV truth was filtered above
    call Helpers.SubsetVcfByLength as SubsetBedtoolsEval {
        input:
            vcf = TruvariMatch.unmatched_vcf,
            vcf_idx = TruvariMatch.unmatched_vcf_idx,
            length_field = length_field_vcf,
            min_length = min_sv_length_bedtools_closest_vcf,
            prefix = "~{prefix}.~{contig}.bedtools_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_bedtools_subset_vcf
    }

    call BedtoolsClosestSV.BedtoolsClosestSV {
        input:
            vcf = SubsetBedtoolsEval.subset_vcf,
            vcf_idx = SubsetBedtoolsEval.subset_vcf_idx,
            truth_sv_vcf = truth_sv_vcf_final,
            truth_sv_vcf_idx = truth_sv_vcf_final_idx,
            prefix = "~{prefix}.~{contig}.bedtools_closest",
            type_field = type_field_vcf,
            length_field = length_field_vcf,
            move_dup_to_origin = move_dup_to_origin,
            source_tag = source_tag_truth_sv_vcf,
            gatk_sv_lr_docker = gatk_sv_lr_docker,
            utils_docker = utils_docker,
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

    output {
        File annotations_tsv_benchmark = BuildBenchmarkAnnotationTsv.merged_tsv
        File annotations_header_benchmark = BuildBenchmarkAnnotationTsv.merged_header
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
            # Resolve every master column to its index in this file once rather than once per row
            col_indices = [col_map.get(col) for col in master_header]
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                row = ['.' if i is None or i >= len(parts) else parts[i] for i in col_indices]
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
