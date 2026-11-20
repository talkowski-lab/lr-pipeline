version 1.0

import "../utils/BedtoolsClosestSV.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "../utils/TruvariMatch.wdl"

workflow AnnotateCallsetOverlap_AF {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Boolean normalize_vcf = false
        Boolean create_variant_attributes = false

        Boolean compare_annotations = true
        Boolean do_exact = true
        Boolean do_truvari = true
        Boolean do_bedtools_closest = true

        Int min_sv_length_truvari
        Int min_sv_length_truth_truvari
        Int min_sv_length_bedtools_closest
        Int min_sv_length_truth_bedtools_closest
        String type_field = "allele_type"
        String length_field = "allele_length"
        String source_tag_truth_snv_indel_vcf = "SNV_indel"
        String source_tag_truth_sv_vcf = "SV"
        String normalize_check_ref = "w"
        String skip_vep_categories = ""
        String af_field_sv_truth = "AF"
        String ac_field_sv_truth = "AC"
        String an_field_sv_truth = "AN"

        String? args_string_vcf
        String? args_string_truth_snv_indel_vcf
        String? args_string_truth_sv_vcf
        String? rename_id_string_vcf
        String? rename_id_string_truth_snv_indel_vcf
        String? rename_id_string_truth_sv_vcf
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_truth_snv_indel_vcf
        Boolean? rename_id_strip_chr_truth_sv_vcf

        String benchmark_annotations_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_subset_sv_truth
        RuntimeAttr? runtime_attr_strip_genotypes
        RuntimeAttr? runtime_attr_normalize_vcf
        RuntimeAttr? runtime_attr_annotate_attributes_vcf
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_rename_sv_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_truvari_run_truvari_09
        RuntimeAttr? runtime_attr_truvari_run_truvari_07
        RuntimeAttr? runtime_attr_truvari_run_truvari_05
        RuntimeAttr? runtime_attr_truvari_concat_matched
        RuntimeAttr? runtime_attr_bedtools_subset_vcf
        RuntimeAttr? runtime_attr_bedtools_subset_truth
        RuntimeAttr? runtime_attr_bedtools_convert_to_symbolic
        RuntimeAttr? runtime_attr_bedtools_split_vcf
        RuntimeAttr? runtime_attr_bedtools_split_truth
        RuntimeAttr? runtime_attr_bedtools_compare
        RuntimeAttr? runtime_attr_bedtools_calculate
        RuntimeAttr? runtime_attr_bedtools_merge_comparisons
        RuntimeAttr? runtime_attr_build_annotation_tsv
        RuntimeAttr? runtime_attr_append_truth_af
        RuntimeAttr? runtime_attr_collect_matched_ids
        RuntimeAttr? runtime_attr_extract_vcf_vep_header
        RuntimeAttr? runtime_attr_extract_truth_vep_header
        RuntimeAttr? runtime_attr_shard_matched_vcf
        RuntimeAttr? runtime_attr_compute_shard_benchmarks
        RuntimeAttr? runtime_attr_merge_shard_benchmarks
        RuntimeAttr? runtime_attr_compute_summary_for_contig
        RuntimeAttr? runtime_attr_merge_annotation_tsvs
        RuntimeAttr? runtime_attr_merge_benchmark_summaries
        RuntimeAttr? runtime_attr_merge_summary_stats
        RuntimeAttr? runtime_attr_merge_plot_tarballs
    }

    Boolean single_contig = length(contigs) == 1
    Boolean any_comparison_enabled = do_exact || do_truvari || do_bedtools_closest
    Boolean do_annotation_summaries = compare_annotations && any_comparison_enabled

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

        call Helpers.StripGenotypes {
            input:
                vcf = vcf_subsetted,
                vcf_idx = vcf_subsetted_idx,
                prefix = "~{prefix}.~{contig}.eval",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_strip_genotypes
        }

        if (normalize_vcf) {
            call Helpers.NormalizeVcf as NormalizeEval {
                input:
                    vcf = StripGenotypes.stripped_vcf,
                    vcf_idx = StripGenotypes.stripped_vcf_idx,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    check_ref = normalize_check_ref,
                    prefix = "~{prefix}.~{contig}.eval.normalized",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_normalize_vcf
            }
        }

        File vcf_normalized = select_first([NormalizeEval.normalized_vcf, StripGenotypes.stripped_vcf])
        File vcf_normalized_idx = select_first([NormalizeEval.normalized_vcf_idx, StripGenotypes.stripped_vcf_idx])

        if (create_variant_attributes) {
            call Helpers.AnnotateVariantAttributes as AnnotateEvalAttributes {
                input:
                    vcf = vcf_normalized,
                    vcf_idx = vcf_normalized_idx,
                    prefix = "~{prefix}.~{contig}.eval.attributes",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_annotate_attributes_vcf
            }
        }

        File vcf_attributed = select_first([AnnotateEvalAttributes.annotated_vcf, vcf_normalized])
        File vcf_attributed_idx = select_first([AnnotateEvalAttributes.annotated_vcf_idx, vcf_normalized_idx])

        if (defined(rename_id_string_vcf)) {
            call Helpers.RenameVariantIds as RenameEvalIds {
                input:
                    vcf = vcf_attributed,
                    vcf_idx = vcf_attributed_idx,
                    prefix = "~{prefix}.~{contig}.eval.renamed",
                    id_format = select_first([rename_id_string_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_vcf
            }
        }

        if (defined(rename_id_string_truth_snv_indel_vcf)) {
            call Helpers.RenameVariantIds as RenameTruthIds {
                input:
                    vcf = truth_snv_indel_vcf_subsetted,
                    vcf_idx = truth_snv_indel_vcf_subsetted_idx,
                    prefix = "~{prefix}.~{contig}.truth.renamed",
                    id_format = select_first([rename_id_string_truth_snv_indel_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_truth_snv_indel_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_truth
            }
        }

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

        File vcf_final = select_first([RenameEvalIds.renamed_vcf, vcf_attributed])
        File vcf_final_idx = select_first([RenameEvalIds.renamed_vcf_idx, vcf_attributed_idx])
        File truth_snv_indel_vcf_final = select_first([RenameTruthIds.renamed_vcf, truth_snv_indel_vcf_subsetted])
        File truth_snv_indel_vcf_final_idx = select_first([RenameTruthIds.renamed_vcf_idx, truth_snv_indel_vcf_subsetted_idx])
        File truth_sv_vcf_final = select_first([RenameSVTruthIds.renamed_vcf, truth_sv_vcf_subsetted])
        File truth_sv_vcf_final_idx = select_first([RenameSVTruthIds.renamed_vcf_idx, truth_sv_vcf_subsetted_idx])

        if (do_exact) {
            call ExactMatch {
                input:
                    vcf = vcf_final,
                    vcf_idx = vcf_final_idx,
                    truth_snv_indel_vcf = truth_snv_indel_vcf_final,
                    truth_snv_indel_vcf_idx = truth_snv_indel_vcf_final_idx,
                    source_tag = source_tag_truth_snv_indel_vcf,
                    prefix = "~{prefix}.~{contig}.exact",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_exact_match
            }
        }

        File vcf_post_exact = select_first([ExactMatch.unmatched_vcf, vcf_final])
        File vcf_post_exact_idx = select_first([ExactMatch.unmatched_vcf_idx, vcf_final_idx])

        if (do_truvari) {
            call TruvariMatch.TruvariMatch {
                input:
                    vcf = vcf_post_exact,
                    vcf_idx = vcf_post_exact_idx,
                    truth_snv_indel_vcf = truth_snv_indel_vcf_final,
                    truth_snv_indel_vcf_idx = truth_snv_indel_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.truvari",
                    min_sv_length = min_sv_length_truvari,
                    min_sv_length_truth = min_sv_length_truth_truvari,
                    length_field = length_field,
                    source_tag = source_tag_truth_snv_indel_vcf,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    utils_docker = utils_docker,
                    runtime_attr_subset_vcf = runtime_attr_truvari_subset_vcf,
                    runtime_attr_subset_truth = runtime_attr_truvari_subset_truth,
                    runtime_attr_run_truvari_09 = runtime_attr_truvari_run_truvari_09,
                    runtime_attr_run_truvari_07 = runtime_attr_truvari_run_truvari_07,
                    runtime_attr_run_truvari_05 = runtime_attr_truvari_run_truvari_05,
                    runtime_attr_concat_matched = runtime_attr_truvari_concat_matched
            }
        }

        File vcf_post_truvari = select_first([TruvariMatch.unmatched_vcf, vcf_post_exact])
        File vcf_post_truvari_idx = select_first([TruvariMatch.unmatched_vcf_idx, vcf_post_exact_idx])

        if (do_bedtools_closest) {
            call BedtoolsClosestSV.BedtoolsClosestSV {
                input:
                    vcf = vcf_post_truvari,
                    vcf_idx = vcf_post_truvari_idx,
                    truth_sv_vcf = truth_sv_vcf_final,
                    truth_sv_vcf_idx = truth_sv_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.bedtools_closest",
                    min_sv_length = min_sv_length_bedtools_closest,
                    min_sv_length_truth = min_sv_length_truth_bedtools_closest,
                    type_field = type_field,
                    length_field = length_field,
                    source_tag = source_tag_truth_sv_vcf,
                    benchmark_annotations_docker = benchmark_annotations_docker,
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
        }

        Array[File] annotation_tsvs_to_merge = select_all([
            ExactMatch.annotation_tsv,
            TruvariMatch.annotation_tsv,
            BedtoolsClosestSV.annotation_tsv,
        ])

        call Helpers.ConcatTsvs as BuildAnnotationTsv {
            input:
                tsvs = annotation_tsvs_to_merge,
                sort_output = true,
                prefix = "~{prefix}.~{contig}.annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_build_annotation_tsv
        }

        call AppendTruthAFACAN {
            input:
                annotation_tsv = BuildAnnotationTsv.concatenated_tsv,
                truth_snv_indel_vcf = truth_snv_indel_vcf_final,
                truth_snv_indel_vcf_idx = truth_snv_indel_vcf_final_idx,
                truth_sv_vcf = truth_sv_vcf_final,
                truth_sv_vcf_idx = truth_sv_vcf_final_idx,
                af_field_sv_truth = af_field_sv_truth,
                ac_field_sv_truth = ac_field_sv_truth,
                an_field_sv_truth = an_field_sv_truth,
                prefix = "~{prefix}.~{contig}.annotations_af",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_append_truth_af
        }

        if (do_annotation_summaries) {
            call ExtractVepHeader as ExtractTruthVepHeader {
                input:
                    vcf = truth_snv_indel_vcf_final,
                    vcf_idx = truth_snv_indel_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.truth",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_truth_vep_header
            }

            call ExtractVepHeader as ExtractEvalVepHeader {
                input:
                    vcf = vcf_final,
                    vcf_idx = vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.eval",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_vcf_vep_header
            }

            call CollectMatchedIDsAndINFO {
                input:
                    annotation_tsv = BuildAnnotationTsv.concatenated_tsv,
                    vcf = vcf_final,
                    vcf_idx = vcf_final_idx,
                    truth_snv_indel_vcf = truth_snv_indel_vcf_final,
                    truth_snv_indel_vcf_idx = truth_snv_indel_vcf_final_idx,
                    truth_sv_vcf = truth_sv_vcf_final,
                    truth_sv_vcf_idx = truth_sv_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.collected",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_collect_matched_ids
            }

            call ComputeSummaryForContig {
                input:
                    vcf = vcf_final,
                    vcf_idx = vcf_final_idx,
                    annotation_tsv = BuildAnnotationTsv.concatenated_tsv,
                    matched_with_info_tsv = CollectMatchedIDsAndINFO.matched_with_info_tsv,
                    eval_vep_header = ExtractEvalVepHeader.vep_header_txt,
                    truth_vep_header = ExtractTruthVepHeader.vep_header_txt,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.summary",
                    docker = benchmark_annotations_docker,
                    runtime_attr_override = runtime_attr_compute_summary_for_contig
            }

            if (defined(records_per_shard)) {
                call ShardedMatchedVariants {
                    input:
                        matched_with_info_tsv = CollectMatchedIDsAndINFO.matched_with_info_tsv,
                        records_per_shard = select_first([records_per_shard]),
                        prefix = "~{prefix}.~{contig}.sharded",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_shard_matched_vcf
                }
            }

            Array[File] matched_tsvs_to_process = select_first([ShardedMatchedVariants.shard_tsvs, [CollectMatchedIDsAndINFO.matched_with_info_tsv]])

            scatter (i in range(length(matched_tsvs_to_process))) {
                call ComputeShardBenchmarks {
                    input:
                        matched_shard_tsv = matched_tsvs_to_process[i],
                        eval_vep_header = ExtractEvalVepHeader.vep_header_txt,
                        truth_vep_header = ExtractTruthVepHeader.vep_header_txt,
                        skip_vep_categories = skip_vep_categories,
                        contig = contig,
                        shard_label = "~{i}",
                        prefix = "~{prefix}.~{contig}.shard_~{i}",
                        docker = benchmark_annotations_docker,
                        runtime_attr_override = runtime_attr_compute_shard_benchmarks
                }
            }

            call MergeShardBenchmarks {
                input:
                    af_pair_tsvs = ComputeShardBenchmarks.af_pairs_tsv,
                    vep_pair_tsvs = ComputeShardBenchmarks.vep_pairs_tsv,
                    skip_vep_categories = skip_vep_categories,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.merged",
                    docker = benchmark_annotations_docker,
                    runtime_attr_override = runtime_attr_merge_shard_benchmarks
            }
        }
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs as MergeAnnotationTsvs {
            input:
                tsvs = AppendTruthAFACAN.annotated_tsv,
                sort_output = false,
                prefix = "~{prefix}.benchmark_annotations",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_annotation_tsvs
        }
    }

    if (!single_contig && do_annotation_summaries) {
        call Helpers.ConcatTsvs as MergeBenchmarkSummaries {
            input:
                tsvs = select_all(ComputeSummaryForContig.benchmark_summary_tsv),
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.benchmark_summary",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_benchmark_summaries
        }

        call Helpers.ConcatTsvs as MergeSummaryStats {
            input:
                tsvs = select_all(ComputeSummaryForContig.summary_stats_tsv),
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.summary_stats",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_summary_stats
        }

        call MergePlotTarballs {
            input:
                tarballs = select_all(MergeShardBenchmarks.plot_tarball),
                prefix = "~{prefix}.plot_tarballs",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_merge_plot_tarballs
        }
    }

    if (do_annotation_summaries) {
        File? benchmark_annotations_summary_tsv_local = if single_contig then ComputeSummaryForContig.benchmark_summary_tsv[0] else MergeBenchmarkSummaries.concatenated_tsv
        File? benchmark_annotations_stats_tsv_local = if single_contig then ComputeSummaryForContig.summary_stats_tsv[0] else MergeSummaryStats.concatenated_tsv
        File? benchmark_annotations_plots_tarball_local = if single_contig then MergeShardBenchmarks.plot_tarball[0] else MergePlotTarballs.merged_tarball
    }

    output {
        File annotations_tsv_benchmark = select_first([MergeAnnotationTsvs.concatenated_tsv, AppendTruthAFACAN.annotated_tsv[0]])
        File? benchmark_annotations_summary_tsv = benchmark_annotations_summary_tsv_local
        File? benchmark_annotations_stats_tsv = benchmark_annotations_stats_tsv_local
        File? benchmark_annotations_plots_tarball = benchmark_annotations_plots_tarball_local
    }
}

task ExactMatch {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String source_tag
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools isec \
            -c none \
            -n=2 \
            -p isec_matched \
            ~{vcf} \
            ~{truth_snv_indel_vcf}

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            isec_matched/0000.vcf \
            > eval_matched.tsv

        bcftools query \
            -f '%ID\t%FILTER\n' \
            isec_matched/0001.vcf \
            | awk -F'\t' 'BEGIN{OFS="\t"} {
                n = split($2, parts, ";")
                out = ""
                for (i = 1; i <= n; i++) {
                    if (parts[i] != "." && parts[i] != "PASS") {
                        out = (out == "" ? parts[i] : out "," parts[i])
                    }
                }
                if (out == "") out = "."
                print $1, out
            }' > truth_matched.tsv

        paste eval_matched.tsv truth_matched.tsv \
            | awk -v src="~{source_tag}" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,"EXACT",$6,src,$7}' \
            > ~{prefix}.tsv

        bcftools isec \
            -C \
            -c none \
            -p isec_unmatched \
            ~{vcf} \
            ~{truth_snv_indel_vcf}

        bgzip -c isec_unmatched/0000.vcf > ~{prefix}.vcf.gz

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotation_tsv = "~{prefix}.tsv"
        File unmatched_vcf = "~{prefix}.vcf.gz"
        File unmatched_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(vcf, "GB") + size(truth_snv_indel_vcf, "GB")) + 5,
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

task CollectMatchedIDsAndINFO {
    input {
        File annotation_tsv
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import subprocess

annotation_tsv = "~{annotation_tsv}"
vcf_eval = "~{vcf}"
truth_snv_indel_vcf = "~{truth_snv_indel_vcf}"
truth_sv_vcf = "~{truth_sv_vcf}"
prefix = "~{prefix}"

eval_to_truth = {}
eval_ids = set()
truth_ids = set()

with open(annotation_tsv) as f:
    for line in f:
        fields = line.strip().split('\t')
        eval_id = fields[4]
        truth_id = fields[6]
        eval_to_truth[eval_id] = truth_id
        eval_ids.add(eval_id)
        truth_ids.add(truth_id)

eval_info = {}
cmd = f"bcftools query -f '%ID\\t%INFO\\n' {vcf_eval}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    parts = line.strip().split('\t', 1)
    if len(parts) == 2 and parts[0] in eval_ids:
        eval_info[parts[0]] = parts[1]
proc.wait()

truth_info = {}
for vcf in [truth_snv_indel_vcf, truth_sv_vcf]:
    cmd = f"bcftools query -f '%ID\\t%INFO\\n' {vcf}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        parts = line.strip().split('\t', 1)
        if len(parts) == 2 and parts[0] in truth_ids:
            truth_info[parts[0]] = parts[1]
    proc.wait()

with open(f"{prefix}.matched_with_info.tsv", 'w') as out:
    for eval_id, truth_id in eval_to_truth.items():
        eval_inf = eval_info.get(eval_id, '.')
        truth_inf = truth_info.get(truth_id, '.')
        out.write(f"{eval_id}\t{truth_id}\t{eval_inf}\t{truth_inf}\n")

EOF
    >>>

    output {
        File matched_with_info_tsv = "~{prefix}.matched_with_info.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(truth_snv_indel_vcf, "GB") + size(truth_sv_vcf, "GB")) + 10,
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

task AppendTruthAFACAN {
    input {
        File annotation_tsv
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        String af_field_sv_truth = "AF"
        String ac_field_sv_truth = "AC"
        String an_field_sv_truth = "AN"
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import subprocess

annotation_tsv = "~{annotation_tsv}"
truth_snv_indel_vcf = "~{truth_snv_indel_vcf}"
truth_sv_vcf = "~{truth_sv_vcf}"
af_field_sv_truth = "~{af_field_sv_truth}"
ac_field_sv_truth = "~{ac_field_sv_truth}"
an_field_sv_truth = "~{an_field_sv_truth}"
prefix = "~{prefix}"

truth_af = {}
for vcf, af_field, ac_field, an_field in [
    (truth_snv_indel_vcf, "AF", "AC", "AN"),
    (truth_sv_vcf, af_field_sv_truth, ac_field_sv_truth, an_field_sv_truth)
]:
    cmd = f"bcftools query -f '%ID\\t%INFO/{af_field}\\t%INFO/{ac_field}\\t%INFO/{an_field}\\n' {vcf}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        parts = line.rstrip('\n').split('\t')
        if len(parts) == 4:
            truth_af[parts[0]] = (parts[1], parts[2], parts[3])
    proc.wait()

with open(annotation_tsv) as f, open(f"{prefix}.tsv", 'w') as out:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        truth_id = fields[6]
        af, ac, an = truth_af.get(truth_id, ('.', '.', '.'))
        out.write('\t'.join(fields + [af, ac, an]) + '\n')

EOF
    >>>

    output {
        File annotated_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(annotation_tsv, "GB") + size(truth_snv_indel_vcf, "GB") + size(truth_sv_vcf, "GB")) + 10,
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

task ExtractVepHeader {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools view \
            -h \
            ~{vcf} \
        | awk 'BEGIN{IGNORECASE=1} /^##INFO=<ID=(vep|csq),/ {print}' > ~{prefix}_vep_header.txt
    >>>

    output {
        File vep_header_txt = "~{prefix}_vep_header.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task ShardedMatchedVariants {
    input {
        File matched_with_info_tsv
        Int records_per_shard
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p shards

        cat ~{matched_with_info_tsv} \
            | awk 'BEGIN{c=0;f=0} {print > sprintf("shards/matched.%06d.tsv", int(c/~{records_per_shard})) ; c++} END{ }'
    >>>

    output {
        Array[File] shard_tsvs = glob("shards/*.tsv")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(matched_with_info_tsv, "GB")) + 5,
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

task ComputeShardBenchmarks {
    input {
        File matched_shard_tsv
        File eval_vep_header
        File truth_vep_header
        String skip_vep_categories
        String contig
        String shard_label
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/compute_benchmarks_shard.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --matched_shard_tsv ~{matched_shard_tsv} \
            --eval_vep_header ~{eval_vep_header} \
            --truth_vep_header ~{truth_vep_header} \
            --shard_label ~{shard_label} \
            ~{if skip_vep_categories != "" then "--skip_vep_categories " + skip_vep_categories else ""}
    >>>

    output {
        File af_pairs_tsv = "~{prefix}.shard_~{shard_label}.af_pairs.tsv"
        File vep_pairs_tsv = "~{prefix}.shard_~{shard_label}.vep_pairs.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(matched_shard_tsv, "GB")) + 5,
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

task MergeShardBenchmarks {
    input {
        Array[File] af_pair_tsvs
        Array[File] vep_pair_tsvs
        String skip_vep_categories
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/merge_benchmarks_from_pairs.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --af_pair_tsvs ~{sep=',' af_pair_tsvs} \
            --vep_pair_tsvs ~{sep=',' vep_pair_tsvs} \
            ~{if skip_vep_categories != "" then "--skip_vep_categories " + skip_vep_categories else ""}
    >>>

    output {
        File plot_tarball = "~{prefix}.benchmarks.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 50 * ceil(size(af_pair_tsvs, "GB")) + 5,
        boot_disk_gb: 50,
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

task ComputeSummaryForContig {
    input {
        File vcf
        File vcf_idx
        File annotation_tsv
        File matched_with_info_tsv
        File eval_vep_header
        File truth_vep_header
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 /opt/gnomad-lr/scripts/benchmark/compute_summary_for_contig.py \
            --prefix ~{prefix} \
            --contig ~{contig} \
            --eval_vcf ~{vcf} \
            --annotation_tsv ~{annotation_tsv} \
            --matched_with_info_tsv ~{matched_with_info_tsv} \
            --eval_vep_header ~{eval_vep_header} \
            --truth_vep_header ~{truth_vep_header}
    >>>

    output {
        File benchmark_summary_tsv = "~{prefix}.benchmark_summary.tsv"
        File summary_stats_tsv = "~{prefix}.summary_stats.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 4,
        disk_gb: 50 * ceil(size(vcf, "GB")) + 5,
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

task MergePlotTarballs {
    input {
        Array[File] tarballs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p final_results/AF_plots
        mkdir -p final_results/VEP_plots

        for tarball in ~{sep=' ' tarballs}; do
            tar -xvf $tarball --strip-components=1 -C final_results
        done

        tar -czf ~{prefix}.plots.tar.gz final_results/
    >>>

    output {
        File merged_tarball = "~{prefix}.plots.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 * ceil(size(tarballs, "GB")) + 5,
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
