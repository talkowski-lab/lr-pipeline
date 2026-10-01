version 1.0

import "../utils/BedtoolsClosestSV.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"
import "../utils/TruvariMatch.wdl"

workflow AnnotateCallsetOverlap {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
        File vcf_sv_truth
        File vcf_sv_truth_idx
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix

        Int? records_per_shard

        Boolean compare_annotations = true
        Boolean do_exact = true
        Boolean do_truvari = true
        Boolean do_bedtools_closest = true
        Int min_sv_length_eval_truvari
        Int min_sv_length_truth_truvari
        Int min_sv_length_eval_bedtools_closest
        Int min_sv_length_truth_bedtools_closest
        String type_field_eval = "SVTYPE"
        String length_field_eval = "SVLEN"
        String skip_vep_categories = ""
        String? args_string_vcf
        String? args_string_vcf_truth
        String? args_string_vcf_sv_truth
        String? rename_id_string_vcf
        String? rename_id_string_vcf_truth
        String? rename_id_string_vcf_sv_truth
        Boolean? rename_id_strip_chr_vcf
        Boolean? rename_id_strip_chr_vcf_truth
        Boolean? rename_id_strip_chr_vcf_sv_truth

        String benchmark_annotations_docker
        String? sv_pipeline_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_eval
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_subset_sv_truth
        RuntimeAttr? runtime_attr_strip_genotypes
        RuntimeAttr? runtime_attr_rename_eval
        RuntimeAttr? runtime_attr_rename_truth
        RuntimeAttr? runtime_attr_rename_sv_truth
        RuntimeAttr? runtime_attr_exact_match
        RuntimeAttr? runtime_attr_truvari_subset_eval
        RuntimeAttr? runtime_attr_truvari_subset_truth
        RuntimeAttr? runtime_attr_truvari_run_truvari
        RuntimeAttr? runtime_attr_bedtools_subset_eval
        RuntimeAttr? runtime_attr_bedtools_subset_truth
        RuntimeAttr? runtime_attr_bedtools_convert_to_symbolic
        RuntimeAttr? runtime_attr_bedtools_match
        RuntimeAttr? runtime_attr_build_annotation_tsv
        RuntimeAttr? runtime_attr_collect_matched_ids
        RuntimeAttr? runtime_attr_extract_eval_vep_header
        RuntimeAttr? runtime_attr_extract_truth_vep_header
        RuntimeAttr? runtime_attr_shard_matched_eval
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
                    vcf = vcf_eval,
                    vcf_idx = vcf_eval_idx,
                    include_args = args_string_vcf,
                    extra_args = if single_contig then "" else "--regions ~{contig}",
                    prefix = "~{prefix}.~{contig}.eval",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_eval
            }
        }

        if (!single_contig || defined(args_string_vcf_truth)) {
            call Helpers.SubsetVcfByArgs as SubsetTruth {
                input:
                    vcf = vcf_truth,
                    vcf_idx = vcf_truth_idx,
                    include_args = args_string_vcf_truth,
                    extra_args = if single_contig then "" else "--regions ~{contig}",
                    prefix = "~{prefix}.~{contig}.truth",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_truth
            }
        }

        if (!single_contig || defined(args_string_vcf_sv_truth)) {
            call Helpers.SubsetVcfByArgs as SubsetSVTruth {
                input:
                    vcf = vcf_sv_truth,
                    vcf_idx = vcf_sv_truth_idx,
                    include_args = args_string_vcf_sv_truth,
                    extra_args = if single_contig then "" else "--regions ~{contig}",
                    prefix = "~{prefix}.~{contig}.sv_truth",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_sv_truth
            }
        }

        File subset_eval_vcf = select_first([SubsetEval.subset_vcf, vcf_eval])
        File subset_eval_vcf_idx = select_first([SubsetEval.subset_vcf_idx, vcf_eval_idx])
        File subset_truth_vcf = select_first([SubsetTruth.subset_vcf, vcf_truth])
        File subset_truth_vcf_idx = select_first([SubsetTruth.subset_vcf_idx, vcf_truth_idx])
        File subset_sv_truth_vcf = select_first([SubsetSVTruth.subset_vcf, vcf_sv_truth])
        File subset_sv_truth_vcf_idx = select_first([SubsetSVTruth.subset_vcf_idx, vcf_sv_truth_idx])

        call Helpers.StripGenotypes {
            input:
                vcf = subset_eval_vcf,
                vcf_idx = subset_eval_vcf_idx,
                prefix = "~{prefix}.~{contig}.eval",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_strip_genotypes
        }

        if (defined(rename_id_string_vcf)) {
            call Helpers.RenameVariantIds as RenameEvalIds {
                input:
                    vcf = StripGenotypes.stripped_vcf,
                    vcf_idx = StripGenotypes.stripped_vcf_idx,
                    prefix = "~{prefix}.~{contig}.eval.renamed",
                    id_format = select_first([rename_id_string_vcf]),
                    strip_chr = select_first([rename_id_strip_chr_vcf, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_eval
            }
        }

        if (defined(rename_id_string_vcf_truth)) {
            call Helpers.RenameVariantIds as RenameTruthIds {
                input:
                    vcf = subset_truth_vcf,
                    vcf_idx = subset_truth_vcf_idx,
                    prefix = "~{prefix}.~{contig}.truth.renamed",
                    id_format = select_first([rename_id_string_vcf_truth]),
                    strip_chr = select_first([rename_id_strip_chr_vcf_truth, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_truth
            }
        }

        if (defined(rename_id_string_vcf_sv_truth)) {
            call Helpers.RenameVariantIds as RenameSVTruthIds {
                input:
                    vcf = subset_sv_truth_vcf,
                    vcf_idx = subset_sv_truth_vcf_idx,
                    prefix = "~{prefix}.~{contig}.sv_truth.renamed",
                    id_format = select_first([rename_id_string_vcf_sv_truth]),
                    strip_chr = select_first([rename_id_strip_chr_vcf_sv_truth, false]),
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_sv_truth
            }
        }

        File eval_vcf_final = select_first([RenameEvalIds.renamed_vcf, StripGenotypes.stripped_vcf])
        File eval_vcf_final_idx = select_first([RenameEvalIds.renamed_vcf_idx, StripGenotypes.stripped_vcf_idx])
        File truth_vcf_final = select_first([RenameTruthIds.renamed_vcf, subset_truth_vcf])
        File truth_vcf_final_idx = select_first([RenameTruthIds.renamed_vcf_idx, subset_truth_vcf_idx])
        File sv_truth_vcf_final = select_first([RenameSVTruthIds.renamed_vcf, subset_sv_truth_vcf])
        File sv_truth_vcf_final_idx = select_first([RenameSVTruthIds.renamed_vcf_idx, subset_sv_truth_vcf_idx])

        if (do_exact) {
            call ExactMatch {
                input:
                    vcf_eval = eval_vcf_final,
                    vcf_eval_idx = eval_vcf_final_idx,
                    vcf_truth = truth_vcf_final,
                    vcf_truth_idx = truth_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.exact",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_exact_match
            }
        }

        File post_exact_vcf = select_first([ExactMatch.unmatched_vcf, eval_vcf_final])
        File post_exact_vcf_idx = select_first([ExactMatch.unmatched_vcf_idx, eval_vcf_final_idx])

        if (do_truvari) {
            call TruvariMatch.TruvariMatch {
                input:
                    vcf_eval = post_exact_vcf,
                    vcf_eval_idx = post_exact_vcf_idx,
                    vcf_truth = truth_vcf_final,
                    vcf_truth_idx = truth_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.truvari",
                    min_sv_length_eval = min_sv_length_eval_truvari,
                    min_sv_length_truth = min_sv_length_truth_truvari,
                    length_field_eval = length_field_eval,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    utils_docker = utils_docker,
                    runtime_attr_subset_eval = runtime_attr_truvari_subset_eval,
                    runtime_attr_subset_truth = runtime_attr_truvari_subset_truth,
                    runtime_attr_run_truvari = runtime_attr_truvari_run_truvari
            }
        }

        File post_truvari_vcf = select_first([TruvariMatch.unmatched_vcf, post_exact_vcf])
        File post_truvari_vcf_idx = select_first([TruvariMatch.unmatched_vcf_idx, post_exact_vcf_idx])

        if (do_bedtools_closest) {
            call BedtoolsClosestSV.BedtoolsClosestSV {
                input:
                    vcf_eval = post_truvari_vcf,
                    vcf_eval_idx = post_truvari_vcf_idx,
                    vcf_sv_truth = sv_truth_vcf_final,
                    vcf_sv_truth_idx = sv_truth_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.bedtools_closest",
                    min_sv_length_eval = min_sv_length_eval_bedtools_closest,
                    min_sv_length_truth = min_sv_length_truth_bedtools_closest,
                    type_field_eval = type_field_eval,
                    length_field_eval = length_field_eval,
                    sv_pipeline_docker = select_first([sv_pipeline_docker, benchmark_annotations_docker]),
                    utils_docker = utils_docker,
                    runtime_attr_subset_eval = runtime_attr_bedtools_subset_eval,
                    runtime_attr_subset_truth = runtime_attr_bedtools_subset_truth,
                    runtime_attr_convert_to_symbolic = runtime_attr_bedtools_convert_to_symbolic,
                    runtime_attr_match = runtime_attr_bedtools_match
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

        if (do_annotation_summaries) {
            call ExtractVepHeader as ExtractTruthVepHeader {
                input:
                    vcf = truth_vcf_final,
                    vcf_idx = truth_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.truth",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_truth_vep_header
            }

            call ExtractVepHeader as ExtractEvalVepHeader {
                input:
                    vcf = eval_vcf_final,
                    vcf_idx = eval_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.eval",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_extract_eval_vep_header
            }

            call CollectMatchedIDsAndINFO {
                input:
                    annotation_tsv = BuildAnnotationTsv.concatenated_tsv,
                    vcf_eval = eval_vcf_final,
                    vcf_eval_idx = eval_vcf_final_idx,
                    vcf_truth_snv = truth_vcf_final,
                    vcf_truth_snv_idx = truth_vcf_final_idx,
                    vcf_truth_sv = sv_truth_vcf_final,
                    vcf_truth_sv_idx = sv_truth_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.collected",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_collect_matched_ids
            }

            call ComputeSummaryForContig {
                input:
                    eval_vcf = eval_vcf_final,
                    eval_vcf_idx = eval_vcf_final_idx,
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
                        runtime_attr_override = runtime_attr_shard_matched_eval
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
                tsvs = BuildAnnotationTsv.concatenated_tsv,
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
        File annotations_tsv_benchmark = select_first([MergeAnnotationTsvs.concatenated_tsv, BuildAnnotationTsv.concatenated_tsv[0]])
        File? benchmark_annotations_summary_tsv = benchmark_annotations_summary_tsv_local
        File? benchmark_annotations_stats_tsv = benchmark_annotations_stats_tsv_local
        File? benchmark_annotations_plots_tarball = benchmark_annotations_plots_tarball_local
    }
}

task ExactMatch {
    input {
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth
        File vcf_truth_idx
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
            ~{vcf_eval} \
            ~{vcf_truth}

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' \
            isec_matched/0000.vcf \
            > eval_matched.tsv

        bcftools query \
            -f '%ID\t%FILTER\t%INFO/AF\n' \
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
                print $1, out, $3
            }' > truth_matched.tsv

        paste eval_matched.tsv truth_matched.tsv \
            | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,"EXACT",$6,"SNV_indel",$7,$8}' \
            > ~{prefix}.tsv

        bcftools isec \
            -C \
            -c none \
            -p isec_unmatched \
            ~{vcf_eval} \
            ~{vcf_truth}

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
        disk_gb: 5 * ceil(size(vcf_eval, "GB") + size(vcf_truth, "GB")) + 5,
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
        File vcf_eval
        File vcf_eval_idx
        File vcf_truth_snv
        File vcf_truth_snv_idx
        File vcf_truth_sv
        File vcf_truth_sv_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import subprocess

annotation_tsv = "~{annotation_tsv}"
vcf_eval = "~{vcf_eval}"
vcf_truth_snv = "~{vcf_truth_snv}"
vcf_truth_sv = "~{vcf_truth_sv}"
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
for vcf in [vcf_truth_snv, vcf_truth_sv]:
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
        disk_gb: 2 * ceil(size(vcf_eval, "GB") + size(vcf_truth_snv, "GB") + size(vcf_truth_sv, "GB")) + 10,
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
        File eval_vcf
        File eval_vcf_idx
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
            --eval_vcf ~{eval_vcf} \
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
        disk_gb: 50 * ceil(size(eval_vcf, "GB")) + 5,
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
