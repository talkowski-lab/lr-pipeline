version 1.0

import "../utils/BedtoolsClosestSV.wdl"
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

        Boolean do_exact = true
        Boolean do_truvari = true
        Boolean do_bedtools_closest = true

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

        String benchmark_annotations_docker
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
        RuntimeAttr? runtime_attr_concat_exact_unmatched
        RuntimeAttr? runtime_attr_append_exact_annotations
        RuntimeAttr? runtime_attr_truvari_subset_vcf
        RuntimeAttr? runtime_attr_truvari_subset_truth
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

    call Helpers.StripGenotypes {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = "~{prefix}.eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_strip_genotypes
    }

    scatter (contig in contigs) {
        if (!single_contig || defined(args_string_vcf)) {
            call Helpers.SubsetVcfByArgs as SubsetEval {
                input:
                    vcf = StripGenotypes.stripped_vcf,
                    vcf_idx = StripGenotypes.stripped_vcf_idx,
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

        File vcf_subsetted = select_first([SubsetEval.subset_vcf, StripGenotypes.stripped_vcf])
        File vcf_subsetted_idx = select_first([SubsetEval.subset_vcf_idx, StripGenotypes.stripped_vcf_idx])
        File truth_snv_indel_vcf_subsetted = select_first([SubsetTruth.subset_vcf, truth_snv_indel_vcf])
        File truth_snv_indel_vcf_subsetted_idx = select_first([SubsetTruth.subset_vcf_idx, truth_snv_indel_vcf_idx])
        File truth_sv_vcf_subsetted = select_first([SubsetSVTruth.subset_vcf, truth_sv_vcf])
        File truth_sv_vcf_subsetted_idx = select_first([SubsetSVTruth.subset_vcf_idx, truth_sv_vcf_idx])

        if (defined(rename_id_string_vcf)) {
            call Helpers.RenameVariantIds as RenameEvalIds {
                input:
                    vcf = vcf_subsetted,
                    vcf_idx = vcf_subsetted_idx,
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

        File vcf_final = select_first([RenameEvalIds.renamed_vcf, vcf_subsetted])
        File vcf_final_idx = select_first([RenameEvalIds.renamed_vcf_idx, vcf_subsetted_idx])
        File truth_snv_indel_vcf_final = select_first([RenameTruthIds.renamed_vcf, truth_snv_indel_vcf_subsetted])
        File truth_snv_indel_vcf_final_idx = select_first([RenameTruthIds.renamed_vcf_idx, truth_snv_indel_vcf_subsetted_idx])
        File truth_sv_vcf_final = select_first([RenameSVTruthIds.renamed_vcf, truth_sv_vcf_subsetted])
        File truth_sv_vcf_final_idx = select_first([RenameSVTruthIds.renamed_vcf_idx, truth_sv_vcf_subsetted_idx])

        if (do_exact) {
            if (defined(shard_bin_size_exact_match)) {
                call Helpers.CreateShardsFromVcfIndex as CreateExactShards {
                    input:
                        vcf_idx = vcf_final_idx,
                        ref_fai = ref_fai,
                        shard_bin_size = select_first([shard_bin_size_exact_match]),
                        prefix = "~{prefix}.~{contig}.exact_shards",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_create_exact_shards
                }

                scatter (k in range(length(CreateExactShards.shard_regions))) {
                    call Helpers.SubsetVcfToRegion as SubsetExactEval {
                        input:
                            vcf = vcf_final,
                            vcf_idx = vcf_final_idx,
                            region = CreateExactShards.shard_regions[k],
                            prefix = "~{prefix}.~{contig}.exact_eval_~{k}",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_subset_exact_vcf
                    }

                    call Helpers.SubsetVcfToRegion as SubsetExactTruth {
                        input:
                            vcf = truth_snv_indel_vcf_final,
                            vcf_idx = truth_snv_indel_vcf_final_idx,
                            region = CreateExactShards.shard_regions[k],
                            prefix = "~{prefix}.~{contig}.exact_truth_~{k}",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_subset_exact_truth
                    }

                    call ExactMatch as ExactMatchShard {
                        input:
                            vcf = SubsetExactEval.subset_vcf,
                            vcf_idx = SubsetExactEval.subset_vcf_idx,
                            truth_snv_indel_vcf = SubsetExactTruth.subset_vcf,
                            truth_snv_indel_vcf_idx = SubsetExactTruth.subset_vcf_idx,
                            source_tag = source_tag_truth_snv_indel_vcf,
                            prefix = "~{prefix}.~{contig}.exact_~{k}",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_exact_match
                    }

                    call AppendAnnotationsFromVcf as AppendExactAnnotationsShard {
                        input:
                            annotation_tsv = ExactMatchShard.annotation_tsv,
                            truth_vcf = SubsetExactTruth.subset_vcf,
                            truth_vcf_idx = SubsetExactTruth.subset_vcf_idx,
                            is_sv_truth = false,
                            prefix = "~{prefix}.~{contig}.exact_annotated_~{k}",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_append_exact_annotations
                    }
                }

                call Helpers.ConcatTsvs as ConcatExactAnnotations {
                    input:
                        tsvs = AppendExactAnnotationsShard.annotated_tsv,
                        sort_output = true,
                        prefix = "~{prefix}.~{contig}.exact_annotations",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_concat_exact_annotations
                }

                call Helpers.ConcatVcfs as ConcatExactUnmatched {
                    input:
                        vcfs = ExactMatchShard.unmatched_vcf,
                        vcf_idxs = ExactMatchShard.unmatched_vcf_idx,
                        allow_overlaps = false,
                        naive = false,
                        prefix = "~{prefix}.~{contig}.exact_unmatched",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_concat_exact_unmatched
                }
            }

            if (!defined(shard_bin_size_exact_match)) {
                call ExactMatch as ExactMatchFull {
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

                call AppendAnnotationsFromVcf as AppendExactAnnotationsFull {
                    input:
                        annotation_tsv = ExactMatchFull.annotation_tsv,
                        truth_vcf = ExactMatchFull.matched_truth_vcf,
                        truth_vcf_idx = ExactMatchFull.matched_truth_vcf_idx,
                        is_sv_truth = false,
                        prefix = "~{prefix}.~{contig}.exact_annotated",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_append_exact_annotations
                }
            }

            File exact_extended_tsv = select_first([ConcatExactAnnotations.concatenated_tsv, AppendExactAnnotationsFull.annotated_tsv])
            File exact_unmatched_vcf = select_first([ConcatExactUnmatched.concat_vcf, ExactMatchFull.unmatched_vcf])
            File exact_unmatched_vcf_idx = select_first([ConcatExactUnmatched.concat_vcf_idx, ExactMatchFull.unmatched_vcf_idx])
        }

        File vcf_post_exact = select_first([exact_unmatched_vcf, vcf_final])
        File vcf_post_exact_idx = select_first([exact_unmatched_vcf_idx, vcf_final_idx])

        if (do_truvari) {
            call TruvariMatch.TruvariMatch {
                input:
                    vcf = vcf_post_exact,
                    vcf_idx = vcf_post_exact_idx,
                    truth_snv_indel_vcf = truth_snv_indel_vcf_final,
                    truth_snv_indel_vcf_idx = truth_snv_indel_vcf_final_idx,
                    prefix = "~{prefix}.~{contig}.truvari",
                    min_sv_length = min_sv_length_truvari_vcf,
                    min_sv_length_truth = min_sv_length_truvari_truth_vcf,
                    length_field = length_field_vcf,
                    source_tag = source_tag_truth_snv_indel_vcf,
                    ref_fa = ref_fa,
                    ref_fai = ref_fai,
                    utils_docker = utils_docker,
                    runtime_attr_subset_vcf = runtime_attr_truvari_subset_vcf,
                    runtime_attr_subset_truth = runtime_attr_truvari_subset_truth,
                    runtime_attr_run_truvari_09 = runtime_attr_truvari_run_truvari_09,
                    runtime_attr_run_truvari_07 = runtime_attr_truvari_run_truvari_07,
                    runtime_attr_run_truvari_05 = runtime_attr_truvari_run_truvari_05,
                    runtime_attr_concat_matched = runtime_attr_truvari_concat_matched,
                    runtime_attr_concat_matched_truth = runtime_attr_truvari_concat_matched_truth
            }

            call AppendAnnotationsFromVcf as AppendTruvariAnnotations {
                input:
                    annotation_tsv = TruvariMatch.annotation_tsv,
                    truth_vcf = TruvariMatch.matched_truth_vcf,
                    truth_vcf_idx = TruvariMatch.matched_truth_vcf_idx,
                    is_sv_truth = false,
                    prefix = "~{prefix}.~{contig}.truvari_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_append_truvari_annotations
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
                    min_sv_length = min_sv_length_bedtools_closest_vcf,
                    min_sv_length_truth = min_sv_length_bedtools_closest_truth_vcf,
                    type_field = type_field_vcf,
                    length_field = length_field_vcf,
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

            call AppendAnnotationsFromVcf as AppendBedtoolsAnnotations {
                input:
                    annotation_tsv = BedtoolsClosestSV.annotation_tsv,
                    truth_vcf = truth_sv_vcf_final,
                    truth_vcf_idx = truth_sv_vcf_final_idx,
                    is_sv_truth = true,
                    prefix = "~{prefix}.~{contig}.bedtools_annotated",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_append_bedtools_annotations
            }
        }

        Array[File] extended_annotation_tsvs = select_all([
            exact_extended_tsv,
            AppendTruvariAnnotations.annotated_tsv,
            AppendBedtoolsAnnotations.annotated_tsv,
        ])

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
                preserve_header = true,
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

        bgzip -c isec_matched/0001.vcf > ~{prefix}.matched_truth.vcf.gz
        tabix -p vcf ~{prefix}.matched_truth.vcf.gz

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
        File matched_truth_vcf = "~{prefix}.matched_truth.vcf.gz"
        File matched_truth_vcf_idx = "~{prefix}.matched_truth.vcf.gz.tbi"
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

task AppendAnnotationsFromVcf {
    input {
        File annotation_tsv
        File truth_vcf
        File truth_vcf_idx
        Boolean is_sv_truth
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'EOF'
import subprocess
import re

annotation_tsv = "~{annotation_tsv}"
truth_vcf = "~{truth_vcf}"
is_sv_truth = ~{true="True" false="False" is_sv_truth}
prefix = "~{prefix}"

def get_ac_af_an_fields(vcf_path):
    cmd = f"bcftools view -h {vcf_path}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    fields = {'AC': [], 'AF': [], 'AN': []}
    for line in proc.stdout:
        m = re.match(r'##INFO=<ID=([^,]+)', line)
        if m:
            fid = m.group(1)
            for p in ['AC', 'AF', 'AN']:
                if fid == p or fid.startswith(p + '_'):
                    fields[p].append(fid)
    proc.wait()
    return fields

vcf_fields = get_ac_af_an_fields(truth_vcf)
dyn_cols = sorted(vcf_fields['AC']) + sorted(vcf_fields['AF']) + sorted(vcf_fields['AN'])

if is_sv_truth:
    extra_fields = ['N_HOMREF', 'N_HET', 'N_HOMALT']
else:
    extra_fields = [
        'dp_hist_all_bin_freq', 'dp_hist_all_n_larger',
        'ab_hist_alt_bin_freq', 'dp_hist_alt_n_larger',
        'age_hist_het_bin_freq', 'age_hist_het_n_smaller', 'age_hist_het_n_larger',
        'nhomalt',
    ]

query_fields = list(vcf_fields['AC']) + list(vcf_fields['AF']) + list(vcf_fields['AN']) + extra_fields
fmt = '%ID\\t' + '\\t'.join(f'%INFO/{f}' for f in query_fields) + '\\n'
cmd = f"bcftools query -f '{fmt}' {truth_vcf}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
truth_info = {}
for line in proc.stdout:
    parts = line.rstrip('\n').split('\t')
    if len(parts) == len(query_fields) + 1:
        truth_info[parts[0]] = {query_fields[i]: parts[i + 1] for i in range(len(query_fields))}
proc.wait()

def parse_hist(val):
    if not val or val == '.':
        return 0
    try:
        return sum(float(x) for x in val.split('|') if x and x != '.')
    except Exception:
        return 0

def to_num(val):
    try:
        return float(val) if val and val != '.' else 0
    except Exception:
        return 0

def compute_genotype_counts(info):
    if is_sv_truth:
        return info.get('N_HOMREF', '.'), info.get('N_HET', '.'), info.get('N_HOMALT', '.')
    homref = str(int(
        parse_hist(info.get('dp_hist_all_bin_freq', '.')) + to_num(info.get('dp_hist_all_n_larger', '.')) -
        parse_hist(info.get('ab_hist_alt_bin_freq', '.')) - to_num(info.get('dp_hist_alt_n_larger', '.'))
    ))
    het = str(int(
        parse_hist(info.get('age_hist_het_bin_freq', '.')) +
        to_num(info.get('age_hist_het_n_smaller', '.')) +
        to_num(info.get('age_hist_het_n_larger', '.'))
    ))
    return homref, het, info.get('nhomalt', '.')

extra_cols = ['match_type', 'truth_ID', 'source_tag', 'filter'] + dyn_cols + ['homref_count', 'het_count', 'homalt_count']
header_row = '\t'.join(['#CHROM', 'POS', 'REF', 'ALT', 'ID'] + extra_cols)

with open(annotation_tsv) as fin, open(f"{prefix}.tsv", 'w') as fout:
    fout.write(header_row + '\n')
    for line in fin:
        fields = line.rstrip('\n').split('\t')
        truth_id = fields[6]
        info = truth_info.get(truth_id, {})
        dyn_vals = [info.get(f, '.') for f in dyn_cols]
        homref, het, homalt = compute_genotype_counts(info)
        fout.write('\t'.join(fields + dyn_vals + [homref, het, homalt]) + '\n')

EOF
    >>>

    output {
        File annotated_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 25,
        disk_gb: 2 * ceil(size(annotation_tsv, "GB") + size(truth_vcf, "GB")) + 10,
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
genotype_cols = ['homref_count', 'het_count', 'homalt_count']
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
    fout.write('\t'.join(master_header) + '\n')
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
