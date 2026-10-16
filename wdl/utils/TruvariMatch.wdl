version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TruvariMatch {
    meta {
        description: [
            "This sub-workflow performs the second comparison round, matching records left unmatched by `ExactMatch` with Truvari at three decreasing sequence-similarity thresholds (0.9, 0.7, 0.5). Each threshold only sees what the previous one failed to match, so a record is annotated with the strictest threshold that matched it."
        ]
    }

    parameter_meta {
        vcf: "Unmatched callset records from `ExactMatch`."
        vcf_idx: "Index for vcf."
        truth_snv_indel_vcf: "Unmatched truth records."
        truth_snv_indel_vcf_idx: "Index for truth_snv_indel_vcf."
        source_tag: "Tag identifying the truth callset in the annotations."
        shard_records_truvari_match: "If set, shards the matching step into consecutive runs of this many callset records each. Each shard is compared against the truth records whose POS falls within its span widened by `shard_flank_truvari_match`."
        shard_flank_truvari_match: "Base pairs added to each side of a shard's span when selecting its truth records. Truvari only pairs records within 500 bp of each other by default, so the default flank keeps matches at shard edges."
        ref_fa: "Reference FASTA and index, when Truvari is run with reference context."
        ref_fai: "Index for ref_fa."
        annotation_tsv: "Truvari match annotations across all three thresholds."
        matched_truth_vcf: "Truth records matched by any threshold."
        matched_truth_vcf_idx: "Index for matched_truth_vcf."
        unmatched_vcf: "Records still unmatched after 0.5, passed to `BedtoolsClosestSV`."
        unmatched_vcf_idx: "Index for unmatched_vcf."
    }

    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String prefix

        String source_tag = "SNV_indel"

        Int? shard_records_truvari_match
        Int shard_flank_truvari_match = 1000

        File? ref_fa
        File? ref_fai

        String utils_docker

        RuntimeAttr? runtime_attr_create_truvari_shards
        RuntimeAttr? runtime_attr_subset_truvari_truth
        RuntimeAttr? runtime_attr_run_truvari_09
        RuntimeAttr? runtime_attr_run_truvari_07
        RuntimeAttr? runtime_attr_run_truvari_05
        RuntimeAttr? runtime_attr_concat_matched
        RuntimeAttr? runtime_attr_concat_matched_truth
        RuntimeAttr? runtime_attr_concat_unmatched
    }

    Boolean shard_truvari = defined(shard_records_truvari_match)

    if (shard_truvari) {
        call Helpers.ShardVcfByRecords as ShardTruvariEval {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                records_per_shard = select_first([shard_records_truvari_match]),
                prefix = "~{prefix}.truvari_eval",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_truvari_shards
        }
    }

    # A callset with no records yields no shards, so it runs whole as a single shard instead
    Array[File] sharded_vcfs = select_first([ShardTruvariEval.shards, []])
    Boolean use_shards = length(sharded_vcfs) > 0
    Array[File] eval_vcfs = if use_shards then sharded_vcfs else [vcf]
    Array[File] eval_vcf_idxs = if use_shards then select_first([ShardTruvariEval.shard_idxs]) else [vcf_idx]

    scatter (k in range(length(eval_vcfs))) {
        if (use_shards) {
            call Helpers.SubsetVcfToRecordSpan as SubsetTruvariTruthSpan {
                input:
                    vcf = truth_snv_indel_vcf,
                    vcf_idx = truth_snv_indel_vcf_idx,
                    span_vcf = eval_vcfs[k],
                    flank = shard_flank_truvari_match,
                    prefix = "~{prefix}.truvari_truth_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_truvari_truth
            }
        }

        File region_vcf = eval_vcfs[k]
        File region_vcf_idx = eval_vcf_idxs[k]
        File region_truth_vcf = select_first([SubsetTruvariTruthSpan.subset_vcf, truth_snv_indel_vcf])
        File region_truth_vcf_idx = select_first([SubsetTruvariTruthSpan.subset_vcf_idx, truth_snv_indel_vcf_idx])

        call RunTruvari as RunTruvari09 {
            input:
                vcf = region_vcf,
                vcf_idx = region_vcf_idx,
                truth_snv_indel_vcf = region_truth_vcf,
                truth_snv_indel_vcf_idx = region_truth_vcf_idx,
                pctseq = 0.9,
                pctsize = 0.9,
                pctovl = 0.9,
                sizemin = 0,
                sizefilt = 0,
                tag_value = "TRUVARI_0.9",
                source_tag = source_tag,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.0.9_~{k}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_run_truvari_09
        }

        call RunTruvari as RunTruvari07 {
            input:
                vcf = RunTruvari09.unmatched_vcf,
                vcf_idx = RunTruvari09.unmatched_vcf_idx,
                truth_snv_indel_vcf = region_truth_vcf,
                truth_snv_indel_vcf_idx = region_truth_vcf_idx,
                pctseq = 0.7,
                pctsize = 0.7,
                pctovl = 0.7,
                sizemin = 0,
                sizefilt = 0,
                tag_value = "TRUVARI_0.7",
                source_tag = source_tag,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.0.7_~{k}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_run_truvari_07
        }

        call RunTruvari as RunTruvari05 {
            input:
                vcf = RunTruvari07.unmatched_vcf,
                vcf_idx = RunTruvari07.unmatched_vcf_idx,
                truth_snv_indel_vcf = region_truth_vcf,
                truth_snv_indel_vcf_idx = region_truth_vcf_idx,
                pctseq = 0.5,
                pctsize = 0.5,
                pctovl = 0.5,
                sizemin = 0,
                sizefilt = 0,
                tag_value = "TRUVARI_0.5",
                source_tag = source_tag,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.0.5_~{k}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_run_truvari_05
        }
    }

    call Helpers.ConcatTsvs as ConcatAnnotationTsvs {
        input:
            tsvs = flatten([RunTruvari09.annotation_tsv, RunTruvari07.annotation_tsv, RunTruvari05.annotation_tsv]),
            sort_output = true,
            prefix = "~{prefix}.truvari_combined",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_matched
    }

    call Helpers.ConcatVcfs as ConcatMatchedTruth {
        input:
            vcfs = flatten([RunTruvari09.matched_truth_vcf, RunTruvari07.matched_truth_vcf, RunTruvari05.matched_truth_vcf]),
            vcf_idxs = flatten([RunTruvari09.matched_truth_vcf_idx, RunTruvari07.matched_truth_vcf_idx, RunTruvari05.matched_truth_vcf_idx]),
            allow_overlaps = true,
            naive = false,
            sort_output = true,
            prefix = "~{prefix}.matched_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_matched_truth
    }

    if (shard_truvari) {
        call Helpers.ConcatVcfs as ConcatUnmatched {
            input:
                vcfs = RunTruvari05.unmatched_vcf,
                vcf_idxs = RunTruvari05.unmatched_vcf_idx,
                allow_overlaps = false,
                naive = false,
                prefix = "~{prefix}.unmatched",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_unmatched
        }
    }

    output {
        File annotation_tsv = ConcatAnnotationTsvs.concatenated_tsv
        File matched_truth_vcf = ConcatMatchedTruth.concat_vcf
        File matched_truth_vcf_idx = ConcatMatchedTruth.concat_vcf_idx
        File unmatched_vcf = select_first([ConcatUnmatched.concat_vcf, RunTruvari05.unmatched_vcf[0]])
        File unmatched_vcf_idx = select_first([ConcatUnmatched.concat_vcf_idx, RunTruvari05.unmatched_vcf_idx[0]])
    }
}

task RunTruvari {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        Float pctseq
        Float pctsize
        Float pctovl
        Int sizemin
        Int sizefilt
        String tag_value
        String source_tag
        File? ref_fa
        File? ref_fai
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ -d "~{prefix}_truvari" ]; then
            rm -r "~{prefix}_truvari"
        fi

        # Truvari looks for each index beside its VCF, but Cromwell localizes separately globbed shard indexes elsewhere
        if [[ "~{vcf_idx}" != "~{vcf}.tbi" ]]; then
            ln -sf "~{vcf_idx}" "~{vcf}.tbi"
        fi
        if [[ "~{truth_snv_indel_vcf_idx}" != "~{truth_snv_indel_vcf}.tbi" ]]; then
            ln -sf "~{truth_snv_indel_vcf_idx}" "~{truth_snv_indel_vcf}.tbi"
        fi

        # Pass --reference only when supplied, because truvari needs it solely to resolve symbolic ALTs
        truvari bench \
            -b ~{truth_snv_indel_vcf} \
            -c ~{vcf} \
            -o "~{prefix}_truvari" \
            ~{if defined(ref_fa) then "--reference " + ref_fa else ""} \
            --pctseq ~{pctseq} \
            --pctsize ~{pctsize} \
            --pctovl ~{pctovl} \
            --sizemin ~{sizemin} \
            --sizefilt ~{sizefilt}

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/MatchId\n' "~{prefix}_truvari/tp-comp.vcf.gz" \
            | awk 'BEGIN{FS=OFS="\t"} {split($6,a,","); print $1,$2,$3,$4,$5,a[1]}' \
            | LC_ALL=C sort -k6,6 > comp.mid.tsv

        bcftools query -f '%ID\t%INFO/MatchId\t%FILTER\n' "~{prefix}_truvari/tp-base.vcf.gz" \
            | awk 'BEGIN{FS=OFS="\t"} {
                split($2,a,",")
                n = split($3, parts, ";")
                out = ""
                for (i = 1; i <= n; i++) {
                    if (parts[i] != "." && parts[i] != "PASS") {
                        out = (out == "" ? parts[i] : out "," parts[i])
                    }
                }
                if (out == "") out = "."
                print a[1],$1,out
            }' \
            | LC_ALL=C sort -k1,1 > base.mid2id.tsv

        LC_ALL=C join -t $'\t' -1 6 -2 1 comp.mid.tsv base.mid2id.tsv \
            | awk -F'\t' -v tag="~{tag_value}" -v src="~{source_tag}" 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6,tag,$7,src,$8}' \
            > ~{prefix}.annotation.tsv
    >>>

    output {
        File annotation_tsv = "~{prefix}.annotation.tsv"
        File matched_truth_vcf = "~{prefix}_truvari/tp-base.vcf.gz"
        File matched_truth_vcf_idx = "~{prefix}_truvari/tp-base.vcf.gz.tbi"
        File unmatched_vcf = "~{prefix}_truvari/fp.vcf.gz"
        File unmatched_vcf_idx = "~{prefix}_truvari/fp.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10 * ceil(size(vcf, "GB") + size(truth_snv_indel_vcf, "GB") + size(ref_fa, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: if tag_value == "TRUVARI_0.9" then "16 GiB" else if tag_value == "TRUVARI_0.7" then "16 GiB" else "16 GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
