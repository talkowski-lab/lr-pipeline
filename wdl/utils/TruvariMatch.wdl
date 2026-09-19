version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TruvariMatch {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        String contig
        String prefix

        String source_tag = "SNV_indel"

        Int? shard_bin_size_truvari_match
        Int min_shard_gap_truvari_match = 10000

        File? ref_fa
        File? ref_fai

        String utils_docker

        RuntimeAttr? runtime_attr_create_truvari_shards
        RuntimeAttr? runtime_attr_subset_truvari_vcf
        RuntimeAttr? runtime_attr_subset_truvari_truth
        RuntimeAttr? runtime_attr_run_truvari_09
        RuntimeAttr? runtime_attr_run_truvari_07
        RuntimeAttr? runtime_attr_run_truvari_05
        RuntimeAttr? runtime_attr_concat_matched
        RuntimeAttr? runtime_attr_concat_matched_truth
        RuntimeAttr? runtime_attr_concat_unmatched
    }

    Boolean shard_truvari = defined(shard_bin_size_truvari_match)

    if (shard_truvari) {
        call Helpers.CreateGapAwareShards as CreateTruvariShards {
            input:
                vcfs = [vcf, truth_snv_indel_vcf],
                vcf_idxs = [vcf_idx, truth_snv_indel_vcf_idx],
                contig = contig,
                shard_bin_size = select_first([shard_bin_size_truvari_match]),
                min_gap = min_shard_gap_truvari_match,
                prefix = "~{prefix}.truvari_shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_truvari_shards
        }
    }

    Array[String] truvari_regions = select_first([CreateTruvariShards.shard_regions, [contig]])

    scatter (k in range(length(truvari_regions))) {
        if (shard_truvari) {
            call Helpers.SubsetVcfToRegionStreaming as SubsetTruvariEvalRegion {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    region = truvari_regions[k],
                    prefix = "~{prefix}.truvari_eval_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_truvari_vcf
            }

            call Helpers.SubsetVcfToRegionStreaming as SubsetTruvariTruthRegion {
                input:
                    vcf = truth_snv_indel_vcf,
                    vcf_idx = truth_snv_indel_vcf_idx,
                    region = truvari_regions[k],
                    prefix = "~{prefix}.truvari_truth_~{k}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_truvari_truth
            }
        }

        File region_vcf = select_first([SubsetTruvariEvalRegion.subset_vcf, vcf])
        File region_vcf_idx = select_first([SubsetTruvariEvalRegion.subset_vcf_idx, vcf_idx])
        File region_truth_vcf = select_first([SubsetTruvariTruthRegion.subset_vcf, truth_snv_indel_vcf])
        File region_truth_vcf_idx = select_first([SubsetTruvariTruthRegion.subset_vcf_idx, truth_snv_indel_vcf_idx])

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
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
