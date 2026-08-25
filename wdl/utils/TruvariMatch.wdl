version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow TruvariMatch {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File ref_fa
        File ref_fai
        String prefix

        String source_tag = "SNV_indel"

        String utils_docker

        RuntimeAttr? runtime_attr_run_truvari_09
        RuntimeAttr? runtime_attr_run_truvari_07
        RuntimeAttr? runtime_attr_run_truvari_05
        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_concat_matched
        RuntimeAttr? runtime_attr_concat_matched_truth
        RuntimeAttr? runtime_attr_concat_unmatched
    }

    call Helpers.SubsetVcfByArgs as SubsetDel {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            include_args = 'INFO/allele_type~"del"',
            prefix = "~{prefix}.del",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_vcf
    }

    call Helpers.SubsetVcfByArgs as SubsetNonDel {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            exclude_args = 'INFO/allele_type~"del"',
            prefix = "~{prefix}.non_del",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_vcf
    }

    call RunTruvari as RunTruvariDel_09 {
        input:
            vcf = SubsetDel.subset_vcf,
            vcf_idx = SubsetDel.subset_vcf_idx,
            truth_snv_indel_vcf = truth_snv_indel_vcf,
            truth_snv_indel_vcf_idx = truth_snv_indel_vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.9,
            pctsize = 0.9,
            pctovl = 0.9,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.9",
            source_tag = source_tag,
            prefix = "~{prefix}.0.9",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_09
    }

    call RunTruvari as RunTruvariNonDel_09 {
        input:
            vcf = SubsetNonDel.subset_vcf,
            vcf_idx = SubsetNonDel.subset_vcf_idx,
            truth_snv_indel_vcf = truth_snv_indel_vcf,
            truth_snv_indel_vcf_idx = truth_snv_indel_vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.9,
            pctsize = 0.9,
            pctovl = 0.0,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.9",
            source_tag = source_tag,
            prefix = "~{prefix}.non_del.0.9",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_09
    }

    call RunTruvari as RunTruvariDel_07 {
        input:
            vcf = RunTruvariDel_09.unmatched_vcf,
            vcf_idx = RunTruvariDel_09.unmatched_vcf_idx,
            truth_snv_indel_vcf = truth_snv_indel_vcf,
            truth_snv_indel_vcf_idx = truth_snv_indel_vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.7,
            pctsize = 0.7,
            pctovl = 0.7,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.7",
            source_tag = source_tag,
            prefix = "~{prefix}.0.7",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_07
    }

    call RunTruvari as RunTruvariNonDel_07 {
        input:
            vcf = RunTruvariNonDel_09.unmatched_vcf,
            vcf_idx = RunTruvariNonDel_09.unmatched_vcf_idx,
            truth_snv_indel_vcf = truth_snv_indel_vcf,
            truth_snv_indel_vcf_idx = truth_snv_indel_vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.7,
            pctsize = 0.7,
            pctovl = 0.0,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.7",
            source_tag = source_tag,
            prefix = "~{prefix}.non_del.0.7",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_07
    }

    call RunTruvari as RunTruvariDel_05 {
        input:
            vcf = RunTruvariDel_07.unmatched_vcf,
            vcf_idx = RunTruvariDel_07.unmatched_vcf_idx,
            truth_snv_indel_vcf = truth_snv_indel_vcf,
            truth_snv_indel_vcf_idx = truth_snv_indel_vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.5,
            pctsize = 0.5,
            pctovl = 0.5,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.5",
            source_tag = source_tag,
            prefix = "~{prefix}.0.5",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_05
    }

    call RunTruvari as RunTruvariNonDel_05 {
        input:
            vcf = RunTruvariNonDel_07.unmatched_vcf,
            vcf_idx = RunTruvariNonDel_07.unmatched_vcf_idx,
            truth_snv_indel_vcf = truth_snv_indel_vcf,
            truth_snv_indel_vcf_idx = truth_snv_indel_vcf_idx,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            pctseq = 0.5,
            pctsize = 0.5,
            pctovl = 0.0,
            sizemin = 0,
            sizefilt = 0,
            tag_value = "TRUVARI_0.5",
            source_tag = source_tag,
            prefix = "~{prefix}.non_del.0.5",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_run_truvari_05
    }

    call Helpers.ConcatTsvs as ConcatAnnotationTsvs {
        input:
            tsvs = [
                RunTruvariDel_09.annotation_tsv,
                RunTruvariNonDel_09.annotation_tsv,
                RunTruvariDel_07.annotation_tsv,
                RunTruvariNonDel_07.annotation_tsv,
                RunTruvariDel_05.annotation_tsv,
                RunTruvariNonDel_05.annotation_tsv
            ],
            sort_output = true,
            prefix = "~{prefix}.truvari_combined",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_matched
    }

    call Helpers.ConcatVcfs as ConcatMatchedTruth {
        input:
            vcfs = [
                RunTruvariDel_09.matched_truth_vcf,
                RunTruvariNonDel_09.matched_truth_vcf,
                RunTruvariDel_07.matched_truth_vcf,
                RunTruvariNonDel_07.matched_truth_vcf,
                RunTruvariDel_05.matched_truth_vcf,
                RunTruvariNonDel_05.matched_truth_vcf
            ],
            vcf_idxs = [
                RunTruvariDel_09.matched_truth_vcf_idx,
                RunTruvariNonDel_09.matched_truth_vcf_idx,
                RunTruvariDel_07.matched_truth_vcf_idx,
                RunTruvariNonDel_07.matched_truth_vcf_idx,
                RunTruvariDel_05.matched_truth_vcf_idx,
                RunTruvariNonDel_05.matched_truth_vcf_idx
            ],
            allow_overlaps = true,
            naive = false,
            sort_output = true,
            prefix = "~{prefix}.matched_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_matched_truth
    }

    call Helpers.ConcatVcfs as ConcatUnmatched {
        input:
            vcfs = [RunTruvariDel_05.unmatched_vcf, RunTruvariNonDel_05.unmatched_vcf],
            vcf_idxs = [RunTruvariDel_05.unmatched_vcf_idx, RunTruvariNonDel_05.unmatched_vcf_idx],
            allow_overlaps = true,
            naive = false,
            sort_output = true,
            prefix = "~{prefix}.unmatched",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_unmatched
    }

    output {
        File annotation_tsv = ConcatAnnotationTsvs.concatenated_tsv
        File matched_truth_vcf = ConcatMatchedTruth.concat_vcf
        File matched_truth_vcf_idx = ConcatMatchedTruth.concat_vcf_idx
        File unmatched_vcf = ConcatUnmatched.concat_vcf
        File unmatched_vcf_idx = ConcatUnmatched.concat_vcf_idx
    }
}

task RunTruvari {
    input {
        File vcf
        File vcf_idx
        File truth_snv_indel_vcf
        File truth_snv_indel_vcf_idx
        File ref_fa
        File ref_fai
        Float pctseq
        Float pctsize
        Float pctovl
        Int sizemin
        Int sizefilt
        String tag_value
        String source_tag
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        if [ -d "~{prefix}_truvari" ]; then
            rm -r "~{prefix}_truvari"
        fi

        truvari bench \
            -b ~{truth_snv_indel_vcf} \
            -c ~{vcf} \
            -o "~{prefix}_truvari" \
            --reference ~{ref_fa} \
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
        disk_gb: 10 * ceil(size(vcf, "GB") + size(truth_snv_indel_vcf, "GB")) + 10,
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
