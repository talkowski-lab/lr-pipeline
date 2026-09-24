version 1.0

import "Helpers.wdl"
import "Structs.wdl"

workflow BedtoolsClosestSV {
    meta {
        description: [
            "This sub-workflow performs the final comparison round, pairing each still-unmatched record with its nearest truth-callset neighbour using `bedtools closest`. Insertions and CNVs are compared separately, since proximity means different things for each, and the two comparisons are merged into a single annotation table."
        ]
    }

    parameter_meta {
        vcf: "Records left unmatched by `TruvariMatch`."
        vcf_idx: "Index for vcf."
        truth_sv_vcf: "Truth SV callset."
        truth_sv_vcf_idx: "Index for truth_sv_vcf."
        min_sv_length: "Minimum SV length applied to the callset."
        min_sv_length_truth: "Minimum SV length applied to the truth callset."
        type_field: "INFO field in the callset VCF holding variant type."
        length_field: "INFO field in the callset VCF holding allele length."
        length_field_truth: "INFO field in the truth VCF holding allele length, used to apply `min_sv_length_truth`. The truth callset arrives here in symbolic form, so this is `SVLEN` unless the caller names it otherwise."
        move_dup_to_origin: "Whether canonical DUPs are repositioned onto their `INFO/ORIGIN` interval before the DUP-vs-DUP reciprocal-overlap comparison. When false each DUP instead spans its own coordinates, from POS over its allele length, and `INFO/ORIGIN` is not required."
        source_tag: "Tag identifying the truth callset in the annotations."
        annotation_tsv: "Nearest-neighbour annotations for the remaining records."
    }

    input {
        File vcf
        File vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        String prefix

        Int min_sv_length
        Int min_sv_length_truth
        String type_field
        String length_field
        String length_field_truth = "SVLEN"
        Boolean move_dup_to_origin = true
        String source_tag = "SV"

        String gatk_sv_lr_docker
        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_truth
        RuntimeAttr? runtime_attr_convert_to_symbolic
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_split_truth
        RuntimeAttr? runtime_attr_compare
        RuntimeAttr? runtime_attr_calculate
        RuntimeAttr? runtime_attr_merge_comparisons
    }

    # Subset the callset to variants at or above the minimum SV length
    call Helpers.SubsetVcfByLength as SubsetEval {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            length_field = length_field,
            min_length = min_sv_length,
            prefix = "~{prefix}.subset_eval",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_vcf
    }

    # Convert the callset to symbolic alleles and split it by type, leaving each DUP on its own coordinates
    call Helpers.ConvertToSymbolic as ConvertEvalUnmoved {
        input:
            vcf = SubsetEval.subset_vcf,
            vcf_idx = SubsetEval.subset_vcf_idx,
            move_dup_to_origin = false,
            type_field = type_field,
            length_field = length_field,
            prefix = "~{prefix}.eval.symbolic.unmoved",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_convert_to_symbolic
    }

    call SplitVcf as SplitEvalUnmoved {
        input:
            vcf = ConvertEvalUnmoved.processed_vcf,
            vcf_idx = ConvertEvalUnmoved.processed_vcf_idx,
            split_cpx = false,
            prefix = "~{prefix}.eval.unmoved",
            docker = gatk_sv_lr_docker,
            runtime_attr_override = runtime_attr_split_vcf
    }

    # Convert and split the callset again with canonical DUPs repositioned onto their ORIGIN coordinates
    if (move_dup_to_origin) {
        call Helpers.ConvertToSymbolic as ConvertEvalMoved {
            input:
                vcf = SubsetEval.subset_vcf,
                vcf_idx = SubsetEval.subset_vcf_idx,
                move_dup_to_origin = true,
                type_field = type_field,
                length_field = length_field,
                prefix = "~{prefix}.eval.symbolic.moved",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_convert_to_symbolic
        }

        call SplitVcf as SplitEvalMoved {
            input:
                vcf = ConvertEvalMoved.processed_vcf,
                vcf_idx = ConvertEvalMoved.processed_vcf_idx,
                split_cpx = false,
                prefix = "~{prefix}.eval.moved",
                docker = gatk_sv_lr_docker,
                runtime_attr_override = runtime_attr_split_vcf
        }
    }

    File eval_dup_bed = select_first([SplitEvalMoved.dup_bed, SplitEvalUnmoved.dup_bed])

    # Subset the truth callset to variants at or above its own minimum SV length
    call Helpers.SubsetVcfByLength as SubsetTruth {
        input:
            vcf = truth_sv_vcf,
            vcf_idx = truth_sv_vcf_idx,
            length_field = length_field_truth,
            min_length = min_sv_length_truth,
            prefix = "~{prefix}.subset_truth",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_subset_truth
    }

    # Split the truth callset by type, breaking complex records into their constituent intervals
    call SplitVcf as SplitTruth {
        input:
            vcf = SubsetTruth.subset_vcf,
            vcf_idx = SubsetTruth.subset_vcf_idx,
            split_cpx = true,
            prefix = "~{prefix}.truth",
            docker = gatk_sv_lr_docker,
            runtime_attr_override = runtime_attr_split_truth
    }

    # Compare DEL in the callset to DEL in the truth callset by reciprocal overlap
    call Helpers.BedtoolsClosest as CompareDEL {
        input:
            bed_a = SplitEvalUnmoved.del_bed,
            bed_b = SplitTruth.del_bed,
            prefix = "~{prefix}.DEL",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedSVs as CalcuDEL {
        input:
            input_bed = CompareDEL.output_bed,
            prefix = "~{prefix}.DEL",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    # Compare INS in the callset to INS in the truth callset by breakpoint proximity and length ratio
    call Helpers.BedtoolsClosest as CompareINS {
        input:
            bed_a = SplitEvalUnmoved.ins_bed,
            bed_b = SplitTruth.ins_bed,
            prefix = "~{prefix}.INS",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedINSs as CalcuINS {
        input:
            input_bed = CompareINS.output_bed,
            prefix = "~{prefix}.INS",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    # Compare DUP in the callset to DUP in the truth callset by reciprocal overlap, on ORIGIN coordinates when moved
    call Helpers.BedtoolsClosest as CompareDUP {
        input:
            bed_a = eval_dup_bed,
            bed_b = SplitTruth.dup_bed,
            prefix = "~{prefix}.DUP",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedSVs as CalcuDUP {
        input:
            input_bed = CompareDUP.output_bed,
            prefix = "~{prefix}.DUP",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    # Compare INS in the callset to truth DUP collapsed to a point, catching insertions the truth callset typed as DUP
    call CollapseRangedToPoint as CollapseTruthDUP {
        input:
            bed = SplitTruth.dup_bed,
            prefix = "~{prefix}.truth.DUP_as_point",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.BedtoolsClosest as CompareINSDUP {
        input:
            bed_a = SplitEvalUnmoved.ins_bed,
            bed_b = CollapseTruthDUP.point_bed,
            prefix = "~{prefix}.INS_DUP",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedINSs as CalcuINSDUP {
        input:
            input_bed = CompareINSDUP.output_bed,
            prefix = "~{prefix}.INS_DUP",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    # Compare callset DUP collapsed to a point at its insertion site to INS in the truth callset
    call CollapseRangedToPoint as CollapseEvalDUP {
        input:
            bed = SplitEvalUnmoved.dup_bed,
            prefix = "~{prefix}.eval.DUP_as_point",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    call Helpers.BedtoolsClosest as CompareDUPINS {
        input:
            bed_a = CollapseEvalDUP.point_bed,
            bed_b = SplitTruth.ins_bed,
            prefix = "~{prefix}.DUP_INS",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_compare
    }

    call SelectMatchedINSs as CalcuDUPINS {
        input:
            input_bed = CompareDUPINS.output_bed,
            prefix = "~{prefix}.DUP_INS",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_calculate
    }

    # Merge the comparisons, keeping same-type matches ahead of the cross-type fallbacks
    call PrioritizedConcatComparisons {
        input:
            primary_tsvs = [CalcuDEL.output_comp, CalcuINS.output_comp, CalcuDUP.output_comp],
            secondary_tsvs = [CalcuINSDUP.output_comp, CalcuDUPINS.output_comp],
            prefix = "~{prefix}.comparison",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    # Join the retained matches back to the callset and truth records to emit the annotation TSV
    call CreateBedtoolsAnnotationTsv {
        input:
            truvari_unmatched_vcf = SubsetEval.subset_vcf,
            truvari_unmatched_vcf_idx = SubsetEval.subset_vcf_idx,
            truth_sv_vcf = SubsetTruth.subset_vcf,
            truth_sv_vcf_idx = SubsetTruth.subset_vcf_idx,
            closest_bed = PrioritizedConcatComparisons.merged_tsv,
            source_tag = source_tag,
            prefix = "~{prefix}.bedtools_closest_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_merge_comparisons
    }

    output {
        File annotation_tsv = CreateBedtoolsAnnotationTsv.annotation_tsv
    }
}

task SplitVcf {
    input {
        File vcf
        File vcf_idx
        Boolean split_cpx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        svtk vcf2bed \
            -i SVTYPE \
            -i SVLEN \
            ~{if split_cpx then "--split-cpx" else ""} \
            ~{vcf} \
            tmp.bed

        cut -f1-4,7-8 tmp.bed > ~{prefix}.bed

        set +o pipefail

        head -1 ~{prefix}.bed > header

        set -o pipefail

        cat header <(awk '$5 == "DEL"' ~{prefix}.bed) > ~{prefix}.DEL.bed
        cat header <(awk '$5 == "DUP"' ~{prefix}.bed) > ~{prefix}.DUP.bed
        cat header <(awk '$5 ~ /^INS/' ~{prefix}.bed) > ~{prefix}.INS.bed
        cat header <(awk '$5 == "INV" || $5 == "CPX"' ~{prefix}.bed) > ~{prefix}.INV_CPX.bed
        cat header <(awk '$5 == "BND" || $5 == "CTX"' ~{prefix}.bed) > ~{prefix}.BND_CTX.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
        File del_bed = "~{prefix}.DEL.bed"
        File dup_bed = "~{prefix}.DUP.bed"
        File ins_bed = "~{prefix}.INS.bed"
        File inv_bed = "~{prefix}.INV_CPX.bed"
        File bnd_bed = "~{prefix}.BND_CTX.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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

task SelectMatchedSVs {
    input {
        File input_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        Rscript /opt/scripts/benchmark/R1.bedtools_closest_CNV.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task SelectMatchedINSs {
    input {
        File input_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        Rscript /opt/scripts/benchmark/R2.bedtools_closest_INS.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task CollapseRangedToPoint {
    input {
        File bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        head -1 ~{bed} > ~{prefix}.bed

        awk 'NR>1 {OFS="\t"; $3=$2+1; print}' ~{bed} >> ~{prefix}.bed
    >>>

    output {
        File point_bed = "~{prefix}.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(bed, "GB")) + 5,
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

task PrioritizedConcatComparisons {
    input {
        Array[File] primary_tsvs
        Array[File] secondary_tsvs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{sep=" " primary_tsvs} > ~{prefix}.comparison

        { grep -v "query_svid" ~{prefix}.comparison || true; } | cut -f1 | sort -u > seen_ids.txt

        for f in ~{sep=" " secondary_tsvs}; do
            awk -F'\t' 'NR==FNR{seen[$1]=1; next} !($1 in seen)' seen_ids.txt "$f" >> ~{prefix}.comparison
        done
    >>>

    output {
        File merged_tsv = "~{prefix}.comparison"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 10,
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

task CreateBedtoolsAnnotationTsv {
    input {
        File truvari_unmatched_vcf
        File truvari_unmatched_vcf_idx
        File truth_sv_vcf
        File truth_sv_vcf_idx
        File closest_bed
        String source_tag
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export LC_ALL=C

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' ~{truvari_unmatched_vcf} \
        | sort -k5,5 > vcf_info.tsv

        bcftools query \
            -f '%ID\t%FILTER\n' ~{truth_sv_vcf} \
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
        }' \
        | sort -k1,1 > truth_filters.tsv

        { grep -v "query_svid" ~{closest_bed} || true; } \
            | awk -F'\t' '$1 != "" {print $1"\t"$2}' \
            | sort -k1,1 > matched_ids.tsv

        join \
            -1 5 \
            -2 1 \
            -t $'\t' \
            vcf_info.tsv \
            matched_ids.tsv \
            > joined.tsv

        sort -k6,6 joined.tsv > joined_sorted.tsv

        join \
            -1 6 \
            -2 1 \
            -t $'\t' \
            -a 1 \
            -e "." \
            -o '1.1,1.2,1.3,1.4,1.5,1.6,2.2' \
            joined_sorted.tsv \
            truth_filters.tsv \
            > joined_with_filter.tsv

        awk -F'\t' -v src="~{source_tag}" 'BEGIN{OFS="\t"} {
            print $2, $3, $4, $5, $1, "BEDTOOLS_CLOSEST", $6, src, $7
        }' joined_with_filter.tsv \
        | sort -k1,1V -k2,2n > ~{prefix}.bedtools_matched.tsv
    >>>

    output {
        File annotation_tsv = "~{prefix}.bedtools_matched.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(truvari_unmatched_vcf, "GB") + size(truth_sv_vcf, "GB") + size(closest_bed, "GB")) + 5,
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
