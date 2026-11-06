version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

# Refine and filter per-sample gCNV depth calls from read depth before clustering
workflow FilterDepthCalls {
    meta {
        description: [
            "This workflow refines the breakpoints of per-sample gCNV depth calls from the binned read depth and removes calls whose depth profile does not support them, before the calls are clustered and genotyped. Each sample's call is refined and judged on that sample's depth alone; calls of other samples are never pooled or merged. Only DEL and DUP calls of at least `min_size` are kept, and calls overlapping `gaps_bed` are dropped.",
            "Calls without segmental duplication (SD) involvement are refined by locating the depth transition within 2 kb of each called breakpoint; when one breakpoint is sharp and the other is not, the other is searched for out to twice the call length. The sharpness of the two breakpoints (edge z-score, copy-number step and edge jump) is scored with a gradient-boosted tree model trained on chr6 of HPRC and HGSVC, with short-read GATK-SV calls as labels, and calls with a probability of at least `gbm_cutoff` are kept. Calls with a breakpoint in, or lying within, an SD are scored separately: each end must show the expected copy-number change against its flank and a sharp transition, the two ends evaluated together under `sd_sharp_rule`.",
            "The passing calls are written as a VCF with one record per call and carrier at the refined breakpoints, together with a table of the decision for every sample-call."
        ]
    }

    parameter_meta {
        depth_vcf: "Merged per-sample gCNV depth VCF, as `DepthPreprocessing` writes it."
        depth_vcf_idx: "Index for depth_vcf."
        contigs: "Contigs to refine and filter; calls on other contigs are dropped."
        min_size: "Minimum call length; shorter calls are dropped."
        gbm_cutoff: "Minimum breakpoint-sharpness model probability for a call without SD involvement to be kept."
        sd_sharp_rule: "How the sharp-transition checks of the two ends of an SD-involved call are combined: 'all' requires every check at both ends; 'majority' requires each end to pass 2 of the 3 checks and each check to pass at one end or more; 'score' requires the mean of the value-to-cutoff ratios to reach 1 with none below 0.8."
        rd_file: "Binned read-depth matrix of the cohort, as `CreateCohortDepthFiles` writes it."
        rd_file_idx: "Index for rd_file."
        median_coverage: "Per-sample median coverage table, whose sample order matches the columns of rd_file."
        seg_dup_bed: "From references."
        gaps_bed: "Reference assembly gaps, centromeres and telomeres; calls overlapping them are dropped."
        ref_fai: "Reference FASTA index, used to flag calls near contig ends."
        sharpness_model: "Breakpoint-sharpness model JSON; defaults to the copy shipped in the utils image."
        filtered_vcf: "Passing calls, one record per call and carrier, at the refined breakpoints."
        filtered_vcf_idx: "Index for filtered_vcf."
        decisions_tsv: "Decision for every sample-call of at least `min_size`: class, kept, original and refined breakpoints, model probability and SD status."
    }

    input {
        File depth_vcf
        File depth_vcf_idx
        Array[String] contigs
        String prefix

        Int min_size = 10000
        Float gbm_cutoff = 0.5
        String sd_sharp_rule = "majority"

        File rd_file
        File rd_file_idx
        File median_coverage
        File seg_dup_bed
        File gaps_bed
        File ref_fai
        File? sharpness_model

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_refine_filter
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_concat_tsvs
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = depth_vcf,
                vcf_idx = depth_vcf_idx,
                contig = contig,
                prefix = prefix + "." + contig,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call RefineFilterDepthCalls {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                contig = contig,
                min_size = min_size,
                gbm_cutoff = gbm_cutoff,
                sd_sharp_rule = sd_sharp_rule,
                rd_file = rd_file,
                rd_file_idx = rd_file_idx,
                median_coverage = median_coverage,
                seg_dup_bed = seg_dup_bed,
                gaps_bed = gaps_bed,
                ref_fai = ref_fai,
                sharpness_model = sharpness_model,
                prefix = prefix + "." + contig + ".depth_filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_refine_filter
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = RefineFilterDepthCalls.filtered_vcf,
            vcf_idxs = RefineFilterDepthCalls.filtered_vcf_idx,
            allow_overlaps = true,
            naive = false,
            sort_output = true,
            prefix = prefix + ".depth_filtered",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = RefineFilterDepthCalls.decisions_tsv,
            sort_output = true,
            preserve_header = true,
            compressed_tsvs = true,
            compressed_output = true,
            prefix = prefix + ".depth_filter_decisions",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    output {
        File filtered_vcf = ConcatVcfs.concat_vcf
        File filtered_vcf_idx = ConcatVcfs.concat_vcf_idx
        File decisions_tsv = ConcatTsvs.concatenated_tsv
    }
}

task RefineFilterDepthCalls {
    input {
        File vcf
        File vcf_idx
        String contig
        Int min_size
        Float gbm_cutoff
        String sd_sharp_rule
        File rd_file
        File rd_file_idx
        File median_coverage
        File seg_dup_bed
        File gaps_bed
        File ref_fai
        File? sharpness_model
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        rd_file: { localization_optional: true }
        rd_file_idx: { localization_optional: true }
    }

    command <<<
        set -euo pipefail

        # Stream this contig's rows of the read-depth matrix, so the cohort matrix is never localized in full
        if [[ "~{rd_file}" == gs://* ]]; then
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        fi
        tabix "~{rd_file}" "~{contig}" | bgzip -c > rd.txt.gz
        tabix -f -s 1 -b 2 -e 3 rd.txt.gz

        # Write one row per call and carrier for calls of at least the minimum size away from gaps
        python3 /opt/scripts/cnv/split_depth_vcf_calls.py \
            --vcf '~{vcf}' \
            --gaps '~{gaps_bed}' \
            --min-size ~{min_size} \
            --out-bed calls.bed.gz

        # Refine and score the calls only when the contig has any, since the scripts need at least one row
        n_calls=$(gzip -cd calls.bed.gz | awk 'NR > 1 {n++} END {print n + 0}')
        if [ "${n_calls}" -gt 0 ]; then
            python3 /opt/scripts/cnv/refine_cnv_breakpoints.py \
                --bed calls.bed.gz \
                --rd rd.txt.gz \
                --samples '~{median_coverage}' \
                --medians '~{median_coverage}' \
                --extend-search \
                --out-prefix calls.bp_refined

            python3 /opt/scripts/cnv/compute_breakpoint_sharpness.py \
                --calls calls.bp_refined.calls.tsv.gz \
                --bed calls.bed.gz \
                --rd rd.txt.gz \
                --samples '~{median_coverage}' \
                --medians '~{median_coverage}' \
                --segdup '~{seg_dup_bed}' \
                --out calls.bp_sharpness.tsv.gz

            python3 /opt/scripts/cnv/score_sharpness_gbm.py \
                --model '~{select_first([sharpness_model, "/opt/scripts/cnv/sharpness_gbm.json"])}' \
                --in calls.bp_sharpness.tsv.gz \
                --out calls.scored.tsv.gz

            python3 /opt/scripts/cnv/refine_sd_cnv_breakpoints.py \
                --bed calls.bed.gz \
                --rd rd.txt.gz \
                --samples '~{median_coverage}' \
                --medians '~{median_coverage}' \
                --segdup '~{seg_dup_bed}' \
                --fai '~{ref_fai}' \
                --no-merge \
                --local-flank \
                --require-sharp \
                --sharp-rule ~{sd_sharp_rule} \
                --out-prefix calls.sd_refined

            python3 /opt/scripts/cnv/filter_depth_cnv_calls.py \
                --vcf '~{vcf}' \
                --bed calls.bed.gz \
                --scored calls.scored.tsv.gz \
                --sd-members calls.sd_refined.members.tsv.gz \
                --sd-events calls.sd_refined.events.tsv.gz \
                --gaps '~{gaps_bed}' \
                --gbm-cutoff ~{gbm_cutoff} \
                --out-prefix '~{prefix}'
        else
            python3 /opt/scripts/cnv/filter_depth_cnv_calls.py \
                --vcf '~{vcf}' \
                --bed calls.bed.gz \
                --out-prefix '~{prefix}'
        fi
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File decisions_tsv = "~{prefix}.decisions.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 32,
        disk_gb: ceil(size(vcf, "GB") * 2) + 50,
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
