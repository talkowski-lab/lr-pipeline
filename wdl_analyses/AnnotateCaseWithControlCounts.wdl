version 1.0

import "Structs.wdl"
import "MatchControlsByAncestry.wdl"

# Annotate a case VCF with allele counts from a control call set, overall and in the ancestry-matched controls (e.g. the
# control_weights.tsv output of MatchControlsFromPca).
# 1. CaseSites: every case variant (multiallelic records split into biallelic ones, left-aligned when ref_fa is given)
#    as CHROM, POS, REF, ALT, and the positions joined into fetch regions (gaps <= site_merge_gap bp).
# 2. ControlStats, scattered over the control files (gs:// paths read in place with remote tabix, n_parallel batches of
#    regions_per_batch regions at a time; set requester_pays_project for requester-pays buckets): control records at
#    the case positions, split the same way, kept when CHROM:POS:REF:ALT matches a case variant, then AF, AC, AN,
#    N_HET, N_HOMREF and N_HOMALT computed from the genotypes for all control samples, the matched controls counted once,
#    and the matched controls counted once per case they were matched to (bcftools +fill-tags sample groups).
# 3. ConcatStats: the per-file counts joined into one tabix-indexed table.
# 4. AnnotateCase, scattered over the case files: the counts added to INFO as <control_name>_*, <control_name>_CTRL_*
#    and <control_name>_CTRL_DUP_*; case variants with no control record are left unannotated.
# 5. ConcatCase: the annotated case shards joined into one VCF.

workflow AnnotateCaseWithControlCounts {
    input {
        Array[String] case_vcfs
        Array[String] case_vcf_idxs
        Array[String] control_vcfs
        Array[String] control_vcf_idxs
        File control_weights
        String prefix

        String case_name = "CDP"
        String control_name = "AoU"
        Boolean control_pass_only = false
        Boolean control_gvcf = false
        String? requester_pays_project
        Int regions_per_batch = 500
        Int n_parallel = 16
        Int site_merge_gap = 1000

        File? ref_fa
        File? ref_fai

        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"

        RuntimeAttr? runtime_attr_case_sites
        RuntimeAttr? runtime_attr_control_stats
        RuntimeAttr? runtime_attr_concat_stats
        RuntimeAttr? runtime_attr_annotate_case
        RuntimeAttr? runtime_attr_concat_case
    }

    call MatchControlsByAncestry.CaseSites {
        input:
            case_vcfs = case_vcfs,
            case_vcf_idxs = case_vcf_idxs,
            site_merge_gap = site_merge_gap,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_case_sites
    }

    scatter (i in range(length(control_vcfs))) {
        call MatchControlsByAncestry.ControlStats {
            input:
                vcf = control_vcfs[i],
                vcf_idx = control_vcf_idxs[i],
                case_regions_bed = CaseSites.regions_bed,
                case_sites_tsv = CaseSites.sites_tsv,
                control_weights = control_weights,
                pass_only = control_pass_only,
                gvcf = control_gvcf,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                requester_pays_project = requester_pays_project,
                regions_per_batch = regions_per_batch,
                n_parallel = n_parallel,
                prefix = "~{prefix}.~{control_name}.shard_~{i}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_control_stats
        }
    }

    call MatchControlsByAncestry.ConcatStats {
        input:
            stats_tsvs = ControlStats.stats_tsv,
            control_name = control_name,
            prefix = prefix,
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_stats
    }

    scatter (i in range(length(case_vcfs))) {
        call MatchControlsByAncestry.AnnotateCase {
            input:
                vcf = case_vcfs[i],
                vcf_idx = case_vcf_idxs[i],
                stats_tsv = ConcatStats.stats_tsv,
                stats_tsv_idx = ConcatStats.stats_tsv_idx,
                control_name = control_name,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                prefix = "~{prefix}.~{case_name}.shard_~{i}",
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_annotate_case
        }
    }

    call MatchControlsByAncestry.ConcatCase {
        input:
            vcfs = AnnotateCase.annotated_vcf,
            vcf_idxs = AnnotateCase.annotated_vcf_idx,
            prefix = "~{prefix}.~{case_name}.~{control_name}_annotated",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_case
    }

    output {
        File annotated_case_vcf = ConcatCase.vcf
        File annotated_case_vcf_idx = ConcatCase.vcf_idx
        File control_stats_tsv = ConcatStats.stats_tsv
        File control_stats_tsv_idx = ConcatStats.stats_tsv_idx
        File case_sites_tsv = CaseSites.sites_tsv
    }
}
