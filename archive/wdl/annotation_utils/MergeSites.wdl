version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow MergeSites {
    meta {
        description: [
            "This utility merges redundant records at the site level within a VCF by collapsing near-identical deletions and insertions. Deletions are collapsed using size-, reciprocal-overlap, sequence- and sample-similarity thresholds plus a breakpoint distance, insertions using size-, sequence- and sample-similarity plus a breakpoint distance, while all other variants pass through untouched. It outputs the merged VCF."
        ]
    }

    parameter_meta {
        vcf: "VCF to merge."
        vcf_idx: "Index for VCF."
        del_sample_similarity: "Minimum sample similarity for collapsing deletions."
        ins_sample_similarity: "Minimum sample similarity for collapsing insertions."
        del_breakpoint_window: "Maximum breakpoint distance, in bp, for collapsing deletions."
        del_reciprocal_overlap: "Minimum reciprocal overlap for collapsing deletions."
        del_sequence_similarity: "Minimum sequence similarity for collapsing deletions."
        del_size_similarity: "Minimum size similarity for collapsing deletions."
        del_size_max: "Maximum deletion size to collapse, or `-1` for no maximum."
        del_size_min: "Minimum deletion size to collapse."
        ins_breakpoint_window: "Maximum breakpoint distance, in bp, for collapsing insertions."
        ins_reciprocal_overlap: "Minimum reciprocal overlap for collapsing insertions."
        ins_sequence_similarity: "Minimum sequence similarity for collapsing insertions."
        ins_size_similarity: "Minimum size similarity for collapsing insertions."
        ins_size_max: "Maximum insertion size to collapse, or `-1` for no maximum."
        ins_size_min: "Minimum insertion size to collapse."
        shard_bin_size: "If set, shards each contig into regions of roughly this many base pairs, run in parallel."
        ref_fai: "From references."
        merged_vcf: "Site-merged VCF."
        merged_vcf_idx: "Index for the merged VCF."
    }

    input {
        File vcf
        File vcf_idx
        String prefix

        Float del_sample_similarity = 0.5
        Float ins_sample_similarity = 0.5

        Int del_breakpoint_window = 500
        Float del_reciprocal_overlap = 0.0
        Float del_sequence_similarity = 0.5
        Float del_size_similarity = 0.5
        Int del_size_max = 50000
        Int del_size_min = 0

        Int ins_breakpoint_window = 200
        Float ins_reciprocal_overlap = 0.0
        Float ins_sequence_similarity = 0.5
        Float ins_size_similarity = 0.5
        Int ins_size_max = 50000
        Int ins_size_min = 0

        Int? shard_bin_size
        File? ref_fai

        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_to_region
        RuntimeAttr? runtime_attr_subset_by_type
        RuntimeAttr? runtime_attr_collapse_dels
        RuntimeAttr? runtime_attr_collapse_ins
        RuntimeAttr? runtime_attr_concat_shard
        RuntimeAttr? runtime_attr_concat
    }

    if (defined(shard_bin_size)) {
        call Helpers.CreateShardsFromVcfIndex {
            input:
                vcf_idx = vcf_idx,
                ref_fai = select_first([ref_fai]),
                shard_bin_size = select_first([shard_bin_size]),
                prefix = "~{prefix}.shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_shards
        }

        scatter (i in range(length(CreateShardsFromVcfIndex.shard_regions))) {
            String shard_region = CreateShardsFromVcfIndex.shard_regions[i]

            call Helpers.SubsetVcfToRegion {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    region = shard_region,
                    prefix = "~{prefix}.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_to_region
            }

            call Helpers.SubsetVcfByArgs as SubsetDelsShard {
                input:
                    vcf = SubsetVcfToRegion.subset_vcf,
                    vcf_idx = SubsetVcfToRegion.subset_vcf_idx,
                    include_args = 'INFO/allele_type~"del"',
                    prefix = "~{prefix}.shard_~{i}.dels",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_by_type
            }

            call Helpers.SubsetVcfByArgs as SubsetInsShard {
                input:
                    vcf = SubsetVcfToRegion.subset_vcf,
                    vcf_idx = SubsetVcfToRegion.subset_vcf_idx,
                    include_args = 'INFO/allele_type~"ins" || INFO/allele_type~"dup"',
                    prefix = "~{prefix}.shard_~{i}.ins",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_by_type
            }

            call Helpers.SubsetVcfByArgs as SubsetOtherShard {
                input:
                    vcf = SubsetVcfToRegion.subset_vcf,
                    vcf_idx = SubsetVcfToRegion.subset_vcf_idx,
                    include_args = 'INFO/allele_type!~"del" && INFO/allele_type!~"ins" && INFO/allele_type!~"dup"',
                    prefix = "~{prefix}.shard_~{i}.other",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_by_type
            }

            call Helpers.ConsolidateCollapsedSites as CollapseDelsShard {
                input:
                    vcf = SubsetDelsShard.subset_vcf,
                    vcf_idx = SubsetDelsShard.subset_vcf_idx,
                    breakpoint_window = del_breakpoint_window,
                    reciprocal_overlap = del_reciprocal_overlap,
                    sample_similarity = del_sample_similarity,
                    sequence_similarity = del_sequence_similarity,
                    size_similarity = del_size_similarity,
                    size_min = del_size_min,
                    size_max = del_size_max,

                    keep_strategy = "common",
                    set_merge_annotations = false,
                    strip_format_to_gt = false,
                    prefix = "~{prefix}.shard_~{i}.dels.collapsed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_collapse_dels
            }

            call Helpers.ConsolidateCollapsedSites as CollapseInsShard {
                input:
                    vcf = SubsetInsShard.subset_vcf,
                    vcf_idx = SubsetInsShard.subset_vcf_idx,
                    breakpoint_window = ins_breakpoint_window,
                    reciprocal_overlap = ins_reciprocal_overlap,
                    sample_similarity = ins_sample_similarity,
                    sequence_similarity = ins_sequence_similarity,
                    size_similarity = ins_size_similarity,
                    size_min = ins_size_min,
                    size_max = ins_size_max,

                    keep_strategy = "common",
                    set_merge_annotations = false,
                    strip_format_to_gt = false,
                    prefix = "~{prefix}.shard_~{i}.ins.collapsed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_collapse_ins
            }

            call Helpers.ConcatVcfs as ConcatShard {
                input:
                    vcfs = [CollapseDelsShard.consolidated_vcf, CollapseInsShard.consolidated_vcf, SubsetOtherShard.subset_vcf],
                    vcf_idxs = [CollapseDelsShard.consolidated_vcf_idx, CollapseInsShard.consolidated_vcf_idx, SubsetOtherShard.subset_vcf_idx],
                    allow_overlaps = true,
                    naive = false,
                    prefix = "~{prefix}.shard_~{i}.merged",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shard
            }
        }

        call Helpers.ConcatVcfs as ConcatShards {
            input:
                vcfs = ConcatShard.concat_vcf,
                vcf_idxs = ConcatShard.concat_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    if (!defined(shard_bin_size)) {
        # Split VCF into DEL, INS and OTHER
        call Helpers.SubsetVcfByArgs as SubsetDels {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                include_args = 'INFO/allele_type~"del"',
                prefix = "~{prefix}.dels",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_by_type
        }

        call Helpers.SubsetVcfByArgs as SubsetIns {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                include_args = 'INFO/allele_type~"ins" || INFO/allele_type~"dup"',
                prefix = "~{prefix}.ins",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_by_type
        }

        call Helpers.SubsetVcfByArgs as SubsetOther {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                include_args = 'INFO/allele_type!~"del" && INFO/allele_type!~"ins" && INFO/allele_type!~"dup"',
                prefix = "~{prefix}.other",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset_by_type
        }

        # Collapse DELs
        call Helpers.ConsolidateCollapsedSites as CollapseDels {
            input:
                vcf = SubsetDels.subset_vcf,
                vcf_idx = SubsetDels.subset_vcf_idx,
                breakpoint_window =del_breakpoint_window,
                reciprocal_overlap = del_reciprocal_overlap,
                sample_similarity = del_sample_similarity,
                sequence_similarity = del_sequence_similarity,
                size_similarity = del_size_similarity,
                size_min = del_size_min,
                size_max = del_size_max,
                keep_strategy = "common",
                set_merge_annotations = false,
                strip_format_to_gt = false,
                prefix = "~{prefix}.dels.collapsed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_collapse_dels
        }

        # Collapse INS
        call Helpers.ConsolidateCollapsedSites as CollapseIns {
            input:
                vcf = SubsetIns.subset_vcf,
                vcf_idx = SubsetIns.subset_vcf_idx,
                breakpoint_window =ins_breakpoint_window,
                reciprocal_overlap = ins_reciprocal_overlap,
                sample_similarity = ins_sample_similarity,
                sequence_similarity = ins_sequence_similarity,
                size_similarity = ins_size_similarity,
                size_min = ins_size_min,
                size_max = ins_size_max,
                keep_strategy = "common",
                set_merge_annotations = false,
                strip_format_to_gt = false,
                prefix = "~{prefix}.ins.collapsed",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_collapse_ins
        }

        call Helpers.ConcatVcfs {
            input:
                vcfs = [CollapseDels.consolidated_vcf, CollapseIns.consolidated_vcf, SubsetOther.subset_vcf],
                vcf_idxs = [CollapseDels.consolidated_vcf_idx, CollapseIns.consolidated_vcf_idx, SubsetOther.subset_vcf_idx],
                allow_overlaps = true,
                naive = false,
                prefix = "~{prefix}.merged",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File merged_vcf = select_first([ConcatShards.concat_vcf, ConcatVcfs.concat_vcf])
        File merged_vcf_idx = select_first([ConcatShards.concat_vcf_idx, ConcatVcfs.concat_vcf_idx])
    }
}
