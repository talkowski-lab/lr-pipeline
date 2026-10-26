version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PostprocessCallset {
    meta {
        description: [
            "This utility bundles every genotype-update and post-processing step applied to a near-final callset into one workflow, with a required `run_` Boolean guarding each step so that the input VCF is left untouched when all are set to `false`. The per-record steps are applied in a single pass over the VCF: each variant is first matched against `transfer_vcf` and has its genotypes transferred (when `run_transfer_genotypes` is set) using its unmodified properties, after which the remaining steps - unphasing, ploidy normalization, allele-count recomputation, TR-ID decrementing, MEI pruning, homopolymer flagging, singleton filtering and same-coordinate sorting - run in order. Some steps require an accompanying field - `run_transfer_genotypes` needs `transfer_vcf`, `run_unphase_samples` needs `unphase_samples`, and `run_normalize_ploidy` needs `ped`. Setting `run_compute_allele_counts` writes `INFO/AC` and `INFO/AN` from the genotypes as they stand after ploidy normalization, which the two singleton filters read, so they can run before `AnnotateAF`. The per-record pass can optionally be region-sharded via `shard_bin_size`. Two contig-wide steps then follow on the post-processed genotypes, in order. `run_resolve_haplotype_overlaps` scatters per sample over the sample's non-ref calls (excluding `allele_type='trv'` and `INFO/TR_ENVELOPED` variants), sweeps each haplotype for overlapping pairs and keeps the variant spanning more reference sequence (larger `len(REF)`), breaking ties by GQ, then `INFO/allele_length`, then type rank (DEL > INS > SNV), then QUAL, then input-file order; the loser's `GT`, `GQ`, `EV`, `BEV`, `AD` and `PL` are cleared and every pair is reported in `overlap_tsv`. `run_normalize_allele_types` then reclassifies `allele_type` and records the original type in a new `allele_subtype` field: `allele_type=dup` variants all get `allele_subtype=tandem_dup` and keep `allele_type=dup` only when tandem against their `INFO/ORIGIN` source - the insertion and ORIGIN lengths must meet `dup_size_similarity` and POS must fall within ORIGIN or within `dup_breakpoint_window` bases of its breakpoints - and otherwise become `ins`; `complex_dup`, `dup_interspersed`, `inv_dup`, `alu_ins`, `line_ins`, `sva_ins` and `numt` become `ins`, and `alu_del`, `line_del` and `sva_del` become `del`. REF/ALT/POS are never modified, and other `allele_type` values pass through unchanged."
        ]
    }

    parameter_meta {
        vcf: "VCF to post-process."
        vcf_idx: "Index for VCF to post-process."
        contigs: "Contigs to process within the input VCF."
        shard_bin_size: "Region-bin size, in bp, used when sharding the per-record pass."
        run_clean_vcf_header: "Whether to run the header-cleaning step."
        run_compute_allele_counts: "Whether to recompute `INFO/AC` and `INFO/AN` from the genotypes after ploidy normalization, adding their header definitions when absent; the singleton filters read `AC`, so set this when `AnnotateAF` has not yet run."
        run_decrement_trv_ids: "Whether to decrement tandem-repeat variant IDs."
        run_drop_filters: "Whether to drop the FILTER values listed in `drop_filters`."
        run_filter_assembly_only_singletons: "Whether to apply the `ASSEMBLY_ONLY_SINGLETON` filter and emit a matching TSV."
        run_filter_single_read_singletons: "Whether to apply the `SINGLE_READ_SUPPORT` filter to singleton calls."
        run_flag_homopolymer_trvs: "Whether to flag tandem repeats with a length-1 shortest motif as `HOMOPOLYMER_TRV`."
        run_normalize_allele_types: "Whether to reclassify `allele_type` values, recording the original type in `allele_subtype`; runs last."
        run_normalize_ploidy: "Whether to normalize ploidy by sex - clearing chrY female calls, making chrX/chrY male calls hemizygous, enforcing diploidy and right-aligning unphased calls (requires `ped`)."
        run_prune_meis: "Whether to reclassify mobile elements whose length falls outside the expected bounds back to plain insertions/deletions."
        run_reassign_suffixes: "Whether to run the variant-ID suffix reassignment step."
        run_resolve_haplotype_overlaps: "Whether to clear the losing call of each pair of non-TR variants overlapping on the same haplotype, after the per-record pass."
        run_sorting: "Whether to sort records sharing a coordinate by absolute allele length and variant ID."
        run_transfer_genotypes: "Whether to transfer genotypes from `transfer_vcf` onto heterozygous calls (run first; requires `transfer_vcf`)."
        run_unphase_samples: "Whether to unphase the samples in `unphase_samples` (requires `unphase_samples`)."
        unphase_samples: "Samples to unphase when `run_unphase_samples` is set (defaults to empty)."
        drop_filters: "FILTER values to drop when `run_drop_filters` is set."
        transfer_vcf: "VCF whose genotypes are transferred when `run_transfer_genotypes` is set."
        transfer_vcf_idx: "Index for `transfer_vcf`."
        ped: "Cohort pedigree file, used for ploidy normalization when `run_normalize_ploidy` is set."
        dup_breakpoint_window: "Maximum distance (bp) between insertion POS and ORIGIN breakpoints to pass the breakpoint check."
        dup_size_similarity: "Minimum size similarity ratio (relative to the larger of the two lengths) between insertion and ORIGIN lengths."
        min_dup_size: "Minimum insertion size (bp) to consider for the tandem check."
        post_processed_vcf: "Post-processed VCF."
        post_processed_vcf_idx: "Index for the post-processed VCF."
        assembly_only_singletons_tsv: "Optional TSV containing one row per assembly-only singleton ALT allele; present only when `run_filter_assembly_only_singletons` is true."
        overlap_tsv: "Optional TSV of all detected overlap pairs, with columns `sample`, `haplotype`, `variant_id_retained`, `var_type_retained`, `size_bin_retained`, `variant_id_cleared`, `var_type_cleared`, `size_bin_cleared`; present only when `run_resolve_haplotype_overlaps` is true."
    }

    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? shard_bin_size

        Boolean run_clean_vcf_header
        Boolean run_compute_allele_counts
        Boolean run_decrement_trv_ids
        Boolean run_drop_filters
        Boolean run_filter_assembly_only_singletons
        Boolean run_filter_single_read_singletons
        Boolean run_flag_homopolymer_trvs
        Boolean run_normalize_allele_types
        Boolean run_normalize_ploidy
        Boolean run_prune_meis
        Boolean run_reassign_suffixes
        Boolean run_resolve_haplotype_overlaps
        Boolean run_sorting
        Boolean run_transfer_genotypes
        Boolean run_unphase_samples

        Array[String] unphase_samples = []
        Array[String] drop_filters = []
        File? transfer_vcf
        File? transfer_vcf_idx
        File? ped

        Int dup_breakpoint_window = 10
        Float dup_size_similarity = 0.9
        Int min_dup_size = 50

        String utils_docker

        RuntimeAttr? runtime_attr_subset_base
        RuntimeAttr? runtime_attr_subset_transfer
        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_post_process
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_get_samples
        RuntimeAttr? runtime_attr_extract_sample
        RuntimeAttr? runtime_attr_detect_overlaps
        RuntimeAttr? runtime_attr_concat_tsvs
        RuntimeAttr? runtime_attr_resolve_and_normalize
        RuntimeAttr? runtime_attr_concat
    }

    # Cromwell rejects an undefined optional entering a nested scatter, so these cross as one-element arrays
    Array[File?] ped_wrapped = [ped]
    Array[RuntimeAttr?] runtime_attr_subset_base_wrapped = [runtime_attr_subset_base]
    Array[RuntimeAttr?] runtime_attr_subset_transfer_wrapped = [runtime_attr_subset_transfer]
    Array[RuntimeAttr?] runtime_attr_post_process_wrapped = [runtime_attr_post_process]
    Array[RuntimeAttr?] runtime_attr_extract_sample_wrapped = [runtime_attr_extract_sample]
    Array[RuntimeAttr?] runtime_attr_detect_overlaps_wrapped = [runtime_attr_detect_overlaps]
    Array[RuntimeAttr?] runtime_attr_resolve_and_normalize_wrapped = [runtime_attr_resolve_and_normalize]

    Boolean single_contig = length(contigs) == 1
    Boolean do_transfer_genotypes = run_transfer_genotypes && defined(transfer_vcf)
    Boolean do_unphase_samples = run_unphase_samples && length(unphase_samples) > 0
    Boolean do_normalize_ploidy = run_normalize_ploidy && defined(ped)
    Boolean do_resolve_or_normalize = run_resolve_haplotype_overlaps || run_normalize_allele_types

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetBase {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.base",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_base
            }
        }

        File contig_base_vcf = select_first([SubsetBase.subset_vcf, vcf])
        File contig_base_vcf_idx = select_first([SubsetBase.subset_vcf_idx, vcf_idx])

        if (do_transfer_genotypes) {
            call Helpers.SubsetVcfToContig as SubsetTransfer {
                input:
                    vcf = select_first([transfer_vcf]),
                    vcf_idx = select_first([transfer_vcf_idx]),
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.transfer",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_transfer
            }
        }

        File? transfer_source_vcf = SubsetTransfer.subset_vcf
        File? transfer_source_vcf_idx = SubsetTransfer.subset_vcf_idx

        if (defined(shard_bin_size)) {
            call Helpers.CreateContigShards {
                input:
                    vcfs = select_all([contig_base_vcf, transfer_source_vcf]),
                    vcf_idxs = select_all([contig_base_vcf_idx, transfer_source_vcf_idx]),
                    contig = contig,
                    shard_bin_size = select_first([shard_bin_size]),
                    prefix = "~{prefix}.~{contig}.shards",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_create_shards
            }

            # Declare as non-optional File so the value can cross into the shard scatter
            File shard_transfer_vcf = select_first([transfer_source_vcf, contig_base_vcf])
            File shard_transfer_vcf_idx = select_first([transfer_source_vcf_idx, contig_base_vcf_idx])

            scatter (i in range(length(CreateContigShards.shard_regions))) {
                String shard_region = CreateContigShards.shard_regions[i]

                call Helpers.SubsetVcfToRegion as SubsetBaseShard {
                    input:
                        vcf = contig_base_vcf,
                        vcf_idx = contig_base_vcf_idx,
                        region = shard_region,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.base",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_base_wrapped[0]
                }

                if (do_transfer_genotypes) {
                    call Helpers.SubsetVcfToRegion as SubsetTransferShard {
                        input:
                            vcf = shard_transfer_vcf,
                            vcf_idx = shard_transfer_vcf_idx,
                            region = shard_region,
                            prefix = "~{prefix}.~{contig}.shard_~{i}.transfer",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_subset_transfer_wrapped[0]
                    }
                }

                call PostProcessTask as PostProcessShard {
                    input:
                        base_vcf = SubsetBaseShard.subset_vcf,
                        base_vcf_idx = SubsetBaseShard.subset_vcf_idx,
                        transfer_vcf = SubsetTransferShard.subset_vcf,
                        transfer_vcf_idx = SubsetTransferShard.subset_vcf_idx,
                        ped = ped_wrapped[0],
                        unphase_samples = unphase_samples,
                        transfer_genotypes = do_transfer_genotypes,
                        unphase = do_unphase_samples,
                        normalize_ploidy = do_normalize_ploidy,
                        clean_vcf_header = run_clean_vcf_header,
                        compute_allele_counts = run_compute_allele_counts,
                        decrement_trv_ids = run_decrement_trv_ids,
                        drop = run_drop_filters,
                        drop_filters = drop_filters,
                        filter_assembly_only_singletons = run_filter_assembly_only_singletons,
                        filter_singletons = run_filter_single_read_singletons,
                        flag_homopolymer_trvs = run_flag_homopolymer_trvs,
                        prune_meis = run_prune_meis,
                        reassign_suffixes = run_reassign_suffixes,
                        sort_records = run_sorting,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.post_processed",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_post_process_wrapped[0]
                }

                File shard_final_vcf = PostProcessShard.processed_vcf
                File shard_final_vcf_idx = PostProcessShard.processed_vcf_idx
            }

            call Helpers.ConcatVcfs as ConcatShards {
                input:
                    vcfs = shard_final_vcf,
                    vcf_idxs = shard_final_vcf_idx,
                    allow_overlaps = false,
                    naive = true,
                    prefix = "~{prefix}.~{contig}.post_processed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }

            if (run_filter_assembly_only_singletons) {
                call Helpers.ConcatTsvs as ConcatShardAssemblyOnlySingletonTsvs {
                    input:
                        tsvs = PostProcessShard.assembly_only_singletons_tsv,
                        sort_output = false,
                        preserve_header = true,
                        prefix = "~{prefix}.~{contig}.assembly_only_singletons",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_concat_shards
                }
            }
        }

        if (!defined(shard_bin_size)) {
            call PostProcessTask {
                input:
                    base_vcf = contig_base_vcf,
                    base_vcf_idx = contig_base_vcf_idx,
                    transfer_vcf = transfer_source_vcf,
                    transfer_vcf_idx = transfer_source_vcf_idx,
                    ped = ped,
                    unphase_samples = unphase_samples,
                    transfer_genotypes = do_transfer_genotypes,
                    unphase = do_unphase_samples,
                    normalize_ploidy = do_normalize_ploidy,
                    clean_vcf_header = run_clean_vcf_header,
                    compute_allele_counts = run_compute_allele_counts,
                    decrement_trv_ids = run_decrement_trv_ids,
                    drop = run_drop_filters,
                    drop_filters = drop_filters,
                    filter_assembly_only_singletons = run_filter_assembly_only_singletons,
                    filter_singletons = run_filter_single_read_singletons,
                    flag_homopolymer_trvs = run_flag_homopolymer_trvs,
                    prune_meis = run_prune_meis,
                    reassign_suffixes = run_reassign_suffixes,
                    sort_records = run_sorting,
                    prefix = "~{prefix}.~{contig}.post_processed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_post_process
            }
        }

        File contig_processed_vcf = select_first([ConcatShards.concat_vcf, PostProcessTask.processed_vcf])
        File contig_processed_vcf_idx = select_first([ConcatShards.concat_vcf_idx, PostProcessTask.processed_vcf_idx])
        File? contig_assembly_only_singletons_tsv = if run_filter_assembly_only_singletons then select_first([ConcatShardAssemblyOnlySingletonTsvs.concatenated_tsv, PostProcessTask.assembly_only_singletons_tsv]) else PostProcessTask.assembly_only_singletons_tsv

        # Overlaps are detected on the whole contig, so a pair split across region shards is still found
        if (run_resolve_haplotype_overlaps) {
            call Helpers.GetSamplesFromVcf {
                input:
                    vcf = contig_processed_vcf,
                    vcf_idx = contig_processed_vcf_idx,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_get_samples
            }

            scatter (sample in GetSamplesFromVcf.samples) {
                call Helpers.ExtractSample {
                    input:
                        vcf = contig_processed_vcf,
                        vcf_idx = contig_processed_vcf_idx,
                        sample = sample,
                        extra_args = "-e 'INFO/allele_type=\"trv\" || INFO/TR_ENVELOPED=1'",
                        prefix = "~{prefix}.~{contig}.~{sample}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_extract_sample_wrapped[0]
                }

                call DetectSampleOverlaps {
                    input:
                        vcf = ExtractSample.subset_vcf,
                        vcf_idx = ExtractSample.subset_vcf_idx,
                        sample = sample,
                        prefix = "~{prefix}.~{contig}.~{sample}",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_detect_overlaps_wrapped[0]
                }
            }

            call Helpers.ConcatTsvs as ConcatContigOverlapTsvs {
                input:
                    tsvs = DetectSampleOverlaps.overlaps_tsv,
                    sort_output = false,
                    preserve_header = true,
                    prefix = "~{prefix}.~{contig}.overlaps",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_tsvs
            }

            call Helpers.ConcatTsvs as ConcatClearedTsvs {
                input:
                    tsvs = DetectSampleOverlaps.cleared_calls_tsv,
                    sort_output = false,
                    preserve_header = false,
                    prefix = "~{prefix}.~{contig}.cleared_calls",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_tsvs
            }
        }

        if (do_resolve_or_normalize) {
            # Reuse the region shards of the per-record pass, which already hold the post-processed records
            if (defined(shard_bin_size)) {
                Array[File] resolve_shard_vcfs = select_first([shard_final_vcf])
                Array[File] resolve_shard_vcf_idxs = select_first([shard_final_vcf_idx])
                Array[File?] cleared_calls_tsv_wrapped = [ConcatClearedTsvs.concatenated_tsv]

                scatter (i in range(length(resolve_shard_vcfs))) {
                    call ResolveAndNormalizeTask as ResolveAndNormalizeShard {
                        input:
                            vcf = resolve_shard_vcfs[i],
                            vcf_idx = resolve_shard_vcf_idxs[i],
                            cleared_calls_tsv = cleared_calls_tsv_wrapped[0],
                            resolve_overlaps = run_resolve_haplotype_overlaps,
                            normalize_allele_types = run_normalize_allele_types,
                            dup_breakpoint_window = dup_breakpoint_window,
                            dup_size_similarity = dup_size_similarity,
                            min_dup_size = min_dup_size,
                            prefix = "~{prefix}.~{contig}.shard_~{i}.resolved",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_resolve_and_normalize_wrapped[0]
                    }
                }

                call Helpers.ConcatVcfs as ConcatResolvedShards {
                    input:
                        vcfs = ResolveAndNormalizeShard.output_vcf,
                        vcf_idxs = ResolveAndNormalizeShard.output_vcf_idx,
                        allow_overlaps = false,
                        naive = true,
                        prefix = "~{prefix}.~{contig}.resolved",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_concat_shards
                }
            }

            if (!defined(shard_bin_size)) {
                call ResolveAndNormalizeTask {
                    input:
                        vcf = contig_processed_vcf,
                        vcf_idx = contig_processed_vcf_idx,
                        cleared_calls_tsv = ConcatClearedTsvs.concatenated_tsv,
                        resolve_overlaps = run_resolve_haplotype_overlaps,
                        normalize_allele_types = run_normalize_allele_types,
                        dup_breakpoint_window = dup_breakpoint_window,
                        dup_size_similarity = dup_size_similarity,
                        min_dup_size = min_dup_size,
                        prefix = "~{prefix}.~{contig}.resolved",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_resolve_and_normalize
                }
            }
        }

        File contig_final_vcf = select_first([ConcatResolvedShards.concat_vcf, ResolveAndNormalizeTask.output_vcf, contig_processed_vcf])
        File contig_final_vcf_idx = select_first([ConcatResolvedShards.concat_vcf_idx, ResolveAndNormalizeTask.output_vcf_idx, contig_processed_vcf_idx])
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = contig_final_vcf,
            vcf_idxs = contig_final_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    if (run_filter_assembly_only_singletons) {
        call Helpers.ConcatTsvs as ConcatAssemblyOnlySingletonTsvs {
            input:
                tsvs = select_all(contig_assembly_only_singletons_tsv),
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.assembly_only_singletons",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    if (run_resolve_haplotype_overlaps) {
        call Helpers.ConcatTsvs as ConcatOverlapTsvs {
            input:
                tsvs = select_all(ConcatContigOverlapTsvs.concatenated_tsv),
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.overlaps",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    output {
        File post_processed_vcf = ConcatVcfs.concat_vcf
        File post_processed_vcf_idx = ConcatVcfs.concat_vcf_idx
        File? assembly_only_singletons_tsv = ConcatAssemblyOnlySingletonTsvs.concatenated_tsv
        File? overlap_tsv = ConcatOverlapTsvs.concatenated_tsv
    }
}

task PostProcessTask {
    input {
        File base_vcf
        File base_vcf_idx
        File? transfer_vcf
        File? transfer_vcf_idx
        File? ped
        Array[String] unphase_samples = []
        Boolean transfer_genotypes
        Boolean unphase
        Boolean normalize_ploidy
        Boolean compute_allele_counts
        Boolean decrement_trv_ids
        Boolean prune_meis
        Boolean flag_homopolymer_trvs
        Boolean sort_records
        Boolean filter_assembly_only_singletons
        Boolean filter_singletons
        Boolean drop
        Array[String] drop_filters = []
        Boolean clean_vcf_header
        Boolean reassign_suffixes
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import re
from collections import defaultdict

import pysam

transfer_genotypes = ~{true="True" false="False" transfer_genotypes}
unphase = ~{true="True" false="False" unphase}
normalize_ploidy = ~{true="True" false="False" normalize_ploidy}
compute_allele_counts = ~{true="True" false="False" compute_allele_counts}
decrement_trv_ids = ~{true="True" false="False" decrement_trv_ids}
prune_meis = ~{true="True" false="False" prune_meis}
flag_homopolymer_trvs = ~{true="True" false="False" flag_homopolymer_trvs}
sort_records = ~{true="True" false="False" sort_records}
filter_assembly_only_singletons = ~{true="True" false="False" filter_assembly_only_singletons}
filter_singletons = ~{true="True" false="False" filter_singletons}
drop = ~{true="True" false="False" drop}
clean_vcf_header = ~{true="True" false="False" clean_vcf_header}
reassign_suffixes = ~{true="True" false="False" reassign_suffixes}

unphase_list = ["~{sep='\", \"' unphase_samples}"] if ~{length(unphase_samples)} > 0 else []
unphase_samples_set = set(unphase_list)
drop_filters_list = ["~{sep='\", \"' drop_filters}"] if ~{length(drop_filters)} > 0 else []
drop_filters_set = set(drop_filters_list)

ASSEMBLY_CALLERS = {"dipcall", "hapdiff"}
CALLER_RE = re.compile(r"^([A-Za-z0-9]+)")

def parse_ped(path):
    sex_by_sample = {}
    with open(path, "r") as handle:
        for line in handle:
            fields = line.strip().split()
            if not fields:
                continue
            sample_id = fields[1]
            sex_code = fields[4]
            if sex_code == "1":
                sex_by_sample[sample_id] = "M"
            elif sex_code == "2":
                sex_by_sample[sample_id] = "F"
            else:
                sex_by_sample[sample_id] = None
    return sex_by_sample

def get_scalar(value):
    if isinstance(value, (list, tuple)):
        for item in value:
            if item is not None:
                return item
        return None
    return value

def is_heterozygous(gt):
    if gt is None:
        return False
    called = [allele for allele in gt if allele is not None]
    return len(called) >= 2 and len(set(called)) > 1

def clear_format_fields(sample_data):
    sample_data["GT"] = (None, None)
    sample_data.phased = False

def right_align_unphased(gt):
    if gt is None:
        return gt
    return tuple(sorted(gt, key=lambda allele: (allele is not None, allele if allele is not None else -1)))

def compute_record_allele_counts(record):
    counts = [0] * len(record.alts)
    allele_number = 0
    for sample_data in record.samples.values():
        gt = sample_data.get("GT")
        if gt is None:
            continue
        for allele in gt:
            if allele is None:
                continue
            allele_number += 1
            if allele > 0:
                counts[allele - 1] += 1
    record.info["AC"] = tuple(counts)
    record.info["AN"] = allele_number

def decrement_trv_id(trv_id):
    if not trv_id:
        return trv_id
    head, _, tail = trv_id.rpartition("-")
    if not head.endswith("-TRV"):
        return trv_id
    try:
        ref_len = int(tail)
    except ValueError:
        return trv_id
    return "{}-{}".format(head, ref_len - 1)

def make_male_hemizygous(gt, phased):
    if gt is None:
        return gt

    alleles = list(gt)
    called_positions = [index for index, allele in enumerate(alleles) if allele is not None]
    if len(called_positions) <= 1:
        return tuple(alleles)

    alt_positions = [index for index, allele in enumerate(alleles) if allele is not None and allele > 0]

    if phased:
        if len(alt_positions) == 1:
            keep_index = alt_positions[0]
        elif alt_positions:
            keep_index = alt_positions[-1]
        else:
            keep_index = called_positions[-1]

        new_gt = [None] * len(alleles)
        new_gt[keep_index] = alleles[keep_index]
        return tuple(new_gt)

    if alt_positions:
        keep_allele = alleles[alt_positions[-1]]
    else:
        keep_allele = alleles[called_positions[-1]]

    new_gt = [None] * len(alleles)
    new_gt[-1] = keep_allele
    return right_align_unphased(tuple(new_gt))

def prune_record_meis(record):
    allele_type = get_scalar(record.info.get("allele_type"))
    allele_length = get_scalar(record.info.get("allele_length"))
    if allele_type is None or allele_length is None:
        return

    length = abs(int(allele_length))
    if allele_type in {"alu_ins", "alu_del"} and (length < 250 or length > 350):
        record.info["allele_type"] = "ins" if "ins" in allele_type else "del"
        if "SUB_FAMILY" in record.info:
            del record.info["SUB_FAMILY"]
    elif allele_type in {"sva_ins", "sva_del"} and (length < 1000 or length > 4000):
        record.info["allele_type"] = "ins" if "ins" in allele_type else "del"
        if "SUB_FAMILY" in record.info:
            del record.info["SUB_FAMILY"]

def shortest_motif_length(record):
    motifs = record.info.get("MOTIFS")
    if motifs is None:
        return None

    motif_values = []
    raw_values = motifs if isinstance(motifs, (list, tuple)) else [motifs]
    for value in raw_values:
        if value is None:
            continue
        motif_values.extend(part for part in str(value).split(",") if part and part != ".")

    if not motif_values:
        return None
    return min(len(motif) for motif in motif_values)

def has_single_read_support(record):
    if 'AC' in record.info and len(record.alts) == 1 and record.info['AC'][0] > 2:
        return False

    alt_depths = []
    for sample_data in record.samples.values():
        ad = sample_data.get("AD")
        if ad is None or len(ad) < 2:
            continue

        alt_depth = 0
        has_alt_depth = False
        for value in ad[1:]:
            if value is not None:
                alt_depth += value
                has_alt_depth = True

        if has_alt_depth and alt_depth > 0:
            alt_depths.append(alt_depth)

    return len(alt_depths) == 1 and alt_depths[0] == 1

def parse_callers(ev_value):
    if ev_value is None:
        return None

    values = ev_value if isinstance(ev_value, (list, tuple)) else [ev_value]
    parts = []
    for value in values:
        if value is not None:
            parts.extend(str(value).split(","))

    callers = set()
    for part in parts:
        match = CALLER_RE.match(part.strip())
        if match:
            callers.add(match.group(1).lower())
    return callers or None

def assembly_only_singleton_alt_idx(record):
    allele_counts = record.info.get("AC")
    if allele_counts is None:
        return []
    if not isinstance(allele_counts, (list, tuple)):
        allele_counts = (allele_counts,)

    matches = []
    for alt_index, allele_count in enumerate(allele_counts):
        if allele_count != 1:
            continue

        allele_number = alt_index + 1
        carrier = None
        for sample in record.samples.values():
            genotype = sample.get("GT")
            if genotype is not None and allele_number in genotype:
                carrier = sample
                break

        if carrier is None:
            continue
        callers = parse_callers(carrier.get("EV"))
        if callers and callers.issubset(ASSEMBLY_CALLERS):
            matches.append(alt_index)
    return matches

def info_string(value):
    if value is None:
        return "."
    if isinstance(value, (list, tuple)):
        values = [str(item) for item in value if item is not None]
        return ",".join(values) if values else "."
    return str(value)

def write_assembly_only_singleton_rows(record, alt_idx, tsv_writer):
    if not alt_idx:
        return
    filters = ",".join(record.filter.keys()) if record.filter.keys() else "."
    for alt_index in alt_idx:
        variant_id = record.id or "{}-{}-{}-{}".format(
            record.chrom, record.pos, record.ref, record.alts[alt_index]
        )
        tsv_writer.write(
            "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                variant_id,
                info_string(record.info.get("allele_type")),
                info_string(record.info.get("allele_length")),
                info_string(record.info.get("dbSNP_ID")),
                info_string(record.info.get("gnomAD_V4_match_ID")),
                filters,
            )
        )

def regroup_header(header):
    group_order = ["fileformat", "OTHER", "contig", "INFO", "FILTER", "FORMAT", "ALT", "SAMPLE", "PEDIGREE"]
    groups = {g: [] for g in group_order}
    key_re = re.compile(r'^##([A-Za-z_]+)=')

    for record in header.records:
        line = str(record).rstrip("\n")
        if line.startswith("##bcftools_"):
            continue
        m = key_re.match(line)
        key = m.group(1) if m else None
        groups[key if key in groups else "OTHER"].append(line)

    new_header = pysam.VariantHeader()
    for group in group_order:
        for line in groups[group]:
            new_header.add_line(line)
    for sample in header.samples:
        new_header.add_sample(sample)
    return new_header

def id_base_and_suffix(rec_id):
    id_val = rec_id if rec_id else ""
    parts = id_val.rsplit('_', 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0], int(parts[1])
    return id_val, 0

def reassign_id_suffixes(buf):
    groups = defaultdict(list)
    for rec in buf:
        base, _ = id_base_and_suffix(rec.id)
        groups[base].append(rec)

    for base, recs in groups.items():
        if len(recs) == 1:
            recs[0].id = base
            continue
        recs.sort(key=lambda r: id_base_and_suffix(r.id)[1])
        for i, r in enumerate(recs, start=1):
            r.id = "{}_{}".format(base, i)

def flush_buffer(buf, out_vcf, tsv_writer, assembly_alt_idx, reassign_suffixes, sort_records):
    def custom_sort_key(rec):
        al_val = get_scalar(rec.info.get("allele_length"))
        abs_al = abs(int(al_val)) if al_val is not None else 0
        return (abs_al, id_base_and_suffix(rec.id))

    if reassign_suffixes:
        reassign_id_suffixes(buf)
    if sort_records:
        buf.sort(key=custom_sort_key)
    for r in buf:
        out_vcf.write(r)
        write_assembly_only_singleton_rows(r, assembly_alt_idx.pop(id(r), []), tsv_writer)

sex_by_sample = parse_ped("~{default="NONE" ped}") if normalize_ploidy else {}

base_reader = pysam.VariantFile("~{base_vcf}")
transfer_reader = pysam.VariantFile("~{transfer_vcf}") if transfer_genotypes else None

header = base_reader.header.copy()
if flag_homopolymer_trvs and "HOMOPOLYMER_TRV" not in header.info:
    header.info.add("HOMOPOLYMER_TRV", 0, "Flag", "Tandem repeat call where the shortest motif has length 1.")
if compute_allele_counts and "AC" not in header.info:
    header.info.add("AC", "A", "Integer", "Number of alleles observed.")
if compute_allele_counts and "AN" not in header.info:
    header.info.add("AN", 1, "Integer", "Total number of alleles genotyped.")
if filter_singletons and "SINGLE_READ_SUPPORT" not in header.filters:
    header.filters.add("SINGLE_READ_SUPPORT", None, None, "Variant supported by a single read in a single sample.")
if drop:
    for drop_filter in drop_filters_set:
        if drop_filter in header.filters:
            header.filters.remove_header(drop_filter)
if filter_assembly_only_singletons and "ASSEMBLY_ONLY_SINGLETON" not in header.filters:
    header.filters.add(
        "ASSEMBLY_ONLY_SINGLETON",
        None,
        None,
        "Singleton variant whose carrier is supported only by assembly-based callers.",
    )
if clean_vcf_header:
    header = regroup_header(header)

output_writer = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)
tsv_writer = open("~{prefix}.assembly_only_singletons.tsv", "w")
tsv_writer.write("ID\tallele_type\tallele_length\tdbSNP_ID\tgnomAD_V4_match_ID\tFILTER\n")
base_samples = list(base_reader.header.samples)
shared_samples = set(base_samples)
if transfer_genotypes:
    shared_samples &= set(transfer_reader.header.samples)

buffer = []
assembly_alt_idx = {}
current_chrom = None
current_pos = None

for record in base_reader:
    if drop and drop_filters_set.intersection(record.filter.keys()):
        continue

    record.translate(output_writer.header)

    # Transfer genotypes first, matching on the unmodified variant properties
    if transfer_genotypes:
        match = None
        for candidate in transfer_reader.fetch(record.chrom, record.start, record.stop):
            if candidate.chrom == record.chrom and candidate.pos == record.pos and candidate.id == record.id and candidate.ref == record.ref and candidate.alts == record.alts:
                match = candidate
                break
        if match is not None:
            for sample in shared_samples:
                base_gt = record.samples[sample].get("GT")
                if is_heterozygous(base_gt):
                    record.samples[sample]["GT"] = match.samples[sample].get("GT")
                    record.samples[sample].phased = match.samples[sample].phased
        else:
            for sample in base_samples:
                record.samples[sample].phased = False

    if unphase or normalize_ploidy:
        for sample in base_samples:
            sample_data = record.samples[sample]

            # Unphase listed samples
            if unphase and sample in unphase_samples_set:
                if sample_data.get("GT") is not None:
                    sample_data.phased = False

            # Normalize ploidy using sample sex
            if normalize_ploidy:
                sample_sex = sex_by_sample.get(sample)

                # Clear FORMAT fields for females on chrY
                if record.chrom == "chrY" and sample_sex == "F":
                    clear_format_fields(sample_data)
                    continue

                # Make male calls hemizygous on chrX and chrY
                if record.chrom in {"chrX", "chrY"} and sample_sex == "M":
                    sample_data["GT"] = make_male_hemizygous(sample_data.get("GT"), sample_data.phased)

                # Ensure diploidy
                current_gt = sample_data.get("GT")
                if current_gt is None:
                    sample_data["GT"] = (None, None)
                elif len(current_gt) == 1:
                    sample_data["GT"] = (None, current_gt[0])

                # Right align unphased calls
                if not sample_data.phased:
                    sample_data["GT"] = right_align_unphased(sample_data.get("GT"))

    # Recompute allele counts from the final genotypes so the singleton filters see current values
    if compute_allele_counts:
        compute_record_allele_counts(record)

    # Decrement TR IDs
    if decrement_trv_ids:
        allele_type = record.info.get("allele_type")
        if allele_type == "trv":
            if record.id:
                record.id = decrement_trv_id(record.id)
        else:
            trid = record.info.get("TRID")
            if trid:
                record.info["TRID"] = decrement_trv_id(trid)

    # Revise MEIs outside the expected size range to indels
    if prune_meis:
        prune_record_meis(record)

    # Flag homopolymer TR variants
    if flag_homopolymer_trvs and get_scalar(record.info.get("allele_type")) == "trv" and shortest_motif_length(record) == 1:
        record.info["HOMOPOLYMER_TRV"] = True

    # Flag variants with single read support
    if filter_singletons and has_single_read_support(record):
        record.filter.add("SINGLE_READ_SUPPORT")

    # Flag singleton ALT alleles supported only by assembly-based callers
    matching_assembly_alt_idx = []
    if filter_assembly_only_singletons:
        matching_assembly_alt_idx = assembly_only_singleton_alt_idx(record)
        if matching_assembly_alt_idx:
            record.filter.add("ASSEMBLY_ONLY_SINGLETON")

    if not sort_records and not reassign_suffixes:
        output_writer.write(record)
        write_assembly_only_singleton_rows(record, matching_assembly_alt_idx, tsv_writer)
        continue

    assembly_alt_idx[id(record)] = matching_assembly_alt_idx

    # Buffer records sharing the exact same coordinate to sort and/or reassign ID suffixes
    if record.chrom != current_chrom or record.pos != current_pos:
        if buffer:
            flush_buffer(buffer, output_writer, tsv_writer, assembly_alt_idx, reassign_suffixes, sort_records)
        buffer = [record]
        current_chrom = record.chrom
        current_pos = record.pos
    else:
        buffer.append(record)

if (sort_records or reassign_suffixes) and buffer:
    flush_buffer(buffer, output_writer, tsv_writer, assembly_alt_idx, reassign_suffixes, sort_records)

base_reader.close()
output_writer.close()
tsv_writer.close()
if transfer_reader is not None:
    transfer_reader.close()
CODE

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File assembly_only_singletons_tsv = "~{prefix}.assembly_only_singletons.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(size(base_vcf, "GB")) + 5,
        disk_gb: 5 * ceil(size(base_vcf, "GB") + size(transfer_vcf, "GB")) + 25,
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

task DetectSampleOverlaps {
    input {
        File vcf
        File vcf_idx
        String sample
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import subprocess

TYPE_RANK = {'DEL': 2, 'INS': 1, 'SNV': 0, 'OTHER': -1}

def classify_type(at, al_str):
    at = str(at).lower()
    try:
        al = abs(int(al_str))
    except (ValueError, TypeError):
        al = None
    if al == 0:
        return 'SNV'
    if 'del' in at:
        return 'DEL'
    if 'ins' in at or 'dup' in at:
        return 'INS'
    return 'OTHER'

def get_length_bin(al_str):
    try:
        al = abs(int(al_str))
    except (ValueError, TypeError):
        return 'unknown'
    if al == 0:
        return 'SNV'
    elif al <= 29:
        return '1-29bp'
    elif al <= 49:
        return '30-49bp'
    elif al <= 500:
        return '50-500bp'
    elif al <= 5000:
        return '500bp-5kbp'
    elif al <= 50000:
        return '5kbp-50kbp'
    else:
        return '50kbp+'

def wins(iv_i, iv_j):
    """Full tie-breaking cascade: ref span → GQ → allele_length → type rank → QUAL → file order."""
    _, _, t_i, _, _, ref_i, gq_i, al_i, qual_i, idx_i = iv_i
    _, _, t_j, _, _, ref_j, gq_j, al_j, qual_j, idx_j = iv_j
    if len(ref_i) != len(ref_j):
        return len(ref_i) > len(ref_j)
    if gq_i != gq_j:
        return gq_i > gq_j
    if al_i != al_j:
        return al_i > al_j
    ri, rj = TYPE_RANK.get(t_i, -1), TYPE_RANK.get(t_j, -1)
    if ri != rj:
        return ri > rj
    if qual_i != qual_j:
        return qual_i > qual_j
    return idx_i < idx_j

sample_name = '~{sample}'

proc = subprocess.Popen(
    ['bcftools', 'query',
     '-f', '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/allele_type\t%INFO/allele_length\t%QUAL[\t%GT\t%GQ]\n',
     '~{vcf}'],
    stdout=subprocess.PIPE,
    text=True,
    bufsize=1
)

hap0 = []
hap1 = []
n_records = 0

for line in proc.stdout:
    row = line.rstrip('\n').split('\t')
    if len(row) < 10:
        continue
    pos1 = int(row[1])
    ref = row[2].upper()
    vid = row[4]
    at = row[5]
    al_str = row[6]
    qual_str = row[7]
    gt_str = row[8]
    gq_str = row[9]
    if '|' not in gt_str:
        continue
    parts = gt_str.split('|')
    try:
        a0 = int(parts[0]) if parts[0] not in ('.', '') else -1
        a1 = int(parts[1]) if parts[1] not in ('.', '') else -1
    except (ValueError, IndexError):
        continue
    try:
        gq = int(gq_str)
    except (ValueError, TypeError):
        gq = 0
    try:
        al_int = abs(int(al_str))
    except (ValueError, TypeError):
        al_int = 0
    try:
        qual = float(qual_str)
    except (ValueError, TypeError):
        qual = 0.0
    iv = (pos1, pos1 + len(ref) - 1, classify_type(at, al_str), get_length_bin(al_str),
          vid, ref, gq, al_int, qual, n_records)
    if a0 > 0:
        hap0.append(iv)
    if a1 > 0:
        hap1.append(iv)
    n_records += 1

proc.wait()

overlaps_rows = []
cleared_set = set()

def find_overlaps(intervals, haplotype):
    sorted_ivs = sorted(intervals, key=lambda x: x[0])
    n = len(sorted_ivs)
    for i in range(n):
        s_i, e_i = sorted_ivs[i][0], sorted_ivs[i][1]
        for j in range(i + 1, n):
            if sorted_ivs[j][0] > e_i:
                break
            vid_i, t_i, lb_i = sorted_ivs[i][4], sorted_ivs[i][2], sorted_ivs[i][3]
            vid_j, t_j, lb_j = sorted_ivs[j][4], sorted_ivs[j][2], sorted_ivs[j][3]
            if wins(sorted_ivs[i], sorted_ivs[j]):
                winner_id, winner_type, winner_lb = vid_i, t_i, lb_i
                loser_id, loser_type, loser_lb = vid_j, t_j, lb_j
            else:
                winner_id, winner_type, winner_lb = vid_j, t_j, lb_j
                loser_id, loser_type, loser_lb = vid_i, t_i, lb_i
            cleared_set.add(loser_id)
            overlaps_rows.append([sample_name, haplotype,
                                   winner_id, winner_type, winner_lb,
                                   loser_id, loser_type, loser_lb])

if len(hap0) >= 2:
    find_overlaps(hap0, 1)
if len(hap1) >= 2:
    find_overlaps(hap1, 2)

with open('~{prefix}.overlaps.tsv', 'w') as f:
    f.write('sample\thaplotype\tvariant_id_retained\tvar_type_retained\tsize_bin_retained\tvariant_id_cleared\tvar_type_cleared\tsize_bin_cleared\n')
    for row in overlaps_rows:
        f.write('\t'.join(str(x) for x in row) + '\n')

with open('~{prefix}.cleared_calls.tsv', 'w') as f:
    for vid in sorted(cleared_set):
        f.write(f'{vid}\t{sample_name}\n')

print(f'{sample_name}: {n_records} records, {len(overlaps_rows)} pairs, {len(cleared_set)} calls to clear')
CODE
    >>>

    output {
        File overlaps_tsv = "~{prefix}.overlaps.tsv"
        File cleared_calls_tsv = "~{prefix}.cleared_calls.tsv"
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

task ResolveAndNormalizeTask {
    input {
        File vcf
        File vcf_idx
        File? cleared_calls_tsv
        Boolean resolve_overlaps
        Boolean normalize_allele_types
        Int dup_breakpoint_window
        Float dup_size_similarity
        Int min_dup_size
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import re
import pysam

resolve_overlaps = ~{true="True" false="False" resolve_overlaps}
normalize_allele_types = ~{true="True" false="False" normalize_allele_types}
dup_breakpoint_window = ~{dup_breakpoint_window}
dup_size_similarity = ~{dup_size_similarity}
min_dup_size = ~{min_dup_size}

INS_SUBTYPES = {"complex_dup", "dup_interspersed", "inv_dup", "alu_ins", "line_ins", "sva_ins", "numt"}
DEL_SUBTYPES = {"alu_del", "line_del", "sva_del"}

def parse_origin(origin):
    m = re.match(r'^(.+):(\d+)-(\d+)_?[+-]$', origin)
    if m is None:
        return None
    return m.group(1), int(m.group(2)), int(m.group(3))

def passes_criteria(pos, dup_length, origin_start, origin_end):
    origin_length = origin_end - origin_start
    larger = max(dup_length, origin_length)
    if larger == 0 or min(dup_length, origin_length) / larger < dup_size_similarity:
        return False
    if origin_start <= pos <= origin_end:
        return True
    return abs(origin_start - pos) <= dup_breakpoint_window or abs(origin_end - pos) <= dup_breakpoint_window

def is_tandem(record):
    allele_length = record.info.get("allele_length")
    if isinstance(allele_length, (list, tuple)):
        allele_length = allele_length[0]
    if allele_length is None or abs(int(allele_length)) < min_dup_size:
        return False

    origin = record.info.get("ORIGIN")
    if origin is None:
        return False
    origin = origin if isinstance(origin, str) else ",".join(origin)
    if "," in origin:
        return False

    parsed = parse_origin(origin)
    if parsed is None:
        return False

    origin_chrom, origin_start, origin_end = parsed
    dup_length = abs(int(allele_length))

    return passes_criteria(record.pos, dup_length, origin_start, origin_end)

def clear_call(s, record):
    s['GT'] = (None, None)
    s.phased = False
    if 'GQ' in record.format:
        s['GQ'] = None
    if 'EV' in record.format:
        s['EV'] = '.'
    if 'BEV' in record.format:
        s['BEV'] = '.'
    for field in ['AD', 'PL']:
        if field in record.format:
            try:
                current = s.get(field)
                if current is not None:
                    s[field] = tuple(None for _ in current)
                else:
                    s[field] = None
            except Exception:
                pass

def normalize(record):
    allele_type = record.info.get("allele_type")
    if isinstance(allele_type, (list, tuple)):
        allele_type = allele_type[0]

    if allele_type == "dup":
        record.info["allele_subtype"] = "tandem_dup"
        if not is_tandem(record):
            record.info["allele_type"] = "ins"
    elif allele_type in INS_SUBTYPES:
        record.info["allele_subtype"] = allele_type
        record.info["allele_type"] = "ins"
    elif allele_type in DEL_SUBTYPES:
        record.info["allele_subtype"] = allele_type
        record.info["allele_type"] = "del"

to_clear = {}
if resolve_overlaps:
    with open("~{cleared_calls_tsv}") as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            vid, sample = line.split('\t')
            to_clear.setdefault(vid, set()).add(sample)

vcf_in = pysam.VariantFile("~{vcf}")
header = vcf_in.header.copy()
if normalize_allele_types:
    header.info.add("allele_subtype", 1, "String", "Original allele_type subtype for reclassified variants: tandem_dup, dup_interspersed, complex_dup, inv_dup, alu_ins, line_ins, sva_ins, numt, alu_del, line_del or sva_del.")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)
n_cleared = 0

for record in vcf_in:
    record.translate(vcf_out.header)

    # Overlaps are resolved before reclassification, since detection typed variants by their original allele_type
    for sample_name in to_clear.get(record.id, ()):
        if sample_name in record.samples:
            clear_call(record.samples[sample_name], record)
            n_cleared += 1

    if normalize_allele_types:
        normalize(record)

    vcf_out.write(record)

vcf_in.close()
vcf_out.close()
print(f'Cleared {n_cleared:,} sample calls')
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: ceil(size(vcf, "GB")) + 4,
        disk_gb: 3 * ceil(size(vcf, "GB")) + 10,
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
