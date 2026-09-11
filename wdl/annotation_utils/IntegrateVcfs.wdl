version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow IntegrateVcfs {
    input {
        # SNV/indel inputs
        File? snv_indel_vcf
        File? snv_indel_vcf_idx
        Boolean? normalize_snv_indel_vcf
        String? snv_indel_vcf_source_tag
        String? size_filter_snv_indel_vcf
        String? size_filter_snv_indel_vcf_description
        File? swap_samples_snv_indel

        # SV inputs
        File? sv_vcf
        File? sv_vcf_idx
        Boolean? normalize_sv_vcf
        String? sv_vcf_source_tag
        String? size_filter_sv_vcf
        String? size_filter_sv_vcf_description
        File? swap_samples_sv

        # Shared inputs
        File ref_fa
        File ref_fai
        Array[String] contigs
        String prefix
        Int min_sv_length
        Int? records_per_shard
        Array[String]? sample_ids
        String utils_docker
        RuntimeAttr? runtime_attr_swap_samples_snv_indel
        RuntimeAttr? runtime_attr_subset_contig_snv_indel
        RuntimeAttr? runtime_attr_shard_snv_indel
        RuntimeAttr? runtime_attr_normalize_snv_indel
        RuntimeAttr? runtime_attr_subset_samples_snv_indel
        RuntimeAttr? runtime_attr_annotate_attributes_snv_indel
        RuntimeAttr? runtime_attr_add_info_snv_indel
        RuntimeAttr? runtime_attr_add_filter_snv_indel
        RuntimeAttr? runtime_attr_concat_snv_indel_shards
        RuntimeAttr? runtime_attr_swap_samples_sv
        RuntimeAttr? runtime_attr_subset_contig_sv
        RuntimeAttr? runtime_attr_shard_sv
        RuntimeAttr? runtime_attr_normalize_sv
        RuntimeAttr? runtime_attr_subset_samples_sv
        RuntimeAttr? runtime_attr_annotate_attributes_sv
        RuntimeAttr? runtime_attr_add_info_sv
        RuntimeAttr? runtime_attr_add_filter_sv
        RuntimeAttr? runtime_attr_concat_sv_shards
        RuntimeAttr? runtime_attr_validate_inputs
        RuntimeAttr? runtime_attr_get_samples
        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_merge
        RuntimeAttr? runtime_attr_rename_and_filter
        RuntimeAttr? runtime_attr_concat
    }

    call ValidateIntegrateVcfsInputs {
        input:
            snv_vcf = defined(snv_indel_vcf), snv_idx = defined(snv_indel_vcf_idx), normalize_snv = normalize_snv_indel_vcf,
            snv_tag = snv_indel_vcf_source_tag, snv_filter = size_filter_snv_indel_vcf, snv_filter_description = size_filter_snv_indel_vcf_description,
            sv_vcf = defined(sv_vcf), sv_idx = defined(sv_vcf_idx), normalize_sv = normalize_sv_vcf,
            sv_tag = sv_vcf_source_tag, sv_filter = size_filter_sv_vcf, sv_filter_description = size_filter_sv_vcf_description,
            docker = utils_docker, runtime_attr_override = runtime_attr_validate_inputs
    }
    Boolean inputs_valid = ValidateIntegrateVcfsInputs.status == "success"
    Boolean has_snv = inputs_valid && defined(snv_indel_vcf) && defined(snv_indel_vcf_idx)
    Boolean has_sv = inputs_valid && defined(sv_vcf) && defined(sv_vcf_idx)

    if (has_snv && defined(swap_samples_snv_indel)) {
        call Helpers.SwapSampleIds as SwapSnvIndel {
            input: vcf = select_first([snv_indel_vcf]), vcf_idx = select_first([snv_indel_vcf_idx]), sample_swap_list = select_first([swap_samples_snv_indel]), prefix = "~{prefix}.snv_indel.swapped", docker = utils_docker, runtime_attr_override = runtime_attr_swap_samples_snv_indel
        }
    }
    if (has_sv && defined(swap_samples_sv)) {
        call Helpers.SwapSampleIds as SwapSv {
            input: vcf = select_first([sv_vcf]), vcf_idx = select_first([sv_vcf_idx]), sample_swap_list = select_first([swap_samples_sv]), prefix = "~{prefix}.sv.swapped", docker = utils_docker, runtime_attr_override = runtime_attr_swap_samples_sv
        }
    }
    File? final_snv_vcf = select_first([SwapSnvIndel.swapped_vcf, snv_indel_vcf])
    File? final_snv_idx = select_first([SwapSnvIndel.swapped_vcf_idx, snv_indel_vcf_idx])
    File? final_sv_vcf = select_first([SwapSv.swapped_vcf, sv_vcf])
    File? final_sv_idx = select_first([SwapSv.swapped_vcf_idx, sv_vcf_idx])
    Array[File] input_vcfs = select_all([final_snv_vcf, final_sv_vcf])
    Array[File] input_idxs = select_all([final_snv_idx, final_sv_idx])

    if (inputs_valid && !defined(sample_ids)) {
        call Helpers.GetSamplesFromVcf {
            input: vcf = select_first([final_snv_vcf, final_sv_vcf]), vcf_idx = select_first([final_snv_idx, final_sv_idx]), docker = utils_docker, runtime_attr_override = runtime_attr_get_samples
        }
        call Helpers.CheckSampleConsistency as CheckInputSampleConsistency {
            input: vcfs = input_vcfs, vcf_idxs = input_idxs, sample_ids = GetSamplesFromVcf.samples, docker = utils_docker, runtime_attr_override = runtime_attr_check_samples
        }
    }
    Array[String] final_sample_ids = select_first([sample_ids, GetSamplesFromVcf.samples])
    Array[String] checked_contigs = if (inputs_valid && select_first([CheckInputSampleConsistency.status, "success"]) == "success") then contigs else []
    Boolean single_contig = length(contigs) == 1

    scatter (contig in checked_contigs) {
        if (has_snv) {
            if (!single_contig) {
                call Helpers.SubsetVcfToContig as SubsetSnv {
                    input: vcf = select_first([final_snv_vcf]), vcf_idx = select_first([final_snv_idx]), contig = contig, prefix = "~{prefix}.~{contig}.snv_indel", docker = utils_docker, runtime_attr_override = runtime_attr_subset_contig_snv_indel
                }
            }
            File snv_contig_vcf = select_first([SubsetSnv.subset_vcf, final_snv_vcf])
            File snv_contig_idx = select_first([SubsetSnv.subset_vcf_idx, final_snv_idx])
            if (defined(records_per_shard)) {
                call Helpers.ShardVcfByRecords as ShardSnv { input: vcf = snv_contig_vcf, vcf_idx = snv_contig_idx, records_per_shard = select_first([records_per_shard]), prefix = "~{prefix}.~{contig}.snv_indel", docker = utils_docker, runtime_attr_override = runtime_attr_shard_snv_indel }
            }
            Array[File] snv_vcfs = select_first([ShardSnv.shards, [snv_contig_vcf]])
            Array[File] snv_idxs = select_first([ShardSnv.shard_idxs, [snv_contig_idx]])
            scatter (i in range(length(snv_vcfs))) {
                if (select_first([normalize_snv_indel_vcf])) {
                    call Helpers.NormalizeVcf as NormalizeSnv { input: vcf = snv_vcfs[i], vcf_idx = snv_idxs[i], ref_fa = ref_fa, ref_fai = ref_fai, prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.normalized", docker = utils_docker, runtime_attr_override = runtime_attr_normalize_snv_indel }
                }
                File normalized_snv = select_first([NormalizeSnv.normalized_vcf, snv_vcfs[i]])
                File normalized_snv_idx = select_first([NormalizeSnv.normalized_vcf_idx, snv_idxs[i]])
                if (defined(sample_ids)) {
                    call Helpers.SubsetVcfToSamples as SubsetSamplesSnv { input: vcf = normalized_snv, vcf_idx = normalized_snv_idx, samples = final_sample_ids, prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.subset", docker = utils_docker, runtime_attr_override = runtime_attr_subset_samples_snv_indel }
                }
                File processed_snv = select_first([SubsetSamplesSnv.subset_vcf, normalized_snv])
                File processed_snv_idx = select_first([SubsetSamplesSnv.subset_vcf_idx, normalized_snv_idx])
                call Helpers.AnnotateVariantAttributes as AnnotateSnv { input: vcf = processed_snv, vcf_idx = processed_snv_idx, prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.annotated", docker = utils_docker, runtime_attr_override = runtime_attr_annotate_attributes_snv_indel }
                call Helpers.AddInfo as AddInfoSnv { input: vcf = AnnotateSnv.annotated_vcf, vcf_idx = AnnotateSnv.annotated_vcf_idx, tag_id = "SOURCE", tag_value = select_first([snv_indel_vcf_source_tag]), tag_description = "Source of variant call", prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.add_info", docker = utils_docker, runtime_attr_override = runtime_attr_add_info_snv_indel }
                call Helpers.AddFilter as FilterSnv { input: vcf = AddInfoSnv.annotated_vcf, vcf_idx = AddInfoSnv.annotated_vcf_idx, filter_name = select_first([size_filter_snv_indel_vcf]), filter_description = select_first([size_filter_snv_indel_vcf_description]), filter_expression = "abs(INFO/allele_length) >= ~{min_sv_length}", prefix = "~{prefix}.~{contig}.snv_indel.shard_~{i}.add_filter", docker = utils_docker, runtime_attr_override = runtime_attr_add_filter_snv_indel }
            }
            if (defined(records_per_shard)) { call Helpers.ConcatVcfs as ConcatSnv { input: vcfs = FilterSnv.flagged_vcf, vcf_idxs = FilterSnv.flagged_vcf_idx, allow_overlaps = false, naive = true, prefix = "~{prefix}.~{contig}.snv_indel.concatenated", docker = utils_docker, runtime_attr_override = runtime_attr_concat_snv_indel_shards } }
            File final_snv_contig_vcf = select_first([ConcatSnv.concat_vcf, FilterSnv.flagged_vcf[0]])
            File final_snv_contig_idx = select_first([ConcatSnv.concat_vcf_idx, FilterSnv.flagged_vcf_idx[0]])
        }
        if (has_sv) {
            if (!single_contig) { call Helpers.SubsetVcfToContig as SubsetSv { input: vcf = select_first([final_sv_vcf]), vcf_idx = select_first([final_sv_idx]), contig = contig, prefix = "~{prefix}.~{contig}.sv", docker = utils_docker, runtime_attr_override = runtime_attr_subset_contig_sv } }
            File sv_contig_vcf = select_first([SubsetSv.subset_vcf, final_sv_vcf])
            File sv_contig_idx = select_first([SubsetSv.subset_vcf_idx, final_sv_idx])
            if (defined(records_per_shard)) { call Helpers.ShardVcfByRecords as ShardSv { input: vcf = sv_contig_vcf, vcf_idx = sv_contig_idx, records_per_shard = select_first([records_per_shard]), prefix = "~{prefix}.~{contig}.sv", docker = utils_docker, runtime_attr_override = runtime_attr_shard_sv } }
            Array[File] sv_vcfs = select_first([ShardSv.shards, [sv_contig_vcf]])
            Array[File] sv_idxs = select_first([ShardSv.shard_idxs, [sv_contig_idx]])
            scatter (i in range(length(sv_vcfs))) {
                if (select_first([normalize_sv_vcf])) { call Helpers.NormalizeVcf as NormalizeSv { input: vcf = sv_vcfs[i], vcf_idx = sv_idxs[i], ref_fa = ref_fa, ref_fai = ref_fai, prefix = "~{prefix}.~{contig}.sv.shard_~{i}.normalized", docker = utils_docker, runtime_attr_override = runtime_attr_normalize_sv } }
                File normalized_sv = select_first([NormalizeSv.normalized_vcf, sv_vcfs[i]])
                File normalized_sv_idx = select_first([NormalizeSv.normalized_vcf_idx, sv_idxs[i]])
                if (defined(sample_ids)) { call Helpers.SubsetVcfToSamples as SubsetSamplesSv { input: vcf = normalized_sv, vcf_idx = normalized_sv_idx, samples = final_sample_ids, prefix = "~{prefix}.~{contig}.sv.shard_~{i}.subset", docker = utils_docker, runtime_attr_override = runtime_attr_subset_samples_sv } }
                File processed_sv = select_first([SubsetSamplesSv.subset_vcf, normalized_sv])
                File processed_sv_idx = select_first([SubsetSamplesSv.subset_vcf_idx, normalized_sv_idx])
                call Helpers.AnnotateVariantAttributes as AnnotateSv { input: vcf = processed_sv, vcf_idx = processed_sv_idx, prefix = "~{prefix}.~{contig}.sv.shard_~{i}.annotated", docker = utils_docker, runtime_attr_override = runtime_attr_annotate_attributes_sv }
                call Helpers.AddInfo as AddInfoSv { input: vcf = AnnotateSv.annotated_vcf, vcf_idx = AnnotateSv.annotated_vcf_idx, tag_id = "SOURCE", tag_value = select_first([sv_vcf_source_tag]), tag_description = "Source of variant call", prefix = "~{prefix}.~{contig}.sv.shard_~{i}.add_info", docker = utils_docker, runtime_attr_override = runtime_attr_add_info_sv }
                call Helpers.AddFilter as FilterSv { input: vcf = AddInfoSv.annotated_vcf, vcf_idx = AddInfoSv.annotated_vcf_idx, filter_name = select_first([size_filter_sv_vcf]), filter_description = select_first([size_filter_sv_vcf_description]), filter_expression = "abs(INFO/allele_length) < ~{min_sv_length}", prefix = "~{prefix}.~{contig}.sv.shard_~{i}.add_filter", docker = utils_docker, runtime_attr_override = runtime_attr_add_filter_sv }
            }
            if (defined(records_per_shard)) { call Helpers.ConcatVcfs as ConcatSv { input: vcfs = FilterSv.flagged_vcf, vcf_idxs = FilterSv.flagged_vcf_idx, allow_overlaps = false, naive = true, prefix = "~{prefix}.~{contig}.sv.concatenated", docker = utils_docker, runtime_attr_override = runtime_attr_concat_sv_shards } }
            File final_sv_contig_vcf = select_first([ConcatSv.concat_vcf, FilterSv.flagged_vcf[0]])
            File final_sv_contig_idx = select_first([ConcatSv.concat_vcf_idx, FilterSv.flagged_vcf_idx[0]])
        }
        Array[File] contig_vcfs = select_all([final_snv_contig_vcf, final_sv_contig_vcf])
        Array[File] contig_idxs = select_all([final_snv_contig_idx, final_sv_contig_idx])
        call Helpers.CheckSampleConsistency as CheckContigSamples { input: vcfs = contig_vcfs, vcf_idxs = contig_idxs, sample_ids = final_sample_ids, docker = utils_docker, runtime_attr_override = runtime_attr_check_samples }
        Array[File] checked_contig_vcfs = if (CheckContigSamples.status == "success") then contig_vcfs else []
        Array[File] checked_contig_idxs = if (CheckContigSamples.status == "success") then contig_idxs else []
        if (length(checked_contig_vcfs) > 1) {
            call Helpers.ConcatVcfs as MergeContigVcfs { input: vcfs = checked_contig_vcfs, vcf_idxs = checked_contig_idxs, allow_overlaps = true, naive = false, prefix = "~{prefix}.~{contig}.integrated", docker = utils_docker, runtime_attr_override = runtime_attr_merge }
        }
        File merged_contig_vcf = select_first([MergeContigVcfs.concat_vcf, final_snv_contig_vcf, final_sv_contig_vcf])
        File merged_contig_idx = select_first([MergeContigVcfs.concat_vcf_idx, final_snv_contig_idx, final_sv_contig_idx])
        call RenameAndFilterVariants { input: vcf = merged_contig_vcf, vcf_idx = merged_contig_idx, size_filter_snv_indel_vcf = size_filter_snv_indel_vcf, snv_indel_vcf_source_tag = snv_indel_vcf_source_tag, sv_vcf_source_tag = sv_vcf_source_tag, sample_consistency_status = CheckContigSamples.status, prefix = "~{prefix}.~{contig}.filtered", docker = utils_docker, runtime_attr_override = runtime_attr_rename_and_filter }
    }
    call Helpers.ConcatVcfs { input: vcfs = RenameAndFilterVariants.filtered_vcf, vcf_idxs = RenameAndFilterVariants.filtered_vcf_idx, allow_overlaps = false, naive = true, prefix = "~{prefix}.integrated", docker = utils_docker, runtime_attr_override = runtime_attr_concat }
    output { File integrated_vcf = ConcatVcfs.concat_vcf File integrated_vcf_idx = ConcatVcfs.concat_vcf_idx }
}

task ValidateIntegrateVcfsInputs {
    input {
        Boolean snv_vcf
        Boolean snv_idx
        Boolean? normalize_snv
        String? snv_tag
        String? snv_filter
        String? snv_filter_description
        Boolean sv_vcf
        Boolean sv_idx
        Boolean? normalize_sv
        String? sv_tag
        String? sv_filter
        String? sv_filter_description
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        require() { [ -n "$2" ] || { echo "ERROR: $1 is required when its VCF is supplied" >&2; exit 1; }; }

        if ! ~{snv_vcf} && ! ~{sv_vcf}; then echo "ERROR: supply at least one VCF" >&2; exit 1; fi
        [ "~{snv_vcf}" = "~{snv_idx}" ] || { echo "ERROR: snv_indel_vcf and snv_indel_vcf_idx must be supplied together" >&2; exit 1; }
        [ "~{sv_vcf}" = "~{sv_idx}" ] || { echo "ERROR: sv_vcf and sv_vcf_idx must be supplied together" >&2; exit 1; }
        if ~{snv_vcf}; then require normalize_snv_indel_vcf "~{normalize_snv}"; require snv_indel_vcf_source_tag "~{snv_tag}"; require size_filter_snv_indel_vcf "~{snv_filter}"; require size_filter_snv_indel_vcf_description "~{snv_filter_description}"; fi
        if ~{sv_vcf}; then require normalize_sv_vcf "~{normalize_sv}"; require sv_vcf_source_tag "~{sv_tag}"; require size_filter_sv_vcf "~{sv_filter}"; require size_filter_sv_vcf_description "~{sv_filter_description}"; fi
    >>>

    output {
        String status = "success"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 5,
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

task RenameAndFilterVariants {
    input {
        File vcf
        File vcf_idx
        String? size_filter_snv_indel_vcf
        String? snv_indel_vcf_source_tag
        String? sv_vcf_source_tag
        String sample_consistency_status
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        [ "~{sample_consistency_status}" = "success" ]

        python3 <<CODE
from pysam import VariantFile
from collections import defaultdict

snv_source = "~{snv_indel_vcf_source_tag}"
sv_source = "~{sv_vcf_source_tag}"
size_filter = "~{size_filter_snv_indel_vcf}"
compare = bool(snv_source and sv_source and size_filter)

def key(record):
    ref = record.ref.upper() if record.ref else record.ref
    alts = tuple(alt.upper() for alt in record.alts) if record.alts else ()
    return record.chrom, record.pos, ref, alts

inv = VariantFile("~{vcf}")
sv_variants = {key(record) for record in inv if compare and record.info.get("SOURCE") == sv_source}
inv.close()

def redundant(record):
    return (
        compare
        and record.info.get("SOURCE") == snv_source
        and size_filter in record.filter
        and key(record) in sv_variants
    )

def variant_id(record):
    allele_type = record.info.get("allele_type").upper()
    if allele_type == "SNV":
        return f"{record.chrom}-{record.pos}-{record.ref}-{record.alts[0]}"
    return f"{record.chrom}-{record.pos}-{allele_type}-{abs(int(record.info.get('allele_length')))}"

inv = VariantFile("~{vcf}")
counts = defaultdict(int)
for record in inv:
    if not redundant(record):
        counts[variant_id(record)] += 1
inv.close()

inv = VariantFile("~{vcf}")
out = VariantFile("~{prefix}.vcf.gz", "w", header=inv.header)
seen = defaultdict(int)
for record in inv:
    if redundant(record):
        continue
    name = variant_id(record)
    if counts[name] > 1:
        seen[name] += 1
        record.id = f"{name}_{seen[name]}"
    else:
        record.id = name
    out.write(record)
inv.close()
out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
