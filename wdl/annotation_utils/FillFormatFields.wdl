version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FillFormatFields {
    input {
        File unfilled_vcf
        File unfilled_vcf_idx
        File filled_vcf
        File filled_vcf_idx
        String contig
        String prefix

        Array[String] format_fields
        String? subset_unfilled_vcf_field
        String? subset_unfilled_vcf_value

        File? ref_fa
        File? ref_fai

        Int? records_per_shard_normalize
        Int? shard_bin_size_fill

        Boolean normalize_unfilled_vcf
        Boolean normalize_filled_vcf
        Boolean fill_alt_gts
        Boolean fill_ref_gts
        Boolean unphase_gts
        Boolean add_pl
        Boolean match_by_id

        String utils_docker

        RuntimeAttr? runtime_attr_shard_normalize_unfilled
        RuntimeAttr? runtime_attr_normalize_unfilled
        RuntimeAttr? runtime_attr_concat_normalize_unfilled
        RuntimeAttr? runtime_attr_shard_normalize_filled
        RuntimeAttr? runtime_attr_normalize_filled
        RuntimeAttr? runtime_attr_concat_normalize_filled

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_unfilled
        RuntimeAttr? runtime_attr_subset_filled
        RuntimeAttr? runtime_attr_fill
        RuntimeAttr? runtime_attr_concat
    }

    # Optionally normalize the unfilled VCF, sharded by record count so normalization
    # never has to run over a whole-contig VCF at once.
    if (normalize_unfilled_vcf) {
        if (defined(records_per_shard_normalize)) {
            call Helpers.ShardVcfByRecords as ShardUnfilledForNormalize {
                input:
                    vcf = unfilled_vcf,
                    vcf_idx = unfilled_vcf_idx,
                    records_per_shard = select_first([records_per_shard_normalize]),
                    prefix = "~{prefix}.unfilled.norm_shard",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard_normalize_unfilled
            }
        }

        Array[File] unfilled_vcfs_to_normalize = select_first([ShardUnfilledForNormalize.shards, [unfilled_vcf]])
        Array[File] unfilled_vcf_idxs_to_normalize = select_first([ShardUnfilledForNormalize.shard_idxs, [unfilled_vcf_idx]])

        scatter (i in range(length(unfilled_vcfs_to_normalize))) {
            call Helpers.NormalizeVcf as NormalizeUnfilled {
                input:
                    vcf = unfilled_vcfs_to_normalize[i],
                    vcf_idx = unfilled_vcf_idxs_to_normalize[i],
                    ref_fa = select_first([ref_fa]),
                    ref_fai = select_first([ref_fai]),
                    prefix = "~{prefix}.unfilled.norm_shard_~{i}.normalized",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_normalize_unfilled
            }
        }

        # Normalization can shift a variant's position (e.g. splitting a multiallelic),
        # so shards must be re-concatenated with sorting before being re-binned for matching.
        call Helpers.ConcatVcfs as ConcatNormalizedUnfilled {
            input:
                vcfs = NormalizeUnfilled.normalized_vcf,
                vcf_idxs = NormalizeUnfilled.normalized_vcf_idx,
                allow_overlaps = true,
                naive = false,
                sort_output = true,
                prefix = "~{prefix}.unfilled.normalized",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_normalize_unfilled
        }
    }

    File final_unfilled_vcf = select_first([ConcatNormalizedUnfilled.concat_vcf, unfilled_vcf])
    File final_unfilled_vcf_idx = select_first([ConcatNormalizedUnfilled.concat_vcf_idx, unfilled_vcf_idx])

    # Same optional, sharded normalization for the filled VCF.
    if (normalize_filled_vcf) {
        if (defined(records_per_shard_normalize)) {
            call Helpers.ShardVcfByRecords as ShardFilledForNormalize {
                input:
                    vcf = filled_vcf,
                    vcf_idx = filled_vcf_idx,
                    records_per_shard = select_first([records_per_shard_normalize]),
                    prefix = "~{prefix}.filled.norm_shard",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard_normalize_filled
            }
        }

        Array[File] filled_vcfs_to_normalize = select_first([ShardFilledForNormalize.shards, [filled_vcf]])
        Array[File] filled_vcf_idxs_to_normalize = select_first([ShardFilledForNormalize.shard_idxs, [filled_vcf_idx]])

        scatter (i in range(length(filled_vcfs_to_normalize))) {
            call Helpers.NormalizeVcf as NormalizeFilled {
                input:
                    vcf = filled_vcfs_to_normalize[i],
                    vcf_idx = filled_vcf_idxs_to_normalize[i],
                    ref_fa = select_first([ref_fa]),
                    ref_fai = select_first([ref_fai]),
                    prefix = "~{prefix}.filled.norm_shard_~{i}.normalized",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_normalize_filled
            }
        }

        call Helpers.ConcatVcfs as ConcatNormalizedFilled {
            input:
                vcfs = NormalizeFilled.normalized_vcf,
                vcf_idxs = NormalizeFilled.normalized_vcf_idx,
                allow_overlaps = true,
                naive = false,
                sort_output = true,
                prefix = "~{prefix}.filled.normalized",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_normalize_filled
        }
    }

    File final_filled_vcf = select_first([ConcatNormalizedFilled.concat_vcf, filled_vcf])
    File final_filled_vcf_idx = select_first([ConcatNormalizedFilled.concat_vcf_idx, filled_vcf_idx])

    if (defined(shard_bin_size_fill)) {
        call Helpers.CreateContigShards {
            input:
                vcfs = [final_unfilled_vcf, final_filled_vcf],
                vcf_idxs = [final_unfilled_vcf_idx, final_filled_vcf_idx],
                contig = contig,
                shard_bin_size = select_first([shard_bin_size_fill]),
                prefix = "~{prefix}.shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_shards
        }

        scatter (i in range(length(CreateContigShards.shard_regions))) {
            String shard_region = CreateContigShards.shard_regions[i]

            call Helpers.SubsetVcfToRegion as SubsetUnfilled {
                input:
                    vcf = final_unfilled_vcf,
                    vcf_idx = final_unfilled_vcf_idx,
                    region = shard_region,
                    prefix = "~{prefix}.shard_~{i}.unfilled",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_unfilled
            }

            call Helpers.SubsetVcfToRegion as SubsetFilled {
                input:
                    vcf = final_filled_vcf,
                    vcf_idx = final_filled_vcf_idx,
                    region = shard_region,
                    prefix = "~{prefix}.shard_~{i}.filled_input",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_filled
            }

            call FillVcfFormatFields {
                input:
                    unfilled_vcf = SubsetUnfilled.subset_vcf,
                    unfilled_vcf_idx = SubsetUnfilled.subset_vcf_idx,
                    filled_vcf = SubsetFilled.subset_vcf,
                    filled_vcf_idx = SubsetFilled.subset_vcf_idx,
                    format_fields = format_fields,
                    match_by_id = match_by_id,
                    subset_unfilled_vcf_field = subset_unfilled_vcf_field,
                    subset_unfilled_vcf_value = subset_unfilled_vcf_value,
                    fill_alt_gts = fill_alt_gts,
                    fill_ref_gts = fill_ref_gts,
                    unphase_gts = unphase_gts,
                    add_pl = add_pl,
                    prefix = "~{prefix}.shard_~{i}.filled",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_fill
            }
        }

        call Helpers.ConcatVcfs {
            input:
                vcfs = FillVcfFormatFields.output_vcf,
                vcf_idxs = FillVcfFormatFields.output_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.filled",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    if (!defined(shard_bin_size_fill)) {
        call FillVcfFormatFields as FillVcfFormatFieldsNoSharding {
            input:
                unfilled_vcf = final_unfilled_vcf,
                unfilled_vcf_idx = final_unfilled_vcf_idx,
                filled_vcf = final_filled_vcf,
                filled_vcf_idx = final_filled_vcf_idx,
                format_fields = format_fields,
                match_by_id = match_by_id,
                subset_unfilled_vcf_field = subset_unfilled_vcf_field,
                subset_unfilled_vcf_value = subset_unfilled_vcf_value,
                fill_alt_gts = fill_alt_gts,
                fill_ref_gts = fill_ref_gts,
                unphase_gts = unphase_gts,
                add_pl = add_pl,
                prefix = prefix,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_fill
        }
    }

    output {
        File refilled_vcf = select_first([ConcatVcfs.concat_vcf, FillVcfFormatFieldsNoSharding.output_vcf])
        File refilled_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, FillVcfFormatFieldsNoSharding.output_vcf_idx])
    }
}

task FillVcfFormatFields {
    input {
        File unfilled_vcf
        File unfilled_vcf_idx
        File filled_vcf
        File filled_vcf_idx
        Array[String] format_fields
        Boolean match_by_id
        String? subset_unfilled_vcf_field
        String? subset_unfilled_vcf_value
        Boolean fill_alt_gts
        Boolean fill_ref_gts
        Boolean unphase_gts
        Boolean add_pl
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        format_fields_file="~{write_lines(format_fields)}"

        if [[ "~{filled_vcf_idx}" != "~{filled_vcf}.tbi" ]]; then
            ln -sf "~{filled_vcf_idx}" "~{filled_vcf}.tbi"
        fi

        python3 <<CODE
import math
import pysam

with open("$format_fields_file") as fh:
    format_fields = [l.strip() for l in fh if l.strip()]

assert "GT" not in format_fields, "GT must not be passed in format_fields; use fill_alt_gts/fill_ref_gts instead"

match_by_id = ~{true="True" false="False" match_by_id}
subset_field = ~{if defined(subset_unfilled_vcf_field) then "'" + subset_unfilled_vcf_field + "'" else "None"}
subset_value = ~{if defined(subset_unfilled_vcf_value) then "'" + subset_unfilled_vcf_value + "'" else "None"}
fill_alt_gts = ~{true="True" false="False" fill_alt_gts}
fill_ref_gts = ~{true="True" false="False" fill_ref_gts}
unphase_gts = ~{true="True" false="False" unphase_gts}
add_pl = ~{true="True" false="False" add_pl}

unfilled_in = pysam.VariantFile("~{unfilled_vcf}")
filled_in = pysam.VariantFile("~{filled_vcf}")

out_header = unfilled_in.header.copy()
for field in format_fields:
    if field in filled_in.header.formats and field not in out_header.formats:
        out_header.add_record(filled_in.header.formats[field].record)
if add_pl and "PL" not in out_header.formats:
    out_header.add_line('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">')

filled_samples = set(filled_in.header.samples)
common_samples = [s for s in out_header.samples if s in filled_samples]
all_samples = list(out_header.samples)

out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=out_header)

def gt_is_alt(gt):
    return any(a is not None and a > 0 for a in gt)

def get_sample_ploidy(sample_data):
    gt = sample_data.get("GT")
    if gt is None:
        return None
    return len(gt)

def calculate_pl(ref_reads, alt_reads, ploidy):
    if ref_reads + alt_reads == 0:
        return (0, 0) if ploidy == 1 else (0, 0, 0)
    genotype_error = 0.05
    if ploidy == 1:
        means = [genotype_error, 1.0 - genotype_error]
    else:
        means = [genotype_error, 1.0 / ploidy, 1.0 - genotype_error]
    log10 = math.log(10)
    ll = [(alt_reads * math.log(p) + ref_reads * math.log(1.0 - p)) / log10 for p in means]
    max_ll = max(ll)
    return tuple(int(round(-10 * (x - max_ll))) for x in ll)

def ad_is_populated(ad):
    return ad is not None and len(ad) == 2 and all(value is not None for value in ad)

def unphase_gt(gt):
    return tuple(sorted(gt, key=lambda allele: (allele is None, allele if allele is not None else 0)))

def in_subset(rec):
    if subset_field is None:
        return True
    info_val = rec.info.get(subset_field)
    if info_val is None:
        return False
    if isinstance(info_val, (list, tuple)):
        return subset_value in [str(v) for v in info_val]
    return str(info_val) == subset_value

def alleles_key(rec):
    ref = rec.ref.upper() if rec.ref else rec.ref
    alts = tuple(a.upper() for a in rec.alts) if rec.alts else ()
    if match_by_id:
        return (rec.chrom, rec.pos, ref, alts, rec.id)
    return (rec.chrom, rec.pos, ref, alts)

for unfilled_rec in unfilled_in:
    # Find matching variant
    unfilled_rec.translate(out_header)
    match = None
    if in_subset(unfilled_rec):
        unfilled_key = alleles_key(unfilled_rec)
        for cand in filled_in.fetch(unfilled_rec.chrom, unfilled_rec.start, unfilled_rec.stop):
            if alleles_key(cand) == unfilled_key:
                match = cand
                break

    if match:
        for sample in common_samples:
            # Set GT field
            if fill_alt_gts or fill_ref_gts:
                src_gt = match.samples[sample].get("GT")
                cur_gt = unfilled_rec.samples[sample].get("GT")
                if src_gt is not None and cur_gt is not None:
                    cur_is_alt = gt_is_alt(cur_gt)
                    if (fill_alt_gts and cur_is_alt) or (fill_ref_gts and not cur_is_alt):
                        unfilled_rec.samples[sample]["GT"] = src_gt
                        unfilled_rec.samples[sample].phased = match.samples[sample].phased

            # Copy over values for format fields
            for field in format_fields:
                if field not in match.format:
                    continue
                value = match.samples[sample].get(field)
                try:
                    unfilled_rec.samples[sample][field] = value
                except Exception:
                    print(f"[{unfilled_rec.id}] Could not set {field} for {sample}.")
                    pass

    # Unphase genotypes
    if unphase_gts:
        for sample in all_samples:
            current_gt = unfilled_rec.samples[sample].get("GT")
            if current_gt is not None:
                unfilled_rec.samples[sample]["GT"] = unphase_gt(current_gt)
                unfilled_rec.samples[sample].phased = False

    # Set PL
    if add_pl and "PL" not in unfilled_rec.format and "AD" in unfilled_rec.format:
        for sample in all_samples:
            ad = unfilled_rec.samples[sample].get("AD")
            if ad_is_populated(ad):
                ploidy = get_sample_ploidy(unfilled_rec.samples[sample])
                if ploidy is not None:
                    unfilled_rec.samples[sample]["PL"] = calculate_pl(ad[0], ad[1], ploidy)

    out.write(unfilled_rec)

out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(unfilled_vcf, "GB") + size(filled_vcf, "GB")) + 10,
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
