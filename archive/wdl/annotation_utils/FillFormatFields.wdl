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

        Int? shard_bin_size

        Boolean fill_alt_gts = false
        Boolean fill_ref_gts = false
        Boolean unphase_gts = false
        Boolean add_pl = false

        String utils_docker

        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_unfilled
        RuntimeAttr? runtime_attr_subset_filled
        RuntimeAttr? runtime_attr_fill
        RuntimeAttr? runtime_attr_concat
    }

    if (defined(shard_bin_size)) {
        call Helpers.CreateContigShards {
            input:
                vcfs = [unfilled_vcf, filled_vcf],
                vcf_idxs = [unfilled_vcf_idx, filled_vcf_idx],
                contig = contig,
                shard_bin_size = select_first([shard_bin_size]),
                prefix = "~{prefix}.shards",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_create_shards
        }

        scatter (i in range(length(CreateContigShards.shard_regions))) {
            String shard_region = CreateContigShards.shard_regions[i]

            call Helpers.SubsetVcfToRegion as SubsetUnfilled {
                input:
                    vcf = unfilled_vcf,
                    vcf_idx = unfilled_vcf_idx,
                    region = shard_region,
                    prefix = "~{prefix}.shard_~{i}.unfilled",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_unfilled
            }

            call Helpers.SubsetVcfToRegion as SubsetFilled {
                input:
                    vcf = filled_vcf,
                    vcf_idx = filled_vcf_idx,
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

    if (!defined(shard_bin_size)) {
        call FillVcfFormatFields as FillVcfFormatFieldsNoSharding {
            input:
                unfilled_vcf = unfilled_vcf,
                unfilled_vcf_idx = unfilled_vcf_idx,
                filled_vcf = filled_vcf,
                filled_vcf_idx = filled_vcf_idx,
                format_fields = format_fields,
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
        Boolean fill_alt_gts = false
        Boolean fill_ref_gts = false
        Boolean unphase_gts = false
        Boolean add_pl = false
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

def alleles_key(rec):
    ref = rec.ref.upper() if rec.ref else rec.ref
    alts = tuple(a.upper() for a in rec.alts) if rec.alts else ()
    return (rec.chrom, rec.pos, ref, alts)

for unfilled_rec in unfilled_in:
    # Find matching variant
    unfilled_rec.translate(out_header)
    match = None
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