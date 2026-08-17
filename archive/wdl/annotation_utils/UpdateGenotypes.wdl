version 1.0

import "../utils/Helpers.wdl"

workflow UpdateGenotypes {
    input {
        File base_vcf
        File base_vcf_idx
        File? phased_vcf
        File? phased_vcf_idx
        Array[String] contigs
        String prefix

        Int? shard_bin_size

        File ped
        Boolean transfer_genotypes = false
        Boolean drop_genotypes = false
        Boolean decrement_trv_ids = false
        Array[String]? unphase_samples
        Array[String]? drop_samples

        String utils_docker

        RuntimeAttr? runtime_attr_drop_samples
        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_subset_base
        RuntimeAttr? runtime_attr_subset_genotyped
        RuntimeAttr? runtime_attr_update_genotypes
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_strip_genotypes
    }

    Boolean single_contig = length(contigs) == 1

    if (defined(drop_samples)) {
        call Helpers.SubsetVcfToSamples as DropBaseSamples {
            input:
                vcf = base_vcf,
                vcf_idx = base_vcf_idx,
                samples = select_first([drop_samples]),
                keep_samples = false,
                filter_to_sample = false,
                prefix = "~{prefix}.drop_samples",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_drop_samples
        }
    }

    File filtered_base_vcf = select_first([DropBaseSamples.subset_vcf, base_vcf])
    File filtered_base_vcf_idx = select_first([DropBaseSamples.subset_vcf_idx, base_vcf_idx])

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig as SubsetBase {
                input:
                    vcf = filtered_base_vcf,
                    vcf_idx = filtered_base_vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.base",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_base
            }
        }

        if (transfer_genotypes && !single_contig) {
            call Helpers.SubsetVcfToContig as SubsetGenotyped {
                input:
                    vcf = select_first([phased_vcf]),
                    vcf_idx = select_first([phased_vcf_idx]),
                    contig = contig,
                    prefix = "~{prefix}.~{contig}.genotyped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_genotyped
            }
        }

        File contig_base_vcf = select_first([SubsetBase.subset_vcf, filtered_base_vcf])
        File contig_base_vcf_idx = select_first([SubsetBase.subset_vcf_idx, filtered_base_vcf_idx])
        File transfer_source_vcf = if transfer_genotypes then select_first([SubsetGenotyped.subset_vcf, select_first([phased_vcf])]) else contig_base_vcf
        File transfer_source_vcf_idx = if transfer_genotypes then select_first([SubsetGenotyped.subset_vcf_idx, select_first([phased_vcf_idx])]) else contig_base_vcf_idx

        if (defined(shard_bin_size)) {
            call Helpers.CreateContigShards {
                input:
                    vcfs = [contig_base_vcf, transfer_source_vcf],
                    vcf_idxs = [contig_base_vcf_idx, transfer_source_vcf_idx],
                    contig = contig,
                    shard_bin_size = select_first([shard_bin_size]),
                    prefix = "~{prefix}.~{contig}.shards",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_create_shards
            }

            scatter (i in range(length(CreateContigShards.shard_regions))) {
                String shard_region = CreateContigShards.shard_regions[i]

                call Helpers.SubsetVcfToRegion as SubsetBaseShard {
                    input:
                        vcf = contig_base_vcf,
                        vcf_idx = contig_base_vcf_idx,
                        region = shard_region,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.base",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_base
                }

                call Helpers.SubsetVcfToRegion as SubsetGenotypedShard {
                    input:
                        vcf = transfer_source_vcf,
                        vcf_idx = transfer_source_vcf_idx,
                        region = shard_region,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.genotyped",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_genotyped
                }

                call UpdateContigGenotypes {
                    input:
                        base_vcf = SubsetBaseShard.subset_vcf,
                        base_vcf_idx = SubsetBaseShard.subset_vcf_idx,
                        phased_vcf = SubsetGenotypedShard.subset_vcf,
                        phased_vcf_idx = SubsetGenotypedShard.subset_vcf_idx,
                        ped = ped,
                        unphase_samples = select_first([unphase_samples, []]),
                        transfer_genotypes = transfer_genotypes,
                        decrement_trv_ids = decrement_trv_ids,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.genotyped",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_update_genotypes
                }
            }

            call Helpers.ConcatVcfs as ConcatContig {
                input:
                    vcfs = UpdateContigGenotypes.genotyped_vcf,
                    vcf_idxs = UpdateContigGenotypes.genotyped_vcf_idx,
                    allow_overlaps = false,
                    naive = true,
                    prefix = "~{prefix}.~{contig}.genotyped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat
            }
        }

        if (!defined(shard_bin_size)) {
            call UpdateContigGenotypes as UpdateContigGenotypesNoSharding {
                input:
                    base_vcf = contig_base_vcf,
                    base_vcf_idx = contig_base_vcf_idx,
                    phased_vcf = transfer_source_vcf,
                    phased_vcf_idx = transfer_source_vcf_idx,
                    ped = ped,
                    unphase_samples = select_first([unphase_samples, []]),
                    transfer_genotypes = transfer_genotypes,
                    decrement_trv_ids = decrement_trv_ids,
                    prefix = "~{prefix}.~{contig}.genotyped",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_update_genotypes
            }
        }

        File genotyped_contig_vcf = select_first([ConcatContig.concat_vcf, UpdateContigGenotypesNoSharding.genotyped_vcf])
        File genotyped_contig_vcf_idx = select_first([ConcatContig.concat_vcf_idx, UpdateContigGenotypesNoSharding.genotyped_vcf_idx])
    }

    if (!single_contig) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = genotyped_contig_vcf,
                vcf_idxs = genotyped_contig_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.genotyped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat
        }
    }

    File merged_genotyped_vcf = select_first([ConcatVcfs.concat_vcf, genotyped_contig_vcf[0]])
    File merged_genotyped_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, genotyped_contig_vcf_idx[0]])

    if (drop_genotypes) {
        call Helpers.StripGenotypes {
            input:
                vcf = merged_genotyped_vcf,
                vcf_idx = merged_genotyped_vcf_idx,
                prefix = "~{prefix}.stripped",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_strip_genotypes
        }
    }

    output {
        File genotyped_vcf = select_first([StripGenotypes.stripped_vcf, merged_genotyped_vcf])
        File genotyped_vcf_idx = select_first([StripGenotypes.stripped_vcf_idx, merged_genotyped_vcf_idx])
    }
}

task UpdateContigGenotypes {
    input {
        File base_vcf
        File base_vcf_idx
        File phased_vcf
        File phased_vcf_idx
        File ped
        Array[String] unphase_samples
        Boolean transfer_genotypes
        Boolean decrement_trv_ids = false
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import sys
import pysam


def parse_ped(path):
    sex_by_sample = {}
    with open(path, "r") as handle:
        for line_number, line in enumerate(handle, start=1):
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


unphase_list = ["~{sep='\", \"' unphase_samples}"] if ~{length(unphase_samples)} > 0 else []
unphase_samples_set = set(unphase_list)

transfer_genotypes = ~{true="True" false="False" transfer_genotypes}
decrement_trv_ids = ~{true="True" false="False" decrement_trv_ids}
sex_by_sample = parse_ped("~{ped}")

base_reader = pysam.VariantFile("~{base_vcf}")
genotyped_reader = pysam.VariantFile("~{phased_vcf}") if transfer_genotypes else None

output_writer = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=base_reader.header)
base_samples = list(base_reader.header.samples)
shared_samples = set(base_samples)
if transfer_genotypes:
    shared_samples &= set(genotyped_reader.header.samples)

for record in base_reader:
    record.translate(output_writer.header)

    # Transfer GT field
    if transfer_genotypes:
        match = None
        for candidate in genotyped_reader.fetch(record.chrom, record.start, record.stop):
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

    for sample in base_samples:
        sample_data = record.samples[sample]
        sample_sex = sex_by_sample.get(sample)

        # Unphase for samples in unphase_samples
        if sample in unphase_samples_set:
            current_gt = sample_data.get("GT")
            if current_gt is not None:
                sample_data.phased = False

        # Clear format fields for females on chrY
        if record.chrom == "chrY" and sample_sex == "F":
            clear_format_fields(sample_data)
            continue

        # Make male calls hemizygous on chrX & chrY
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

    output_writer.write(record)

base_reader.close()
output_writer.close()
if genotyped_reader is not None:
    genotyped_reader.close()
CODE

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File genotyped_vcf = "~{prefix}.vcf.gz"
        File genotyped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 5 * ceil(size(base_vcf, "GB") + size(phased_vcf, "GB")) + 25,
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
