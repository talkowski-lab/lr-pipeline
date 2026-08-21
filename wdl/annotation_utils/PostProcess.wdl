version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow PostProcess {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? shard_bin_size
        
        Boolean run_clean_vcf_header
        Boolean run_decrement_trv_ids
        Boolean run_drop_filters
        Boolean run_drop_genotypes
        Boolean run_filter_singletons
        Boolean run_flag_homopolymer_trvs
        Boolean run_normalize_ploidy
        Boolean run_prune_meis
        Boolean run_reassign_suffixes
        Boolean run_sorting
        Boolean run_transfer_genotypes
        Boolean run_unphase_samples

        Array[String] unphase_samples = []
        Array[String] drop_filters = []
        File? transfer_vcf
        File? transfer_vcf_idx
        File? ped
        
        String utils_docker

        RuntimeAttr? runtime_attr_subset_base
        RuntimeAttr? runtime_attr_subset_transfer
        RuntimeAttr? runtime_attr_create_shards
        RuntimeAttr? runtime_attr_post_process
        RuntimeAttr? runtime_attr_strip_genotypes
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat
    }

    Boolean single_contig = length(contigs) == 1
    Boolean do_transfer_genotypes = run_transfer_genotypes && defined(transfer_vcf)
    Boolean do_unphase_samples = run_unphase_samples && length(unphase_samples) > 0
    Boolean do_normalize_ploidy = run_normalize_ploidy && defined(ped)

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

            # Non-optional File so the value can cross into the shard scatter; the base placeholder is unused unless transferring
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
                        runtime_attr_override = runtime_attr_subset_base
                }

                if (do_transfer_genotypes) {
                    call Helpers.SubsetVcfToRegion as SubsetTransferShard {
                        input:
                            vcf = shard_transfer_vcf,
                            vcf_idx = shard_transfer_vcf_idx,
                            region = shard_region,
                            prefix = "~{prefix}.~{contig}.shard_~{i}.transfer",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_subset_transfer
                    }
                }

                call PostProcessTask as PostProcessShard {
                    input:
                        base_vcf = SubsetBaseShard.subset_vcf,
                        base_vcf_idx = SubsetBaseShard.subset_vcf_idx,
                        transfer_vcf = SubsetTransferShard.subset_vcf,
                        transfer_vcf_idx = SubsetTransferShard.subset_vcf_idx,
                        ped = ped,
                        unphase_samples = unphase_samples,
                        transfer_genotypes = do_transfer_genotypes,
                        unphase = do_unphase_samples,
                        normalize_ploidy = do_normalize_ploidy,
                        clean_vcf_header = run_clean_vcf_header,
                        decrement_trv_ids = run_decrement_trv_ids,
                        drop = run_drop_filters,
                        drop_filters = drop_filters,
                        filter_singletons = run_filter_singletons,
                        flag_homopolymer_trvs = run_flag_homopolymer_trvs,
                        prune_meis = run_prune_meis,
                        reassign_suffixes = run_reassign_suffixes,
                        sort_records = run_sorting,
                        prefix = "~{prefix}.~{contig}.shard_~{i}.post_processed",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_post_process
                }

                if (run_drop_genotypes) {
                    call Helpers.StripGenotypes as StripShardGenotypes {
                        input:
                            vcf = PostProcessShard.processed_vcf,
                            vcf_idx = PostProcessShard.processed_vcf_idx,
                            prefix = "~{prefix}.~{contig}.shard_~{i}.dropped",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_strip_genotypes
                    }
                }

                File shard_final_vcf = select_first([StripShardGenotypes.stripped_vcf, PostProcessShard.processed_vcf])
                File shard_final_vcf_idx = select_first([StripShardGenotypes.stripped_vcf_idx, PostProcessShard.processed_vcf_idx])
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
                    decrement_trv_ids = run_decrement_trv_ids,
                    drop = run_drop_filters,
                    drop_filters = drop_filters,
                    filter_singletons = run_filter_singletons,
                    flag_homopolymer_trvs = run_flag_homopolymer_trvs,
                    prune_meis = run_prune_meis,
                    reassign_suffixes = run_reassign_suffixes,
                    sort_records = run_sorting,
                    prefix = "~{prefix}.~{contig}.post_processed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_post_process
            }

            if (run_drop_genotypes) {
                call Helpers.StripGenotypes as StripContigGenotypes {
                    input:
                        vcf = PostProcessTask.processed_vcf,
                        vcf_idx = PostProcessTask.processed_vcf_idx,
                        prefix = "~{prefix}.~{contig}.dropped",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_strip_genotypes
                }
            }
        }

        File contig_processed_vcf = select_first([ConcatShards.concat_vcf, StripContigGenotypes.stripped_vcf, PostProcessTask.processed_vcf])
        File contig_processed_vcf_idx = select_first([ConcatShards.concat_vcf_idx, StripContigGenotypes.stripped_vcf_idx, PostProcessTask.processed_vcf_idx])
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = contig_processed_vcf,
            vcf_idxs = contig_processed_vcf_idx,
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File post_processed_vcf = ConcatVcfs.concat_vcf
        File post_processed_vcf_idx = ConcatVcfs.concat_vcf_idx
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
        Boolean decrement_trv_ids
        Boolean prune_meis
        Boolean flag_homopolymer_trvs
        Boolean sort_records
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
decrement_trv_ids = ~{true="True" false="False" decrement_trv_ids}
prune_meis = ~{true="True" false="False" prune_meis}
flag_homopolymer_trvs = ~{true="True" false="False" flag_homopolymer_trvs}
sort_records = ~{true="True" false="False" sort_records}
filter_singletons = ~{true="True" false="False" filter_singletons}
drop = ~{true="True" false="False" drop}
clean_vcf_header = ~{true="True" false="False" clean_vcf_header}
reassign_suffixes = ~{true="True" false="False" reassign_suffixes}

unphase_list = ["~{sep='\", \"' unphase_samples}"] if ~{length(unphase_samples)} > 0 else []
unphase_samples_set = set(unphase_list)
drop_filters_list = ["~{sep='\", \"' drop_filters}"] if ~{length(drop_filters)} > 0 else []
drop_filters_set = set(drop_filters_list)


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


def flush_buffer(buf, out_vcf, reassign_suffixes, sort_records):
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


sex_by_sample = parse_ped("~{default="NONE" ped}") if normalize_ploidy else {}

base_reader = pysam.VariantFile("~{base_vcf}")
transfer_reader = pysam.VariantFile("~{transfer_vcf}") if transfer_genotypes else None

header = base_reader.header.copy()
if flag_homopolymer_trvs and "HOMOPOLYMER_TRV" not in header.info:
    header.info.add("HOMOPOLYMER_TRV", 0, "Flag", "Tandem repeat call where the shortest motif has length 1.")
if filter_singletons and "SINGLE_READ_SUPPORT" not in header.filters:
    header.filters.add("SINGLE_READ_SUPPORT", None, None, "Variant supported by a single read in a single sample.")
if drop:
    for drop_filter in drop_filters_set:
        if drop_filter in header.filters:
            header.filters.remove_header(drop_filter)
if clean_vcf_header:
    header = regroup_header(header)

output_writer = pysam.VariantFile("~{prefix}.vcf.gz", "wz", header=header)
base_samples = list(base_reader.header.samples)
shared_samples = set(base_samples)
if transfer_genotypes:
    shared_samples &= set(transfer_reader.header.samples)

buffer = []
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

    # Revise MEIs outside the expected size range to indels
    if prune_meis:
        prune_record_meis(record)

    # Flag homopolymer TR variants
    if flag_homopolymer_trvs and get_scalar(record.info.get("allele_type")) == "trv" and shortest_motif_length(record) == 1:
        record.info["HOMOPOLYMER_TRV"] = True

    # Flag variants with single read support
    if filter_singletons and has_single_read_support(record):
        record.filter.add("SINGLE_READ_SUPPORT")

    if not sort_records and not reassign_suffixes:
        output_writer.write(record)
        continue

    # Buffer records sharing the exact same coordinate to sort and/or reassign ID suffixes
    if record.chrom != current_chrom or record.pos != current_pos:
        if buffer:
            flush_buffer(buffer, output_writer, reassign_suffixes, sort_records)
        buffer = [record]
        current_chrom = record.chrom
        current_pos = record.pos
    else:
        buffer.append(record)

if (sort_records or reassign_suffixes) and buffer:
    flush_buffer(buffer, output_writer, reassign_suffixes, sort_records)

base_reader.close()
output_writer.close()
if transfer_reader is not None:
    transfer_reader.close()
CODE

        tabix -p vcf -f ~{prefix}.vcf.gz
    >>>

    output {
        File processed_vcf = "~{prefix}.vcf.gz"
        File processed_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
