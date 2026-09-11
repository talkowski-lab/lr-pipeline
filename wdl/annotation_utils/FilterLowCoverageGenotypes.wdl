version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FilterLowCoverageGenotypes {
    input {
        File vcf
        File vcf_idx
        String prefix

        File sample_cutoffs_tsv
        File ped
        String? subset_unfilled_vcf_field
        String? subset_unfilled_vcf_value
        Boolean filter_non_ref

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_concat_tsvs
    }

    if (defined(records_per_shard)) {
        call Helpers.ShardVcfByRecords {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                records_per_shard = select_first([records_per_shard]),
                prefix = "~{prefix}.sharded",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_shard
        }
    }

    Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [vcf]])
    Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [vcf_idx]])

    scatter (i in range(length(vcfs_to_process))) {
        call FilterLowCoverageGenotypesShard {
            input:
                vcf = vcfs_to_process[i],
                vcf_idx = vcf_idxs_to_process[i],
                sample_cutoffs_tsv = sample_cutoffs_tsv,
                ped = ped,
                subset_unfilled_vcf_field = subset_unfilled_vcf_field,
                subset_unfilled_vcf_value = subset_unfilled_vcf_value,
                filter_non_ref = filter_non_ref,
                prefix = "~{prefix}.shard_~{i}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_filter
        }
    }

    if (defined(records_per_shard)) {
        call Helpers.ConcatVcfs {
            input:
                vcfs = FilterLowCoverageGenotypesShard.filtered_vcf,
                vcf_idxs = FilterLowCoverageGenotypesShard.filtered_vcf_idx,
                allow_overlaps = false,
                naive = true,
                prefix = "~{prefix}.low_coverage_filtered",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcfs
        }

        call Helpers.ConcatTsvs {
            input:
                tsvs = FilterLowCoverageGenotypesShard.filtered_genotypes_tsv,
                sort_output = false,
                preserve_header = true,
                prefix = "~{prefix}.low_coverage_filtered_genotypes",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_tsvs
        }
    }

    output {
        File filtered_vcf = select_first([ConcatVcfs.concat_vcf, FilterLowCoverageGenotypesShard.filtered_vcf[0]])
        File filtered_vcf_idx = select_first([ConcatVcfs.concat_vcf_idx, FilterLowCoverageGenotypesShard.filtered_vcf_idx[0]])
        File filtered_genotypes_tsv = select_first([ConcatTsvs.concatenated_tsv, FilterLowCoverageGenotypesShard.filtered_genotypes_tsv[0]])
    }
}

task FilterLowCoverageGenotypesShard {
    input {
        File vcf
        File vcf_idx
        File sample_cutoffs_tsv
        File ped
        String? subset_unfilled_vcf_field
        String? subset_unfilled_vcf_value
        Boolean filter_non_ref
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        # Set selected low-coverage GTs to missing while preserving all other FORMAT values.
        python3 <<'PYCODE'
import csv

import pysam


VCF = "~{vcf}"
CUTOFFS = "~{sample_cutoffs_tsv}"
PED = "~{ped}"
SUBSET_FIELD = ~{if defined(subset_unfilled_vcf_field) then "'" + subset_unfilled_vcf_field + "'" else "None"}
SUBSET_VALUE = ~{if defined(subset_unfilled_vcf_value) then "'" + subset_unfilled_vcf_value + "'" else "None"}
FILTER_NON_REF = ~{true="True" false="False" filter_non_ref}
OUTPUT_VCF = "~{prefix}.vcf.gz"
OUTPUT_TSV = "~{prefix}.filtered_genotypes.tsv"
SEX_CHROMS = {"chrX", "chrY"}


def read_cutoffs(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"sample_id", "cutoff"}
        if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
            raise ValueError("sample_cutoffs_tsv must contain sample_id and cutoff columns")
        cutoffs = {}
        for row in reader:
            sample_id = row["sample_id"]
            if not sample_id:
                raise ValueError("sample_cutoffs_tsv contains an empty sample_id")
            if sample_id in cutoffs:
                raise ValueError(f"sample_cutoffs_tsv contains duplicate sample_id: {sample_id}")
            try:
                cutoffs[sample_id] = float(row["cutoff"])
            except (TypeError, ValueError) as error:
                raise ValueError(f"Invalid cutoff for sample {sample_id}") from error
    return cutoffs


def read_ped_sexes(path):
    sexes = {}
    with open(path, "r") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"PED line {line_number} has fewer than 5 columns")
            sample_id = fields[1]
            if fields[4] == "1":
                sex = "male"
            elif fields[4] == "2":
                sex = "female"
            elif fields[4] == "0":
                sex = None
            else:
                raise ValueError(f"Sample {sample_id} has unsupported PED sex code {fields[4]}")
            if sample_id in sexes and sexes[sample_id] != sex:
                raise ValueError(f"Sample {sample_id} has conflicting PED sex entries")
            sexes[sample_id] = sex
    return sexes


def is_called(gt):
    return gt is not None and any(allele is not None for allele in gt)


def is_non_ref(gt):
    return any(allele is not None and allele > 0 for allele in gt)


def allele_counts(record):
    counts = [0] * len(record.alts)
    for sample in record.samples.values():
        gt = sample.get("GT")
        if gt is None:
            continue
        for allele in gt:
            if allele is not None and allele > 0:
                counts[allele - 1] += 1
    return ",".join(str(count) for count in counts)


def in_subset(record):
    if SUBSET_FIELD is None:
        return True
    info_val = record.info.get(SUBSET_FIELD)
    if info_val is None:
        return False
    if isinstance(info_val, (list, tuple)):
        return SUBSET_VALUE in [str(v) for v in info_val]
    return str(info_val) == SUBSET_VALUE


cutoffs = read_cutoffs(CUTOFFS)
sexes = read_ped_sexes(PED)
vcf_in = pysam.VariantFile(VCF)
vcf_samples = set(vcf_in.header.samples)

cutoff_samples = set(cutoffs)
if vcf_samples != cutoff_samples:
    missing_cutoffs = vcf_samples - cutoff_samples
    extra_cutoffs = cutoff_samples - vcf_samples
    details = []
    if missing_cutoffs:
        details.append("missing samples: " + ", ".join(sorted(missing_cutoffs)))
    if extra_cutoffs:
        details.append("unexpected samples: " + ", ".join(sorted(extra_cutoffs)))
    raise ValueError(
        "sample_cutoffs_tsv does not match VCF samples (" + "; ".join(details) + ")"
    )

missing_ped_samples = vcf_samples - set(sexes)
if missing_ped_samples:
    raise ValueError("Samples missing from PED: " + ", ".join(sorted(missing_ped_samples)))


def sample_cutoff(sample_id, chrom):
    base_cutoff = cutoffs[sample_id]
    if sexes[sample_id] == "male" and chrom in SEX_CHROMS:
        return base_cutoff / 2
    return base_cutoff


vcf_out = pysam.VariantFile(OUTPUT_VCF, "wz", header=vcf_in.header)
with open(OUTPUT_TSV, "w", newline="") as report_handle:
    report_writer = csv.writer(report_handle, delimiter="\t", lineterminator="\n")
    report_writer.writerow([
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "ID",
        "AC_before_filtering",
        "AC_after_filtering",
        "sample_count_filtered",
        "filtered_sample_ids",
    ])

    for record in vcf_in:
        ac_before = allele_counts(record)
        filtered_samples = []
        if in_subset(record) and "GT" in record.format and "DP" in record.format:
            for sample_id, sample in record.samples.items():
                gt = sample.get("GT")
                dp = sample.get("DP")
                if (not is_called(gt) or (not FILTER_NON_REF and is_non_ref(gt))
                        or dp is None or dp > sample_cutoff(sample_id, record.chrom)):
                    continue
                filtered_samples.append(sample_id)
                sample["GT"] = tuple(None for _ in gt)
                sample.phased = False

        if filtered_samples:
            report_writer.writerow([
                record.chrom,
                record.pos,
                record.ref,
                ",".join(record.alts),
                record.id or ".",
                ac_before,
                allele_counts(record),
                len(filtered_samples),
                ",".join(filtered_samples),
            ])
        vcf_out.write(record)

vcf_in.close()
vcf_out.close()
PYCODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File filtered_vcf = "~{prefix}.vcf.gz"
        File filtered_vcf_idx = "~{prefix}.vcf.gz.tbi"
        File filtered_genotypes_tsv = "~{prefix}.filtered_genotypes.tsv"
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
