version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow FilterLowCoverageGenotypes {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        File sample_cutoffs_tsv
        String prefix

        Int? records_per_shard

        String utils_docker

        RuntimeAttr? runtime_attr_subset
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_concat_tsvs
    }

    scatter (contig in contigs) {
        call Helpers.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_subset
        }

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = SubsetVcfToContig.subset_vcf,
                    vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.sharded",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] shard_vcfs = select_first([ShardVcfByRecords.shards, [SubsetVcfToContig.subset_vcf]])
        Array[File] shard_vcf_idxs = select_first([ShardVcfByRecords.shard_idxs, [SubsetVcfToContig.subset_vcf_idx]])

        scatter (i in range(length(shard_vcfs))) {
            call FilterLowCoverageGenotypesShard {
                input:
                    vcf = shard_vcfs[i],
                    vcf_idx = shard_vcf_idxs[i],
                    sample_cutoffs_tsv = sample_cutoffs_tsv,
                    prefix = "~{prefix}.~{contig}.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_filter
            }
        }
    }

    call Helpers.ConcatVcfs {
        input:
            vcfs = flatten(FilterLowCoverageGenotypesShard.filtered_vcf),
            vcf_idxs = flatten(FilterLowCoverageGenotypesShard.filtered_vcf_idx),
            allow_overlaps = false,
            naive = true,
            prefix = "~{prefix}.low_coverage_filtered",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    call Helpers.ConcatTsvs {
        input:
            tsvs = flatten(FilterLowCoverageGenotypesShard.filtered_genotypes_tsv),
            sort_output = false,
            preserve_header = true,
            prefix = "~{prefix}.low_coverage_filtered_genotypes",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_tsvs
    }

    output {
        File filtered_vcf = ConcatVcfs.concat_vcf
        File filtered_vcf_idx = ConcatVcfs.concat_vcf_idx
        File filtered_genotypes_tsv = ConcatTsvs.concatenated_tsv
    }
}

task FilterLowCoverageGenotypesShard {
    input {
        File vcf
        File vcf_idx
        File sample_cutoffs_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'PYCODE'
import csv

import pysam


VCF = "~{vcf}"
CUTOFFS = "~{sample_cutoffs_tsv}"
OUTPUT_VCF = "~{prefix}.vcf.gz"
OUTPUT_TSV = "~{prefix}.filtered_genotypes.tsv"


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


def is_nonref(gt):
    return gt is not None and any(allele is not None and allele > 0 for allele in gt)


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


cutoffs = read_cutoffs(CUTOFFS)
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
        if "GT" in record.format and "DP" in record.format:
            for sample_id, sample in record.samples.items():
                gt = sample.get("GT")
                dp = sample.get("DP")
                if not is_nonref(gt) or dp is None or dp > cutoffs[sample_id]:
                    continue
                filtered_samples.append(sample_id)
                for field in record.format.keys():
                    value = sample.get(field)
                    sample[field] = (
                        tuple(None for _ in value)
                        if isinstance(value, tuple)
                        else None
                    )
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
