version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateAgeMetrics {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        File age_data
        Array[Int] age_bins
        String reference_date

        String utils_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_generate_tsv
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_vcf
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        if (!single_contig) {
            call Helpers.SubsetVcfToContig {
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx,
                    contig = contig,
                    prefix = "~{prefix}.~{contig}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_vcf
            }
        }

        File contig_vcf = select_first([SubsetVcfToContig.subset_vcf, vcf])
        File contig_vcf_idx = select_first([SubsetVcfToContig.subset_vcf_idx, vcf_idx])

        if (defined(records_per_shard)) {
            call Helpers.ShardVcfByRecords {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    records_per_shard = select_first([records_per_shard]),
                    prefix = "~{prefix}.~{contig}.age",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_shard
            }
        }

        Array[File] vcfs_to_process = select_first([ShardVcfByRecords.shards, [contig_vcf]])
        Array[File] vcf_idxs_to_process = select_first([ShardVcfByRecords.shard_idxs, [contig_vcf_idx]])

        scatter (i in range(length(vcfs_to_process))) {
            call GenerateAgeAnnotationTsv {
                input:
                    vcf = vcfs_to_process[i],
                    vcf_idx = vcf_idxs_to_process[i],
                    age_data = age_data,
                    age_bins = age_bins,
                    reference_date = reference_date,
                    prefix = "~{prefix}.~{contig}.age.shard_~{i}",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_generate_tsv
            }
        }

        if (defined(records_per_shard)) {
            call Helpers.ConcatTsvs as ConcatShards {
                input:
                    tsvs = GenerateAgeAnnotationTsv.annotations_tsv,
                    sort_output = false,
                    prefix = "~{prefix}.~{contig}.age",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_concat_shards
            }
        }

        File final_annotations_tsv = select_first([ConcatShards.concatenated_tsv, GenerateAgeAnnotationTsv.annotations_tsv[0]])
    }

    if (!single_contig) {
        call Helpers.ConcatTsvs {
            input:
                tsvs = final_annotations_tsv,
                sort_output = false,
                prefix = "~{prefix}.age_annotated",
                docker = utils_docker,
                runtime_attr_override = runtime_attr_concat_vcf
        }
    }

    output {
        File annotations_tsv_age = select_first([ConcatTsvs.concatenated_tsv, final_annotations_tsv[0]])
    }
}

task GenerateAgeAnnotationTsv {
    input {
        File vcf
        File vcf_idx
        File age_data
        Array[Int] age_bins
        String reference_date
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import bisect
import csv
import json
import pysam
from datetime import datetime

with open("~{write_json(age_bins)}") as f:
    bins = json.load(f)
ref_date = datetime.strptime("~{reference_date}", "%Y-%m-%d")
prefix = "~{prefix}"

# Build age map: person_id -> age (float, in years)
age_map = {}
with open("~{age_data}") as f:
    reader = csv.DictReader(f)
    for row in reader:
        pid = str(row["person_id"]).strip()
        dob_str = str(row["date_of_birth"]).strip()
        try:
            dob = datetime.fromisoformat(dob_str.replace("+00:00", ""))
            age_map[pid] = (ref_date - dob).days / 365.25
        except (ValueError, KeyError):
            pass

bin_edges_str = "|".join(map(str, bins))
n_bins = len(bins) - 1

def get_bin_index(age):
    if age < bins[0]:
        return "smaller"
    if age >= bins[-1]:
        return "larger"
    idx = bisect.bisect_right(bins, age) - 1
    return min(idx, n_bins - 1)

vcf_in = pysam.VariantFile("~{vcf}")

with open(f"{prefix}.tsv", "w") as out:
    for record in vcf_in:
        n_alts = len(record.alts)
        het_bins = [[0] * n_bins for _ in range(n_alts)]
        het_n_smaller = [0] * n_alts
        het_n_larger = [0] * n_alts
        hom_bins = [[0] * n_bins for _ in range(n_alts)]
        hom_n_smaller = [0] * n_alts
        hom_n_larger = [0] * n_alts

        for sample_name, sample_data in record.samples.items():
            age = age_map.get(sample_name)
            if age is None:
                continue
            gt = sample_data.get("GT")
            if gt is None or None in gt:
                continue

            b = get_bin_index(age)

            for ai in range(n_alts):
                alt_allele = ai + 1
                allele_counts = gt.count(alt_allele)
                ref_counts = gt.count(0)

                is_het = allele_counts == 1 and ref_counts == len(gt) - 1
                is_hom = allele_counts == len(gt)

                if is_het:
                    if b == "smaller":
                        het_n_smaller[ai] += 1
                    elif b == "larger":
                        het_n_larger[ai] += 1
                    else:
                        het_bins[ai][b] += 1
                if is_hom:
                    if b == "smaller":
                        hom_n_smaller[ai] += 1
                    elif b == "larger":
                        hom_n_larger[ai] += 1
                    else:
                        hom_bins[ai][b] += 1

        row = [
            record.chrom,
            str(record.pos),
            record.ref,
            ",".join(record.alts),
            record.id if record.id else ".",
            ",".join("|".join(map(str, het_bins[ai])) for ai in range(n_alts)),
            ",".join(str(het_n_smaller[ai]) for ai in range(n_alts)),
            ",".join(str(het_n_larger[ai]) for ai in range(n_alts)),
            ",".join("|".join(map(str, hom_bins[ai])) for ai in range(n_alts)),
            ",".join(str(hom_n_smaller[ai]) for ai in range(n_alts)),
            ",".join(str(hom_n_larger[ai]) for ai in range(n_alts)),
        ]
        out.write("\t".join(row) + "\n")
CODE
    >>>

    output {
        File annotations_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
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
