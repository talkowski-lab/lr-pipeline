version 1.0

import "../utils/Helpers.wdl"

workflow PreprocessVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix

        Array[Boolean] normalize_vcfs
        Array[String] source_tags
        Array[File] swap_sample_lists
        Array[Int] min_length_cutoffs
        Array[Int] max_length_cutoffs
        Int? records_per_shard
        Array[String] sample_ids

        File ref_fa
        File ref_fai

        String utils_docker

        RuntimeAttr? runtime_attr_validate_inputs
        RuntimeAttr? runtime_attr_swap_samples
        RuntimeAttr? runtime_attr_get_samples
        RuntimeAttr? runtime_attr_check_samples
        RuntimeAttr? runtime_attr_shard
        RuntimeAttr? runtime_attr_normalize
        RuntimeAttr? runtime_attr_subset_samples
        RuntimeAttr? runtime_attr_annotate_attributes
        RuntimeAttr? runtime_attr_add_info
        RuntimeAttr? runtime_attr_add_length_filters
        RuntimeAttr? runtime_attr_concat_shards
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_rename_variants
    }

    call ValidatePreprocessVcfsInputs {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            normalize_vcfs = normalize_vcfs,
            source_tags = source_tags,
            swap_sample_lists = swap_sample_lists,
            min_length_cutoffs = min_length_cutoffs,
            max_length_cutoffs = max_length_cutoffs,
            docker = utils_docker,
            runtime_attr_override = runtime_attr_validate_inputs
    }

    Boolean inputs_valid = ValidatePreprocessVcfsInputs.status == "success"

    if (inputs_valid) {
        if (length(swap_sample_lists) > 0) {
            scatter (vcf_index in range(length(vcfs))) {
                if (size(swap_sample_lists[vcf_index], "B") > 0) {
                    call Helpers.SwapSampleIds as SwapSampleIds {
                        input:
                            vcf = vcfs[vcf_index],
                            vcf_idx = vcf_idxs[vcf_index],
                            sample_swap_list = swap_sample_lists[vcf_index],
                            prefix = "~{prefix}.vcf_~{vcf_index}.swapped",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_swap_samples
                    }
                }

                File swapped_vcf = select_first([
                    SwapSampleIds.swapped_vcf,
                    vcfs[vcf_index]
                ])
                File swapped_vcf_idx = select_first([
                    SwapSampleIds.swapped_vcf_idx,
                    vcf_idxs[vcf_index]
                ])
            }
        }

        Array[File] swapped_vcfs = select_first([swapped_vcf, vcfs])
        Array[File] swapped_vcf_idxs = select_first([swapped_vcf_idx, vcf_idxs])

        if (length(sample_ids) > 0) {
            scatter (vcf_index in range(length(vcfs))) {
                call Helpers.SubsetVcfToSamples as SubsetInputVcfToSamples {
                    input:
                        vcf = swapped_vcfs[vcf_index],
                        vcf_idx = swapped_vcf_idxs[vcf_index],
                        samples = sample_ids,
                        prefix = "~{prefix}.vcf_~{vcf_index}.subset",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_subset_samples
                }
            }
        }

        Array[File] final_input_vcfs = select_first([
            SubsetInputVcfToSamples.subset_vcf,
            swapped_vcfs
        ])
        Array[File] final_input_vcf_idxs = select_first([
            SubsetInputVcfToSamples.subset_vcf_idx,
            swapped_vcf_idxs
        ])

        call Helpers.GetSamplesFromVcf {
            input:
                vcf = final_input_vcfs[0],
                vcf_idx = final_input_vcf_idxs[0],
                docker = utils_docker,
                runtime_attr_override = runtime_attr_get_samples
        }

        call Helpers.CheckSampleConsistency as CheckInputSamples {
            input:
                vcfs = final_input_vcfs,
                vcf_idxs = final_input_vcf_idxs,
                sample_ids = GetSamplesFromVcf.samples,
                docker = utils_docker,
                runtime_attr_override = runtime_attr_check_samples
        }

        Boolean input_samples_valid = CheckInputSamples.status == "success"

        if (input_samples_valid) {
            scatter (vcf_index in range(length(final_input_vcfs))) {
                if (defined(records_per_shard)) {
                    call Helpers.ShardVcfByRecords as ShardVcf {
                        input:
                            vcf = final_input_vcfs[vcf_index],
                            vcf_idx = final_input_vcf_idxs[vcf_index],
                            records_per_shard = select_first([records_per_shard]),
                            prefix = "~{prefix}.vcf_~{vcf_index}",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_shard
                    }
                }

                Array[File] shard_vcfs = select_first([
                    ShardVcf.shards,
                    [final_input_vcfs[vcf_index]]
                ])
                Array[File] shard_vcf_idxs = select_first([
                    ShardVcf.shard_idxs,
                    [final_input_vcf_idxs[vcf_index]]
                ])

                scatter (shard_index in range(length(shard_vcfs))) {
                    if (length(normalize_vcfs) > 0 && normalize_vcfs[vcf_index]) {
                        call Helpers.NormalizeVcf {
                            input:
                                vcf = shard_vcfs[shard_index],
                                vcf_idx = shard_vcf_idxs[shard_index],
                                ref_fa = ref_fa,
                                ref_fai = ref_fai,
                                prefix = "~{prefix}.vcf_~{vcf_index}.shard_~{shard_index}.normalized",
                                docker = utils_docker,
                                runtime_attr_override = runtime_attr_normalize
                        }
                    }

                    File normalized_vcf = select_first([
                        NormalizeVcf.normalized_vcf,
                        shard_vcfs[shard_index]
                    ])
                    File normalized_vcf_idx = select_first([
                        NormalizeVcf.normalized_vcf_idx,
                        shard_vcf_idxs[shard_index]
                    ])

                    call Helpers.AnnotateVariantAttributes {
                        input:
                            vcf = normalized_vcf,
                            vcf_idx = normalized_vcf_idx,
                            prefix = "~{prefix}.vcf_~{vcf_index}.shard_~{shard_index}.annotated",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_annotate_attributes
                    }

                    if (length(source_tags) > 0) {
                        call Helpers.AddInfo {
                            input:
                                vcf = AnnotateVariantAttributes.annotated_vcf,
                                vcf_idx = AnnotateVariantAttributes.annotated_vcf_idx,
                                tag_id = "SOURCE",
                                tag_value = source_tags[vcf_index],
                                tag_description = "Source of variant call",
                                prefix = "~{prefix}.vcf_~{vcf_index}.shard_~{shard_index}.source",
                                docker = utils_docker,
                                runtime_attr_override = runtime_attr_add_info
                        }
                    }

                    File source_annotated_vcf = select_first([
                        AddInfo.annotated_vcf,
                        AnnotateVariantAttributes.annotated_vcf
                    ])
                    File source_annotated_vcf_idx = select_first([
                        AddInfo.annotated_vcf_idx,
                        AnnotateVariantAttributes.annotated_vcf_idx
                    ])

                    if (length(min_length_cutoffs) > 0 && length(max_length_cutoffs) > 0) {
                        call AddLengthFilters as AddMinAndMaxLengthFilters {
                            input:
                                vcf = source_annotated_vcf,
                                vcf_idx = source_annotated_vcf_idx,
                                source_tag = source_tags[vcf_index],
                                min_length_cutoff = min_length_cutoffs[vcf_index],
                                max_length_cutoff = max_length_cutoffs[vcf_index],
                                source_tags = source_tags,
                                min_length_cutoffs = min_length_cutoffs,
                                max_length_cutoffs = max_length_cutoffs,
                                prefix = "~{prefix}.vcf_~{vcf_index}.shard_~{shard_index}.length_filtered",
                                docker = utils_docker,
                                runtime_attr_override = runtime_attr_add_length_filters
                        }
                    }

                    if (length(min_length_cutoffs) > 0 && length(max_length_cutoffs) == 0) {
                        call AddLengthFilters as AddMinLengthFilters {
                            input:
                                vcf = source_annotated_vcf,
                                vcf_idx = source_annotated_vcf_idx,
                                source_tag = source_tags[vcf_index],
                                min_length_cutoff = min_length_cutoffs[vcf_index],
                                source_tags = source_tags,
                                min_length_cutoffs = min_length_cutoffs,
                                max_length_cutoffs = max_length_cutoffs,
                                prefix = "~{prefix}.vcf_~{vcf_index}.shard_~{shard_index}.length_filtered",
                                docker = utils_docker,
                                runtime_attr_override = runtime_attr_add_length_filters
                        }
                    }

                    if (length(min_length_cutoffs) == 0 && length(max_length_cutoffs) > 0) {
                        call AddLengthFilters as AddMaxLengthFilters {
                            input:
                                vcf = source_annotated_vcf,
                                vcf_idx = source_annotated_vcf_idx,
                                source_tag = source_tags[vcf_index],
                                max_length_cutoff = max_length_cutoffs[vcf_index],
                                source_tags = source_tags,
                                min_length_cutoffs = min_length_cutoffs,
                                max_length_cutoffs = max_length_cutoffs,
                                prefix = "~{prefix}.vcf_~{vcf_index}.shard_~{shard_index}.length_filtered",
                                docker = utils_docker,
                                runtime_attr_override = runtime_attr_add_length_filters
                        }
                    }

                    File filtered_vcf = select_first([
                        AddMinAndMaxLengthFilters.filtered_vcf,
                        AddMinLengthFilters.filtered_vcf,
                        AddMaxLengthFilters.filtered_vcf,
                        source_annotated_vcf
                    ])
                    File filtered_vcf_idx = select_first([
                        AddMinAndMaxLengthFilters.filtered_vcf_idx,
                        AddMinLengthFilters.filtered_vcf_idx,
                        AddMaxLengthFilters.filtered_vcf_idx,
                        source_annotated_vcf_idx
                    ])
                }

                if (defined(records_per_shard)) {
                    call Helpers.ConcatVcfs as ConcatShards {
                        input:
                            vcfs = filtered_vcf,
                            vcf_idxs = filtered_vcf_idx,
                            allow_overlaps = false,
                            naive = true,
                            prefix = "~{prefix}.vcf_~{vcf_index}.concatenated",
                            docker = utils_docker,
                            runtime_attr_override = runtime_attr_concat_shards
                    }
                }

                File processed_vcf = select_first([
                    ConcatShards.concat_vcf,
                    filtered_vcf[0]
                ])
                File processed_vcf_idx = select_first([
                    ConcatShards.concat_vcf_idx,
                    filtered_vcf_idx[0]
                ])
            }
            call Helpers.CheckSampleConsistency as CheckProcessedSamples {
                input:
                    vcfs = processed_vcf,
                    vcf_idxs = processed_vcf_idx,
                    sample_ids = GetSamplesFromVcf.samples,
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_check_samples
            }

            if (length(processed_vcf) > 1) {
                call Helpers.ConcatVcfs as ConcatVcfs {
                    input:
                        vcfs = processed_vcf,
                        vcf_idxs = processed_vcf_idx,
                        allow_overlaps = true,
                        naive = false,
                        sort_output = true,
                        prefix = "~{prefix}.preprocessed",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_concat_vcfs
                }
            }

            File concatenated_vcf = select_first([
                ConcatVcfs.concat_vcf,
                processed_vcf[0]
            ])
            File concatenated_vcf_idx = select_first([
                ConcatVcfs.concat_vcf_idx,
                processed_vcf_idx[0]
            ])

            call RenameVariants {
                input:
                    vcf = concatenated_vcf,
                    vcf_idx = concatenated_vcf_idx,
                    sample_consistency_status = CheckProcessedSamples.status,
                    prefix = "~{prefix}.preprocessed",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_rename_variants
            }
        }
    }

    output {
        File preprocessed_vcf = select_first([RenameVariants.renamed_vcf])
        File preprocessed_vcf_idx = select_first([RenameVariants.renamed_vcf_idx])
    }
}

task ValidatePreprocessVcfsInputs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[Boolean] normalize_vcfs
        Array[String] source_tags
        Array[File] swap_sample_lists
        Array[Int] min_length_cutoffs
        Array[Int] max_length_cutoffs
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import json
import re

inputs = {
    "vcfs": json.load(open("~{write_json(vcfs)}")),
    "vcf_idxs": json.load(open("~{write_json(vcf_idxs)}")),
    "normalize_vcfs": json.load(open("~{write_json(normalize_vcfs)}")),
    "source_tags": json.load(open("~{write_json(source_tags)}")),
    "swap_sample_lists": json.load(open("~{write_json(swap_sample_lists)}")),
    "min_length_cutoffs": json.load(open("~{write_json(min_length_cutoffs)}")),
    "max_length_cutoffs": json.load(open("~{write_json(max_length_cutoffs)}")),
}

if not inputs["vcfs"]:
    raise ValueError("vcfs must not be empty")

expected_length = len(inputs["vcfs"])

if len(inputs["vcf_idxs"]) != expected_length:
    raise ValueError(
        f"vcf_idxs must contain one index per VCF: expected {expected_length}, got {len(inputs['vcf_idxs'])}"
    )

for name in (
    "normalize_vcfs",
    "source_tags",
    "swap_sample_lists",
    "min_length_cutoffs",
    "max_length_cutoffs",
):
    values = inputs[name]
    if values and len(values) != expected_length:
        raise ValueError(
            f"{name} is enabled, so it must contain one value per VCF: expected {expected_length}, got {len(values)}"
        )

source_tags = inputs["source_tags"]
cutoffs_enabled = bool(inputs["min_length_cutoffs"] or inputs["max_length_cutoffs"])
if cutoffs_enabled and not source_tags:
    raise ValueError("source_tags must be provided when min_length_cutoffs or max_length_cutoffs is enabled")

if source_tags:
    if len(set(source_tags)) != len(source_tags):
        raise ValueError("source_tags must be unique")
    for index, source_tag in enumerate(source_tags):
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.]*", source_tag):
            raise ValueError(f"source_tags[{index}] is not a valid VCF FILTER identifier")

for name in ("min_length_cutoffs", "max_length_cutoffs"):
    for index, cutoff in enumerate(inputs[name]):
        if cutoff < 0:
            raise ValueError(f"{name}[{index}] must be nonnegative")

if inputs["min_length_cutoffs"] and inputs["max_length_cutoffs"]:
    for index, (min_length, max_length) in enumerate(zip(
        inputs["min_length_cutoffs"],
        inputs["max_length_cutoffs"],
    )):
        if min_length > max_length:
            raise ValueError(
                f"min_length_cutoffs[{index}] must be less than or equal to max_length_cutoffs[{index}]"
            )
CODE
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

task AddLengthFilters {
    input {
        File vcf
        File vcf_idx
        String source_tag
        Int? min_length_cutoff
        Int? max_length_cutoff
        Array[String] source_tags
        Array[Int] min_length_cutoffs
        Array[Int] max_length_cutoffs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<'CODE'
import json

source_tags = json.load(open("~{write_json(source_tags)}"))
min_length_cutoffs = json.load(open("~{write_json(min_length_cutoffs)}"))
max_length_cutoffs = json.load(open("~{write_json(max_length_cutoffs)}"))

with open("length_filter_headers.txt", "w") as output:
    if min_length_cutoffs:
        for source_tag, min_length_cutoff in zip(source_tags, min_length_cutoffs):
            output.write(
                f'##FILTER=<ID=SMALL_{source_tag},Description="Allele length is below {min_length_cutoff}">\\n'
            )
    if max_length_cutoffs:
        for source_tag, max_length_cutoff in zip(source_tags, max_length_cutoffs):
            output.write(
                f'##FILTER=<ID=LARGE_{source_tag},Description="Allele length is above {max_length_cutoff}">\\n'
            )
CODE

        bcftools view -h ~{vcf} | grep "^##" > header.txt
        cat length_filter_headers.txt >> header.txt
        bcftools view -h ~{vcf} | grep "^#CHROM" >> header.txt

        bcftools reheader -h header.txt ~{vcf} | bcftools view -Oz -o reheader.vcf.gz
        input_vcf=reheader.vcf.gz

        if [ "~{defined(min_length_cutoff)}" = "true" ]; then
            bcftools filter --mode + -s SMALL_~{source_tag} -e 'abs(INFO/allele_length) < ~{select_first([min_length_cutoff, 0])}' -Oz -o min_filtered.vcf.gz "$input_vcf"
            input_vcf=min_filtered.vcf.gz
        fi

        if [ "~{defined(max_length_cutoff)}" = "true" ]; then
            bcftools filter --mode + -s LARGE_~{source_tag} -e 'abs(INFO/allele_length) > ~{select_first([max_length_cutoff, 0])}' -Oz -o max_filtered.vcf.gz "$input_vcf"
            input_vcf=max_filtered.vcf.gz
        fi

        mv "$input_vcf" ~{prefix}.vcf.gz
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

task RenameVariants {
    input {
        File vcf
        File vcf_idx
        String sample_consistency_status
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        [ "~{sample_consistency_status}" = "success" ]

        python3 <<'CODE'
from collections import defaultdict
from pysam import VariantFile

def variant_id(record):
    allele_type = record.info.get("allele_type").upper()
    if allele_type == "SNV":
        return f"{record.chrom}-{record.pos}-{record.ref}-{record.alts[0]}"
    return f"{record.chrom}-{record.pos}-{allele_type}-{abs(int(record.info.get('allele_length')))}"

input_vcf = VariantFile("~{vcf}")
counts = defaultdict(int)
for record in input_vcf:
    counts[variant_id(record)] += 1
input_vcf.close()

input_vcf = VariantFile("~{vcf}")
output_vcf = VariantFile("~{prefix}.vcf.gz", "w", header=input_vcf.header)
seen = defaultdict(int)
for record in input_vcf:
    name = variant_id(record)
    if counts[name] > 1:
        seen[name] += 1
        record.id = f"{name}_{seen[name]}"
    else:
        record.id = name
    output_vcf.write(record)
input_vcf.close()
output_vcf.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.vcf.gz"
        File renamed_vcf_idx = "~{prefix}.vcf.gz.tbi"
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
