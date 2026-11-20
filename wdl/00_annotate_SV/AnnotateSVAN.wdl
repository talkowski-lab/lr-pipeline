version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow AnnotateSVAN {
    input {
        File vcf
        File vcf_idx
        Array[String] contigs
        String prefix

        Int? records_per_shard

        String type_field = "allele_type"
        String type_ins = "ins"
        String type_del = "del"
        String length_field = "allele_length"
        Int min_length = 0
        Boolean annotate_ins = true
        Boolean annotate_del = true

        File vntr_bed
        File exons_bed
        File repeats_bed
        File ref_fa
        Array[File] ref_fa_indices
        File mei_fa
        Array[File] mei_fa_indices

        String utils_docker
        String svan_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_ins
        RuntimeAttr? runtime_attr_subset_del
        RuntimeAttr? runtime_attr_shard_ins
        RuntimeAttr? runtime_attr_shard_del
        RuntimeAttr? runtime_attr_reset_filters_ins
        RuntimeAttr? runtime_attr_reset_filters_del
        RuntimeAttr? runtime_attr_generate_trf_ins
        RuntimeAttr? runtime_attr_generate_trf_del
        RuntimeAttr? runtime_attr_annotate_ins
        RuntimeAttr? runtime_attr_annotate_del
        RuntimeAttr? runtime_attr_reformat_dup_coord_ins
        RuntimeAttr? runtime_attr_extract_ins
        RuntimeAttr? runtime_attr_extract_del
        RuntimeAttr? runtime_attr_concat_ins
        RuntimeAttr? runtime_attr_concat_del
        RuntimeAttr? runtime_attr_concat_final
    }

    Boolean single_contig = length(contigs) == 1

    scatter (contig in contigs) {
        # Preprocessing
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

        # Insertions
        if (annotate_ins) {
            call Helpers.SubsetVcfByArgs as SubsetIns {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    include_args = "INFO/~{type_field}=\"~{type_ins}\" && abs(INFO/~{length_field}) >= ~{min_length}",
                    prefix = "~{prefix}.~{contig}.ins_subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_ins
            }

            if (defined(records_per_shard)) {
                call Helpers.ShardVcfByRecords as ShardIns {
                    input:
                        vcf = SubsetIns.subset_vcf,
                        vcf_idx = SubsetIns.subset_vcf_idx,
                        records_per_shard = select_first([records_per_shard]),
                        prefix = "~{prefix}.~{contig}.ins",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_shard_ins
                }
            }

            Array[File] ins_vcfs_to_process = select_first([ShardIns.shards, [SubsetIns.subset_vcf]])
            Array[File] ins_vcf_idxs_to_process = select_first([ShardIns.shard_idxs, [SubsetIns.subset_vcf_idx]])

            scatter (i in range(length(ins_vcfs_to_process))) {
                call Helpers.ResetVcfFilters as ResetIns {
                    input:
                        vcf = ins_vcfs_to_process[i],
                        vcf_idx = ins_vcf_idxs_to_process[i],
                        prefix = "~{prefix}.~{contig}.ins_shard_~{i}.reset_filters",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_reset_filters_ins
                }

                call GenerateTRF as GenerateTRFIns {
                    input:
                        vcf = ResetIns.reset_vcf,
                        vcf_idx = ResetIns.reset_vcf_idx,
                        mode = "ins",
                        prefix = "~{prefix}.~{contig}.ins_shard_~{i}.trf",
                        docker = svan_docker,
                        runtime_attr_override = runtime_attr_generate_trf_ins
                }

                call RunSvanAnnotate as SvanAnnotateIns {
                    input:
                        vcf = ResetIns.reset_vcf,
                        vcf_idx = ResetIns.reset_vcf_idx,
                        trf_output = GenerateTRFIns.trf_output,
                        vntr_bed = vntr_bed,
                        exons_bed = exons_bed,
                        repeats_bed = repeats_bed,
                        mei_fa = mei_fa,
                        mei_fa_indices = mei_fa_indices,
                        ref_fa = ref_fa,
                        ref_fa_indices = ref_fa_indices,
                        mode = "ins",
                        prefix = "~{prefix}.~{contig}.ins_shard_~{i}.svan",
                        docker = svan_docker,
                        runtime_attr_override = runtime_attr_annotate_ins
                }

                call ReformatDupCoord {
                    input:
                        vcf = SvanAnnotateIns.annotated_vcf,
                        vcf_idx = SvanAnnotateIns.annotated_vcf_idx,
                        prefix = "~{prefix}.~{contig}.ins_shard_~{i}.dup_coord_reformatted",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_reformat_dup_coord_ins
                }

                call Helpers.ExtractVcfAnnotations as ExtractIns {
                    input:
                        vcf = ReformatDupCoord.reformatted_vcf,
                        vcf_idx = ReformatDupCoord.reformatted_vcf_idx,
                        original_vcf = ResetIns.reset_vcf,
                        original_vcf_idx = ResetIns.reset_vcf_idx,
                        prefix = "~{prefix}.~{contig}.ins_shard_~{i}.annotations",
                        docker = utils_docker,
                        add_header_row = true,
                        runtime_attr_override = runtime_attr_extract_ins
                }
            }

            if (defined(records_per_shard)) {
                call Helpers.ConcatTsvs as ConcatInsShards {
                    input:
                        tsvs = ExtractIns.annotations_tsv,
                        sort_output = true,
                        preserve_header = true,
                        prefix = "~{prefix}.~{contig}.ins_annotations",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_concat_ins
                }
            }

            File final_ins_annotations = select_first([ConcatInsShards.concatenated_tsv, ExtractIns.annotations_tsv[0]])
        }

        # Deletions
        if (annotate_del) {
            call Helpers.SubsetVcfByArgs as SubsetDel {
                input:
                    vcf = contig_vcf,
                    vcf_idx = contig_vcf_idx,
                    include_args = "INFO/~{type_field}=\"~{type_del}\" && abs(INFO/~{length_field}) >= ~{min_length}",
                    prefix = "~{prefix}.~{contig}.del_subset",
                    docker = utils_docker,
                    runtime_attr_override = runtime_attr_subset_del
            }

            if (defined(records_per_shard)) {
                call Helpers.ShardVcfByRecords as ShardDel {
                    input:
                        vcf = SubsetDel.subset_vcf,
                        vcf_idx = SubsetDel.subset_vcf_idx,
                        records_per_shard = select_first([records_per_shard]),
                        prefix = "~{prefix}.~{contig}.del",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_shard_del
                }
            }

            Array[File] del_vcfs_to_process = select_first([ShardDel.shards, [SubsetDel.subset_vcf]])
            Array[File] del_vcf_idxs_to_process = select_first([ShardDel.shard_idxs, [SubsetDel.subset_vcf_idx]])

            scatter (i in range(length(del_vcfs_to_process))) {
                call Helpers.ResetVcfFilters as ResetDel {
                    input:
                        vcf = del_vcfs_to_process[i],
                        vcf_idx = del_vcf_idxs_to_process[i],
                        prefix = "~{prefix}.~{contig}.del_shard_~{i}.reset_filters",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_reset_filters_del
                }

                call GenerateTRF as GenerateTRFDel {
                    input:
                        vcf = ResetDel.reset_vcf,
                        vcf_idx = ResetDel.reset_vcf_idx,
                        mode = "del",
                        prefix = "~{prefix}.~{contig}.del_shard_~{i}.trf",
                        docker = svan_docker,
                        runtime_attr_override = runtime_attr_generate_trf_del
                }

                call RunSvanAnnotate as SvanAnnotateDel {
                    input:
                        vcf = ResetDel.reset_vcf,
                        vcf_idx = ResetDel.reset_vcf_idx,
                        trf_output = GenerateTRFDel.trf_output,
                        vntr_bed = vntr_bed,
                        exons_bed = exons_bed,
                        repeats_bed = repeats_bed,
                        mei_fa = mei_fa,
                        mei_fa_indices = mei_fa_indices,
                        ref_fa = ref_fa,
                        ref_fa_indices = ref_fa_indices,
                        mode = "del",
                        prefix = "~{prefix}.~{contig}.del_shard_~{i}.svan",
                        docker = svan_docker,
                        runtime_attr_override = runtime_attr_annotate_del
                }

                call Helpers.ExtractVcfAnnotations as ExtractDel {
                    input:
                        vcf = SvanAnnotateDel.annotated_vcf,
                        vcf_idx = SvanAnnotateDel.annotated_vcf_idx,
                        original_vcf = ResetDel.reset_vcf,
                        original_vcf_idx = ResetDel.reset_vcf_idx,
                        prefix = "~{prefix}.~{contig}.del_shard_~{i}.annotations",
                        docker = utils_docker,
                        add_header_row = true,
                        runtime_attr_override = runtime_attr_extract_del
                }
            }

            if (defined(records_per_shard)) {
                call Helpers.ConcatTsvs as ConcatDelShards {
                    input:
                        tsvs = ExtractDel.annotations_tsv,
                        sort_output = true,
                        preserve_header = true,
                        prefix = "~{prefix}.~{contig}.del_annotations",
                        docker = utils_docker,
                        runtime_attr_override = runtime_attr_concat_del
                }
            }

            File final_del_annotations = select_first([ConcatDelShards.concatenated_tsv, ExtractDel.annotations_tsv[0]])
        }
    }

    # Postprocessing
    call Helpers.ConcatAlignedTsvs {
        input:
            tsvs = flatten([select_all(final_ins_annotations), select_all(final_del_annotations)]),
            prefix = "~{prefix}.svan_annotations",
            docker = utils_docker,
            runtime_attr_override = runtime_attr_concat_final
    }

    output {
        File annotations_tsv_svan = ConcatAlignedTsvs.merged_tsv
        File annotations_header_svan = ConcatAlignedTsvs.merged_header
    }
}

task GenerateTRF {
    input {
        File vcf
        File vcf_idx
        String mode
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        gunzip -c ~{vcf} > tmp.vcf

        mkdir -p work_dir

        if [[ "~{mode}" == "ins" ]]; then
            python3 /app/SVAN/scripts/ins2fasta.py tmp.vcf work_dir
            trf work_dir/insertions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.out
        elif [[ "~{mode}" == "del" ]]; then
            python3 /app/SVAN/scripts/del2fasta.py tmp.vcf work_dir
            trf work_dir/deletions_seq.fa 2 7 7 80 10 10 500 -h -d -ngs > ~{prefix}.out
        fi
    >>>

    output {
        File trf_output = "~{prefix}.out"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 16,
        disk_gb: ceil(size(vcf, "GB")) * 3 + 20,
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

task RunSvanAnnotate {
    input {
        File vcf
        File vcf_idx
        File trf_output
        File vntr_bed
        File exons_bed
        File repeats_bed
        File mei_fa
        Array[File] mei_fa_indices
        File ref_fa
        Array[File] ref_fa_indices
        String mode
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p work_dir

        if [[ ~{vcf} == *.gz ]]; then
            gunzip -c ~{vcf} > work_dir/input.vcf
            vcf_input="work_dir/input.vcf"
        else
            vcf_input="~{vcf}"
        fi

        svan_script_name=""
        if [[ "~{mode}" == "ins" ]]; then
            svan_script_name="SVAN-INS.py"
        elif [[ "~{mode}" == "del" ]]; then
            svan_script_name="SVAN-DEL.py"
        else
            echo "Invalid mode provided."
            exit 1
        fi

        python3 /app/SVAN/$svan_script_name \
            "$vcf_input" \
            ~{trf_output} \
            ~{vntr_bed} \
            ~{exons_bed} \
            ~{repeats_bed} \
            ~{mei_fa} \
            ~{ref_fa} \
            svan_annotated \
            -o work_dir

        bcftools sort \
            --max-mem ~{select_first([runtime_attr.mem_gb, default_attr.mem_gb]) - 1}G \
            -T . \
            -Oz -o ~{prefix}.vcf.gz \
            work_dir/svan_annotated.vcf

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.vcf.gz"
        File annotated_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16,
        disk_gb: 2 * ceil(size(vcf, "GB") + size(mei_fa, "GB") + size(mei_fa_indices, "GB") + size(ref_fa, "GB")  + size(ref_fa_indices, "GB")) + 20,
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

task ReformatDupCoord {
    input {
        File vcf
        File vcf_idx
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python3 <<CODE
import pysam

def flank_to_absolute(flank_val, record_pos, alt_len):
    body, strand = flank_val.rsplit("_", 1)
    _, coord_part = body.rsplit("_", 1)
    abs_chrom, _, local_range = coord_part.split(":")
    local_start, local_end = map(int, local_range.split("-"))
    flank_start = max(1, record_pos - alt_len - 100)
    return f"{abs_chrom}:{flank_start + local_start}-{flank_start + local_end}_{strand}"

vcf_in = pysam.VariantFile("~{vcf}")
vcf_out = pysam.VariantFile("~{prefix}.vcf.gz", "w", header=vcf_in.header)
for record in vcf_in:
    dup_coord = record.info.get("DUP_COORD")
    if dup_coord is not None:
        dup_coord_str = dup_coord if isinstance(dup_coord, str) else ",".join(dup_coord)
        all_values = [v.strip() for v in dup_coord_str.split(",") if v.strip()]
        has_flank = any(v.startswith("flank_") for v in all_values)
        if has_flank or len(all_values) > 1:
            if has_flank:
                alt_len = abs(len(record.ref) - len(record.alts[0]))
            abs_values = []
            for v in all_values:
                if v.startswith("flank_"):
                    abs_values.append(flank_to_absolute(v, record.pos, alt_len))
                else:
                    abs_values.append(v)
            record.info["DUP_COORD"] = ",".join(abs_values)
    vcf_out.write(record)
vcf_in.close()
vcf_out.close()
CODE

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File reformatted_vcf = "~{prefix}.vcf.gz"
        File reformatted_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
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
