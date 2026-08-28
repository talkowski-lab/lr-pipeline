version 1.0

import "Structs.wdl"
import "Helpers.wdl"

# Defragment GATK-gCNV CNVs per sample and merge

workflow DepthPreprocessing {
    input {
        Array[String]+ sample_ids
        Array[File]+ genotyped_segments_vcfs
        Array[File]+ genotyped_segments_vcf_idxs
        String prefix
        File contig_ploidy_calls_tar
        File primary_contigs_list
        File ref_fai
        File pedigree
        String batch_id

        String? chr_x
        String? chr_y

        Int gcnv_qs_cutoff
        Float? defragment_max_dist

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_gcnv_vcf_to_bed
        RuntimeAttr? runtime_attr_merge_sample
        RuntimeAttr? runtime_attr_merge_set
        RuntimeAttr? runtime_attr_make_ploidy_table
        RuntimeAttr? runtime_attr_cnv_bed_to_vcf
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (i in range(length(sample_ids))) {
        call GcnvVcfToBed {
            input:
                sample_id = sample_ids[i],
                prefix = prefix + "." + sample_ids[i],
                sample_index = i,
                vcf = genotyped_segments_vcfs[i],
                vcf_idx = genotyped_segments_vcf_idxs[i],
                contig_ploidy_calls_tar = contig_ploidy_calls_tar,
                qs_cutoff = gcnv_qs_cutoff,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_gcnv_vcf_to_bed
        }
    }

    scatter (i in range(length(sample_ids))) {
        call MergeSample as merge_sample_del {
            input:
                gcnv = GcnvVcfToBed.del_bed[i],
                prefix = prefix + "." + sample_ids[i] + ".del",
                max_dist = defragment_max_dist,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_merge_sample
        }
    }

    scatter (i in range(length(sample_ids))) {
        call MergeSample as merge_sample_dup {
            input:
                gcnv = GcnvVcfToBed.dup_bed[i],
                prefix = prefix + "." + sample_ids[i] + ".dup",
                max_dist = defragment_max_dist,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_merge_sample
        }
    }

    call MergeSet as merge_set_del {
        input:
            beds = merge_sample_del.sample_bed,
            svtype = "DEL",
            batch_id = batch_id,
            prefix = prefix + ".del",
            docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_merge_set
    }

    call MergeSet as merge_set_dup {
        input:
            beds = merge_sample_dup.sample_bed,
            svtype = "DUP",
            batch_id = batch_id,
            prefix = prefix + ".dup",
            docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_merge_set
    }

    call MakePloidyTable {
        input:
            pedigree = pedigree,
            contigs_list = primary_contigs_list,
            chr_x = chr_x,
            chr_y = chr_y,
            prefix = prefix + ".ploidy",
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_make_ploidy_table
    }

    call CNVBEDToVCF as make_del_vcf {
        input:
            bed = merge_set_del.out,
            sample_list = write_lines(sample_ids),
            contig_list = primary_contigs_list,
            ploidy_table = MakePloidyTable.ploidy_table,
            ref_fai = ref_fai,
            vid_prefix = "~{batch_id}_DEL",
            prefix = prefix + ".del",
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_cnv_bed_to_vcf
    }

    call CNVBEDToVCF as make_dup_vcf {
        input:
            bed = merge_set_dup.out,
            sample_list = write_lines(sample_ids),
            contig_list = primary_contigs_list,
            ploidy_table = MakePloidyTable.ploidy_table,
            ref_fai = ref_fai,
            vid_prefix = "~{batch_id}_DUP",
            prefix = prefix + ".dup",
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_cnv_bed_to_vcf
    }

    Array[File] concat_vcfs = [make_del_vcf.vcf, make_dup_vcf.vcf]
    RuntimeAttr default_attr_concat_vcfs = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(concat_vcfs, "GB") * 3) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr concat_attr = select_first([runtime_attr_concat_vcfs, default_attr_concat_vcfs])
    Int concat_sort_mem_mb = ceil(select_first([concat_attr.mem_gb, default_attr_concat_vcfs.mem_gb]) * 0.8 * 1024 * 1.04)

    call Helpers.ConcatVcfs as ConcatVCFs {
        input:
            vcfs = concat_vcfs,
            vcf_idxs = [make_del_vcf.vcf_idx, make_dup_vcf.vcf_idx],
            allow_overlaps = true,
            naive = false,
            sort_output = true,
            no_version = true,
            no_address = true,
            sort_mem_mb = concat_sort_mem_mb,
            prefix = prefix + ".raw_depth_cnvs",
            docker = sv_base_mini_docker,
            runtime_attr_override = concat_attr
    }

    output {
        File del_bed = merge_set_del.out
        File del_bed_idx = merge_set_del.out_idx
        File dup_bed = merge_set_dup.out
        File dup_bed_idx = merge_set_dup.out_idx
        File merged_vcf = ConcatVCFs.concat_vcf
        File merged_vcf_idx = ConcatVCFs.concat_vcf_idx
        File ploidy_table = MakePloidyTable.ploidy_table
    }
}

task GcnvVcfToBed {
    input {
        File vcf
        File vcf_idx
        File contig_ploidy_calls_tar
        Int sample_index
        String sample_id
        String prefix
        Int qs_cutoff
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        tar xzf ~{contig_ploidy_calls_tar}
        # The tar file contains one directory per sample given to GATK DetermineGermlineContigPloidy,
        # with the naming scheme SAMPLE_0 to SAMPLE_N-1, presumably in the order the samples were given
        # to the tool.
        calls_dir='~{"SAMPLE_" + sample_index}'
        expected_sample_id='~{sample_id}'
        actual_sample_id="$(cat "${calls_dir}/sample_name.txt")"
        if [[ "${expected_sample_id}" != "${actual_sample_id}" ]]; then
            printf 'Expected sample ID does not match actual sample ID\n' >&2
            printf 'Expected: %s\n' "${expected_sample_id}" >&2
            printf 'Actual: %s\n' "${actual_sample_id}" >&2
            printf 'Likely that sample order for this task differs from order given to GATK DetermineGermlineContigPloidy\n' >&2
            exit 1
        fi
        cp "${calls_dir}/contig_ploidy.tsv" contig_ploidy.tsv

        python /opt/WGD/bin/convert_gcnv.py \
            --cutoff ~{qs_cutoff} \
            contig_ploidy.tsv \
            ~{vcf} \
            ~{sample_id}

        mv ~{sample_id}.del.bed ~{prefix}.del.bed
        mv ~{sample_id}.dup.bed ~{prefix}.dup.bed
    >>>

    output {
        File del_bed = "~{prefix}.del.bed"
        File dup_bed = "~{prefix}.dup.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([vcf, vcf_idx, contig_ploidy_calls_tar], "GB") * 2) + 50,
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
        noAddress: true
    }
}

task MergeSample {
    input {
        File gcnv
        String prefix
        Float? max_dist
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        sort ~{gcnv} -k1,1V -k2,2n > ~{prefix}.bed
        bedtools merge -i ~{prefix}.bed -d 0 -c 4,5,6,7 -o distinct > ~{prefix}.merged.bed
        /opt/sv-pipeline/00_preprocessing/scripts/defragment_cnvs.py \
            --max-dist ~{if defined(max_dist) then max_dist else "0.25"} \
            ~{prefix}.merged.bed \
            ~{prefix}.merged.defrag.bed
        sort -k1,1V -k2,2n ~{prefix}.merged.defrag.bed > ~{prefix}.merged.defrag.sorted.bed
    >>>

    output {
        File sample_bed = "~{prefix}.merged.defrag.sorted.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(gcnv, "GB") * 2) + 50,
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
        noAddress: true
    }
}

task MergeSet {
    input {
        Array[File] beds
        String svtype
        String batch_id
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        cat ~{write_lines(beds)} \
            | xargs cat \
            | sort -k1,1V -k2,2n \
            | awk -v OFS="\t" -v svtype=~{svtype} -v batch=~{batch_id} '{$4=batch"_"svtype"_"NR; print}' \
            | cat <(echo -e "#chr\\tstart\\tend\\tname\\tsample\\tsvtype\\tsources") - \
            | bgzip -c > ~{prefix}.bed.gz
        tabix -p bed ~{prefix}.bed.gz
    >>>

    output {
        File out = "~{prefix}.bed.gz"
        File out_idx = "~{prefix}.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(beds, "GB") * 2) + 50,
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
        noAddress: true
    }
}

task MakePloidyTable {
    input {
        File pedigree
        File contigs_list
        String? chr_x
        String? chr_y
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python /opt/sv-pipeline/scripts/ploidy_table_from_ped.py \
            --ped '~{pedigree}' \
            --out '~{prefix}.tsv' \
            --contigs '~{contigs_list}' \
            ~{"--chr-x " + chr_x} \
            ~{"--chr-y " + chr_y}
    >>>

    output {
        File ploidy_table = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(pedigree, "GB") * 3) + 50,
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
        noAddress: true
    }
}

task CNVBEDToVCF {
    input {
        File bed
        File sample_list
        File contig_list
        File ploidy_table
        File ref_fai
        String vid_prefix
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        python /opt/sv-pipeline/scripts/convert_bed_to_gatk_vcf.py \
            --bed '~{bed}' \
            --out '~{prefix}.vcf.gz' \
            --sample '~{sample_list}' \
            --contigs '~{contig_list}' \
            --vid-prefix '~{vid_prefix}' \
            --ploidy-table '~{ploidy_table}' \
            --fai '~{ref_fai}'

        tabix '~{prefix}.vcf.gz'
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([bed, sample_list, contig_list, ploidy_table, ref_fai], "GB") * 2) + 50,
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
        noAddress: true
    }
}
