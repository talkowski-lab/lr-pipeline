version 1.0

import "Structs.wdl"

# Defragment GATK-gCNV CNVs per sample and merge

workflow DepthPreprocessing {
    input {
        Array[String]+ sample_ids
        Array[File]+ genotyped_segments_vcfs
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
                sample_index = i,
                vcf = genotyped_segments_vcfs[i],
                contig_ploidy_calls_tar = contig_ploidy_calls_tar,
                qs_cutoff = gcnv_qs_cutoff,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_gcnv_vcf_to_bed
        }
    }

    scatter (i in range(length(sample_ids))) {
        call MergeSample as merge_sample_del {
            input:
                sample_id = sample_ids[i],
                gcnv = GcnvVcfToBed.del_bed[i],
                max_dist = defragment_max_dist,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_merge_sample
        }
    }

    scatter (i in range(length(sample_ids))) {
        call MergeSample as merge_sample_dup {
            input:
                sample_id = sample_ids[i],
                gcnv = GcnvVcfToBed.dup_bed[i],
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
            docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_merge_set
    }

    call MergeSet as merge_set_dup {
        input:
            beds = merge_sample_dup.sample_bed,
            svtype = "DUP",
            batch_id = batch_id,
            docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_merge_set
    }

    call MakePloidyTable {
        input:
            pedigree = pedigree,
            contigs_list = primary_contigs_list,
            chr_x = chr_x,
            chr_y = chr_y,
            prefix = "~{batch_id}-ploidy",
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
            prefix = "merged_del",
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
            prefix = "merged_dup",
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_cnv_bed_to_vcf
    }

    call ConcatVCFs {
        input:
            vcfs = [make_del_vcf.vcf, make_dup_vcf.vcf],
            vcf_idxs = [make_del_vcf.vcf_index, make_dup_vcf.vcf_index],
            prefix = "~{batch_id}_raw_depth_CNVs",
            docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    output {
        File del_bed = merge_set_del.out
        File del_bed_index = merge_set_del.out_idx
        File dup_bed = merge_set_dup.out
        File dup_bed_index = merge_set_dup.out_idx
        File merged_vcf = ConcatVCFs.concat_vcf
        File merged_vcf_index = ConcatVCFs.concat_vcf_index
        File ploidy_table = MakePloidyTable.ploidy_table
    }
}

task GcnvVcfToBed {
    input {
        File vcf
        File contig_ploidy_calls_tar
        Int sample_index
        String sample_id
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

        tabix ~{vcf}
        python /opt/WGD/bin/convert_gcnv.py \
            --cutoff ~{qs_cutoff} \
            contig_ploidy.tsv \
            ~{vcf} \
            ~{sample_id}
    >>>

    output {
        File del_bed = "~{sample_id}.del.bed"
        File dup_bed = "~{sample_id}.dup.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([vcf, contig_ploidy_calls_tar], "GB") * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
        String sample_id
        Float? max_dist
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        sort ~{gcnv} -k1,1V -k2,2n > ~{sample_id}.bed
        bedtools merge -i ~{sample_id}.bed -d 0 -c 4,5,6,7 -o distinct > ~{sample_id}.merged.bed
        /opt/sv-pipeline/00_preprocessing/scripts/defragment_cnvs.py \
            --max-dist ~{if defined(max_dist) then max_dist else "0.25"} \
            ~{sample_id}.merged.bed \
            ~{sample_id}.merged.defrag.bed
        sort -k1,1V -k2,2n ~{sample_id}.merged.defrag.bed > ~{sample_id}.merged.defrag.sorted.bed
    >>>

    output {
        File sample_bed = "~{sample_id}.merged.defrag.sorted.bed"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(gcnv, "GB") * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
            | bgzip -c > ~{batch_id}.~{svtype}.bed.gz;
        tabix -p bed ~{batch_id}.~{svtype}.bed.gz
    >>>

    output {
        File out = "~{batch_id}.~{svtype}.bed.gz"
        File out_idx = "~{batch_id}.~{svtype}.bed.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(beds, "GB") * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
        preemptible_tries: 3,
        max_retries: 1
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
        File vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size([bed, sample_list, contig_list, ploidy_table, ref_fai], "GB") * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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

task ConcatVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int sort_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024 * 1.04)

    command <<<
        set -euo pipefail

        bcftools concat --no-version --allow-overlaps --output-type u \
            --file-list '~{write_lines(vcfs)}' \
            | bcftools sort --max-mem '~{sort_mem_mb}' --output-type z \
                --output '~{prefix}.vcf.gz'
        tabix '~{prefix}.vcf.gz'
    >>>

    output {
        File concat_vcf = "~{prefix}.vcf.gz"
        File concat_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 2,
        mem_gb: 8,
        disk_gb: ceil(size(vcfs, "GB") * 3) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
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
