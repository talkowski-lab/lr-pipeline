version 1.0

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
    }

    scatter (i in range(length(sample_ids))) {
        call GcnvVcfToBed {
            input:
                sample_id = sample_ids[i],
                sample_index = i,
                vcf = genotyped_segments_vcfs[i],
                contig_ploidy_calls_tar = contig_ploidy_calls_tar,
                sv_pipeline_docker = sv_pipeline_docker,
                qs_cutoff = gcnv_qs_cutoff
        }
    }

    scatter (i in range(length(sample_ids))) {
        call MergeSample as merge_sample_del {
            input:
                sample_id = sample_ids[i],
                gcnv = GcnvVcfToBed.del_bed[i],
                max_dist = defragment_max_dist,
                sv_pipeline_docker = sv_pipeline_docker
        }
    }

    scatter (i in range(length(sample_ids))) {
        call MergeSample as merge_sample_dup {
            input:
                sample_id = sample_ids[i],
                gcnv = GcnvVcfToBed.dup_bed[i],
                max_dist = defragment_max_dist,
                sv_pipeline_docker = sv_pipeline_docker
        }
    }

    call MergeSet as merge_set_del {
        input:
            beds = merge_sample_del.sample_bed,
            svtype = "DEL",
            batch_id = batch_id,
            sv_base_mini_docker = sv_base_mini_docker
    }

    call MergeSet as merge_set_dup {
        input:
            beds = merge_sample_dup.sample_bed,
            svtype = "DUP",
            batch_id = batch_id,
            sv_base_mini_docker = sv_base_mini_docker
    }

    call MakePloidyTable {
        input:
            pedigree = pedigree,
            contigs_list = primary_contigs_list,
            chr_x = chr_x,
            chr_y = chr_y,
            output_prefix = "~{batch_id}-ploidy",
            sv_pipeline_docker = sv_pipeline_docker
    }

    call CNVBEDToVCF as make_del_vcf {
        input:
            bed = merge_set_del.out,
            sample_list = write_lines(sample_ids),
            contig_list = primary_contigs_list,
            ploidy_table = MakePloidyTable.ploidy_table,
            ref_fai = ref_fai,
            vid_prefix = "~{batch_id}_DEL",
            output_prefix = "merged_del",
            sv_pipeline_docker = sv_pipeline_docker
    }

    call CNVBEDToVCF as make_dup_vcf {
        input:
            bed = merge_set_dup.out,
            sample_list = write_lines(sample_ids),
            contig_list = primary_contigs_list,
            ploidy_table = MakePloidyTable.ploidy_table,
            ref_fai = ref_fai,
            vid_prefix = "~{batch_id}_DUP",
            output_prefix = "merged_dup",
            sv_pipeline_docker = sv_pipeline_docker
    }

    call ConcatVCFs {
        input:
            vcfs = [make_del_vcf.vcf, make_dup_vcf.vcf],
            vcf_idxs = [make_del_vcf.vcf_index, make_dup_vcf.vcf_index],
            output_prefix = "~{batch_id}_raw_depth_CNVs",
            sv_base_mini_docker = sv_base_mini_docker
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

        String sv_pipeline_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }

    Int default_disk_gb = ceil(size([vcf, contig_ploidy_calls_tar], "GB") * 2) + 50

    runtime {
        cpu: select_first([cpu, 1])
        memory: select_first([mem_gib, 4]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: sv_pipeline_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
        noAddress: true
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
}

task MergeSample {
    input {
        File gcnv
        String sample_id
        Float? max_dist
  
        String sv_pipeline_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }
  
    Int default_disk_gb = ceil(size(gcnv, "GB") * 2) + 50
  
    runtime {
        cpu: select_first([cpu, 1])
        memory: select_first([mem_gib, 4]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: sv_pipeline_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
        noAddress: true
    }
  
    command <<<
        set -euo pipefail
  
        sort ~{gcnv} -k1,1V -k2,2n > ~{sample_id}.bed
        bedtools merge -i ~{sample_id}.bed -d 0 -c 4,5,6,7 -o distinct > ~{sample_id}.merged.bed
        /opt/sv-pipeline/00_preprocessing/scripts/defragment_cnvs.py \
        --max-dist ~{if defined(max_dist) then max_dist else "0.25"} ~{sample_id}.merged.bed ~{sample_id}.merged.defrag.bed
        sort -k1,1V -k2,2n ~{sample_id}.merged.defrag.bed > ~{sample_id}.merged.defrag.sorted.bed
    >>>
  
    output {
        File sample_bed = "~{sample_id}.merged.defrag.sorted.bed"
    }
}

task MergeSet {
    input {
        Array[File] beds
        String svtype
        String batch_id

        String sv_base_mini_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }

    Int default_disk_gb = ceil(size(beds, "GB") * 2) + 50

    runtime {
        cpu: select_first([cpu, 1])
        memory: select_first([mem_gib, 4]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: sv_base_mini_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
        noAddress: true
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
}

task MakePloidyTable {
    input {
        File pedigree
        File contigs_list
        String? chr_x
        String? chr_y
        String output_prefix
  
        String sv_pipeline_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }
  
    Int default_disk_gb = ceil(size(pedigree, "GB") * 3) + 50
  
    runtime {
        cpu: select_first([cpu, 1])
        memory: select_first([mem_gib, 4]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: sv_pipeline_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
        noAddress: true
    }
  
    command <<<
        set -euo pipefail
  
        python /opt/sv-pipeline/scripts/ploidy_table_from_ped.py \
        --ped '~{pedigree}' \
        --out '~{output_prefix}.tsv' \
        --contigs '~{contigs_list}' \
        ~{"--chr-x " + chr_x} \
        ~{"--chr-y " + chr_y}
    >>>
  
    output {
        File ploidy_table = "~{output_prefix}.tsv"
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
        String output_prefix

        String sv_pipeline_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }

    Int default_disk_gb = ceil(size([bed, sample_list, contig_list, ploidy_table, ref_fai], "GB") * 2) + 50

    runtime {
        cpu: select_first([cpu, 1])
        memory: select_first([mem_gib, 4]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: sv_pipeline_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
        noAddress: true
    }

    command <<<
        set -euo pipefail

        python /opt/sv-pipeline/scripts/convert_bed_to_gatk_vcf.py \
        --bed '~{bed}' \
        --out '~{output_prefix}.vcf.gz' \
        --sample '~{sample_list}' \
        --contigs '~{contig_list}' \
        --vid-prefix '~{vid_prefix}' \
        --ploidy-table '~{ploidy_table}' \
        --fai '~{ref_fai}'

        tabix '~{output_prefix}.vcf.gz'
    >>>

    output {
        File vcf = "~{output_prefix}.vcf.gz"
        File vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }
}

task ConcatVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix

        String sv_base_mini_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }

    Int default_disk_gb = ceil(size(vcfs, "GB") * 3) + 50
    Int sort_mem_mb = ceil(select_first([mem_gib, 8]) * 0.8 * 1024 * 1.04)

    runtime {
        cpu: select_first([cpu, 2])
        memory: select_first([mem_gib, 8]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: sv_base_mini_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
        noAddress: true
    }

    command <<<
        set -euo pipefail

        bcftools concat --no-version --allow-overlaps --output-type u \
        --file-list '~{write_lines(vcfs)}' \
        | bcftools sort --max-mem '~{sort_mem_mb}' --output-type z \
        --output '~{output_prefix}.vcf.gz'
        tabix '~{output_prefix}.vcf.gz'
    >>>

    output {
        File concat_vcf = "~{output_prefix}.vcf.gz"
        File concat_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }
}
