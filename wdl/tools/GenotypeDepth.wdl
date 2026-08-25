version 1.0

workflow GenotypeDepth {
    input {
        String batch_id
        File vcf

        File training_intervals
        File median_coverage
        File rd_file
        File ref_dict
        File ploidy_table

        File contig_list

        String chr_x = "chrX"
        String chr_y = "chrY"

        String gatk_docker
        String sv_base_mini_docker
    }

    call TrainSVGenotyping {
        input:
            vcf = vcf,
            vcf_index = vcf + ".tbi",
            output_prefix = batch_id,
            training_intervals = training_intervals,
            median_coverage = median_coverage,
            chr_x = chr_x,
            chr_y = chr_y,
            rd_file = rd_file,
            rd_file_index = rd_file + ".tbi",
            ref_dict = ref_dict,
            ploidy_table = ploidy_table,
            gatk_docker = gatk_docker
    }

    scatter (contig in read_lines(contig_list)) {
        call GenotypeSVs {
            input:
                vcf = vcf,
                vcf_index = vcf + ".tbi",
                output_prefix = "~{batch_id}.genotype_batch.~{contig}",
                contig = contig,
                median_coverage = median_coverage,
                rd_file = rd_file,
                rd_file_index = rd_file + ".tbi",
                ref_dict = ref_dict,
                ploidy_table = ploidy_table,
                rd_table = TrainSVGenotyping.rd_table,
                gatk_docker = gatk_docker
        }
    }

    call ConcatVCFs {
        input:
            vcfs = GenotypeSVs.out,
            vcf_idxs = GenotypeSVs.out_index,
            output_prefix = batch_id + ".genotype_batch",
            sv_base_mini_docker = sv_base_mini_docker
    }

    output {
        File genotyped_depth_vcf = ConcatVCFs.concat_vcf
        File genotyped_depth_vcf_index = ConcatVCFs.concat_vcf_index
        File genotyping_rd_table = TrainSVGenotyping.rd_table
    }
}

task TrainSVGenotyping {
    input {
        File vcf
        File vcf_index
        File training_intervals
        File median_coverage
        File rd_file
        File rd_file_index
        String chr_x
        String chr_y
        File ref_dict
        File ploidy_table
        String output_prefix

        String gatk_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }

    parameter_meta {
        rd_file: { localization_optional: true }
    }

    Int default_disk_gb = ceil(size([vcf, rd_file], "GB") + 50)
    Int java_mem_mib = ceil(select_first([mem_gib, 16]) * 0.8 * 1024)

    runtime {
        cpu: select_first([cpu, 1])
        memory: select_first([mem_gib, 16]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: gatk_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
    }

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{java_mem_mib}M" TrainSVGenotyping \
        -XL '~{chr_x}' -XL '~{chr_y}' \
        -V '~{vcf}' \
        --training-intervals '~{training_intervals}' \
        -O '~{output_prefix}.vcf.gz' \
        --median-coverage '~{median_coverage}' \
        --rd-file '~{rd_file}' \
        --sequence-dictionary '~{ref_dict}' \
        --ploidy-table '~{ploidy_table}' \
        --output-dir ./ \
        --output-name ~{output_prefix}
    >>>

    output {
        File rd_table = "~{output_prefix}.rd_geno_params.tsv"
    }
}

task GenotypeSVs {
    input {
        File vcf
        File vcf_index
        String output_prefix
        File median_coverage
        File rd_file
        File rd_file_index
        File ref_dict
        File ploidy_table
        File rd_table
        String? contig

        String gatk_docker
        Float? mem_gib
        Int? disk_gb
        Int? cpu
        Int? boot_disk_gb
        Int? preemptible_tries
        Int? max_retries
    }

    parameter_meta {
        rd_file: { localization_optional: true }
    }

    Int default_disk_gb = ceil(size([vcf, rd_file], "GB") + 50)
    Int java_mem_mib = ceil(select_first([mem_gib, 8]) * 0.8 * 1024)

    runtime {
        cpu: select_first([cpu, 1])
        memory: select_first([mem_gib, 8]) + " GiB"
        disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([boot_disk_gb, 10])
        docker: gatk_docker
        preemptible: select_first([preemptible_tries, 3])
        maxRetries: select_first([max_retries, 1])
    }

    command <<<
        set -euo pipefail

        printf 'pe_count\tmedian_hom\tsd_het\n0\t0\t0\n' > pe_table.tsv
        printf 'sr_count\tmedian_hom\tsd_het\trare_min\trare_max\tcommon_min\tcommon_max\trare_pass\trare_fail\tcommon_pass\tcommon_fail\trare_single\trare_both\tcommon_single\tcommon_both\n' > sr_table.tsv
        printf '0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n' >> sr_table.tsv

        gatk --java-options '-Xmx~{java_mem_mib}M' PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{rd_file} \
        ~{"-L " + contig} \
        -O local.rd.txt.gz

        gatk --java-options '-Xmx~{java_mem_mib}M' GenotypeSVs \
        -V '~{vcf}' \
        -O '~{output_prefix}.vcf.gz' \
        ~{"-L " + contig} \
        --median-coverage '~{median_coverage}' \
        --rd-file local.rd.txt.gz \
        --sequence-dictionary '~{ref_dict}' \
        --ploidy-table '~{ploidy_table}' \
        --rd-table '~{rd_table}' \
        --pe-table pe_table.tsv \
        --sr-table sr_table.tsv
    >>>

    output {
        File out = "~{output_prefix}.vcf.gz"
        File out_index = "~{output_prefix}.vcf.gz.tbi"
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

        bcftools concat --no-version --allow-overlaps --output-type z --output '~{output_prefix}.vcf.gz' \
        --file-list '~{write_lines(vcfs)}'
        tabix '~{output_prefix}.vcf.gz'
    >>>

    output {
        File concat_vcf = "~{output_prefix}.vcf.gz"
        File concat_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }
}
