version 1.0

import "Structs.wdl"
import "Helpers.wdl"

workflow GenotypeDepth {
    input {
        String prefix
        File vcf
        File vcf_idx

        File training_intervals
        File median_coverage
        File rd_file
        File rd_file_idx
        File ref_dict
        File ploidy_table

        File contig_list
        File? contig_subset_list

        String chr_x = "chrX"
        String chr_y = "chrY"

        String gatk_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_train_sv_genotyping
        RuntimeAttr? runtime_attr_genotype_svs
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    call TrainSVGenotyping {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            prefix = prefix + ".train_sv_genotyping",
            training_intervals = training_intervals,
            median_coverage = median_coverage,
            chr_x = chr_x,
            chr_y = chr_y,
            rd_file = rd_file,
            rd_file_idx = rd_file_idx,
            ref_dict = ref_dict,
            ploidy_table = ploidy_table,
            docker = gatk_docker,
            runtime_attr_override = runtime_attr_train_sv_genotyping
    }

    Array[String] contigs = read_lines(select_first([contig_subset_list, contig_list]))
    scatter (contig in contigs) {
        call GenotypeSVs {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                prefix = "~{prefix}.genotype_batch.~{contig}",
                contig = contig,
                median_coverage = median_coverage,
                rd_file = rd_file,
                rd_file_idx = rd_file_idx,
                ref_dict = ref_dict,
                ploidy_table = ploidy_table,
                rd_table = TrainSVGenotyping.rd_table,
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_genotype_svs
        }
    }

    RuntimeAttr default_attr_concat_vcfs = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: ceil(size(GenotypeSVs.out, "GB") * 3) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr concat_attr = select_first([runtime_attr_concat_vcfs, default_attr_concat_vcfs])

    call Helpers.ConcatVcfs as ConcatVCFs {
        input:
            vcfs = GenotypeSVs.out,
            vcf_idxs = GenotypeSVs.out_idx,
            allow_overlaps = true,
            naive = false,
            no_version = true,
            no_address = true,
            prefix = prefix + ".genotype_batch",
            docker = sv_base_mini_docker,
            runtime_attr_override = concat_attr
    }

    output {
        File genotyped_depth_vcf = ConcatVCFs.concat_vcf
        File genotyped_depth_vcf_idx = ConcatVCFs.concat_vcf_idx
        File genotyping_rd_table = TrainSVGenotyping.rd_table
    }
}

task TrainSVGenotyping {
    input {
        File vcf
        File vcf_idx
        File training_intervals
        File median_coverage
        File rd_file
        File rd_file_idx
        String chr_x
        String chr_y
        File ref_dict
        File ploidy_table
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        rd_file: { localization_optional: true }
    }

    Int java_mem_mib = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{java_mem_mib}M" TrainSVGenotyping \
            -XL '~{chr_x}' \
            -XL '~{chr_y}' \
            -V '~{vcf}' \
            --training-intervals '~{training_intervals}' \
            -O '~{prefix}.vcf.gz' \
            --median-coverage '~{median_coverage}' \
            --rd-file '~{rd_file}' \
            --sequence-dictionary '~{ref_dict}' \
            --ploidy-table '~{ploidy_table}' \
            --output-dir ./ \
            --output-name ~{prefix}
    >>>

    output {
        File rd_table = "~{prefix}.rd_geno_params.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: ceil(size([vcf, rd_file], "GB") + 50),
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

task GenotypeSVs {
    input {
        File vcf
        File vcf_idx
        String prefix
        File median_coverage
        File rd_file
        File rd_file_idx
        File ref_dict
        File ploidy_table
        File rd_table
        String? contig
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        rd_file: { localization_optional: true }
    }

    Int java_mem_mib = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

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
            -O '~{prefix}.vcf.gz' \
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
        File out = "~{prefix}.vcf.gz"
        File out_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size([vcf, rd_file], "GB") + 50),
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
