version 1.0

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow Vamos {
    input {
        File? read_bam
        File? read_bai
        Array[File]? assembly_bams
        Array[File]? assembly_bais
        File repeat_catalog_vamos
        String prefix

        String sample_id

        String vamos_docker

        RuntimeAttr? runtime_attr_run_vamos_reads
        RuntimeAttr? runtime_attr_run_vamos_contig
    }

    Boolean run_contig_mode = defined(assembly_bams)
    Boolean run_read_mode = !run_contig_mode && defined(read_bam)

    if (run_contig_mode) {
        Array[File] local_assembly_bams = select_first([assembly_bams])
        Array[File] local_assembly_bais = select_first([assembly_bais])

        scatter (i in range(length(local_assembly_bams))) {
            call RunVamos as RunVamosAssembly {
                input:
                    bam = local_assembly_bams[i],
                    bai = local_assembly_bais[i],
                    sample_id = sample_id,
                    mode = "--contig",
                    repeat_catalog_vamos = repeat_catalog_vamos,
                    prefix = "~{prefix}.assembly_~{i}",
                    docker = vamos_docker,
                    runtime_attr_override = runtime_attr_run_vamos_contig
            }
        }
    }

    if (run_read_mode) {
        call RunVamos as RunVamosReads {
            input:
                bam = select_first([read_bam]),
                bai = select_first([read_bai]),
                sample_id = sample_id,
                mode = "--read",
                repeat_catalog_vamos = repeat_catalog_vamos,
                prefix = "~{prefix}.reads",
                docker = vamos_docker,
                runtime_attr_override = runtime_attr_run_vamos_reads
        }
    }

    output {
        Array[File] vamos_assembly_vcfs = select_first([RunVamosAssembly.vamos_vcf, []])
        Array[File] vamos_assembly_vcf_idxs = select_first([RunVamosAssembly.vamos_vcf_idx, []])
        File? vamos_reads_vcf = RunVamosReads.vamos_vcf
        File? vamos_reads_vcf_idx = RunVamosReads.vamos_vcf_idx
    }
}

task RunVamos {
    input {
        File bam
        File bai
        String sample_id
        String mode
        File repeat_catalog_vamos
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        samtools quickcheck ~{bam}

        if [[ "~{repeat_catalog_vamos}" == *.gz ]]; then
            gzip -dc ~{repeat_catalog_vamos} > repeat_catalog.tsv
        else
            cp ~{repeat_catalog_vamos} repeat_catalog.tsv
        fi

        vamos \
            ~{mode} \
            -b ~{bam} \
            -r repeat_catalog.tsv \
            -s ~{sample_id} \
            -o ~{prefix}.vcf \
            -t ~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}

        bgzip -f ~{prefix}.vcf

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vamos_vcf = "~{prefix}.vcf.gz"
        File vamos_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 32,
        mem_gb: 64,
        disk_gb: 5 * ceil(size(bam, "GB") + size(repeat_catalog_vamos, "GB")) + 25,
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
