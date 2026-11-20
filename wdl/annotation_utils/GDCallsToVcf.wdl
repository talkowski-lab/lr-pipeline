version 1.0

import "../utils/Structs.wdl"
import "../utils/Helpers.wdl"

workflow GDCallsToVcf {
    input {
        Array[File] gd_output_tarballs
        Array[File] ploidy_tables
        File header_vcf
        File header_vcf_idx
        Array[String] contigs
        String prefix

        File gd_table
        File par_bed

        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_extract_gd_calls
        RuntimeAttr? runtime_attr_combine_gd_calls
        RuntimeAttr? runtime_attr_write_gd_vcf
    }

    scatter (i in range(length(gd_output_tarballs))) {
        call ExtractGDCalls {
            input:
                gd_output_tarball = gd_output_tarballs[i],
                prefix = "~{prefix}.batch_~{i}",
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_extract_gd_calls
        }
    }

    call CombineGDCalls {
        input:
            gd_calls = ExtractGDCalls.gd_calls,
            ploidy_tables = ploidy_tables,
            contigs = contigs,
            prefix = prefix,
            docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_combine_gd_calls
    }

    scatter (i in range(length(contigs))) {
        call WriteGDVcf {
            input:
                header_vcf = header_vcf,
                header_vcf_idx = header_vcf_idx,
                gd_calls = CombineGDCalls.contig_gd_calls[i],
                ploidy_table = CombineGDCalls.ploidy_table,
                contig = contigs[i],
                gd_table = gd_table,
                par_bed = par_bed,
                prefix = prefix,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_write_gd_vcf
        }
    }

    output {
        Array[File] gd_vcfs = WriteGDVcf.gd_vcf
        Array[File] gd_vcf_idxs = WriteGDVcf.gd_vcf_idx
        File gd_calls = CombineGDCalls.combined_gd_calls
    }
}

task ExtractGDCalls {
    input {
        File gd_output_tarball
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        tar -xzf ~{gd_output_tarball} --wildcards '*/call/gd_cnv_calls.tsv.gz' -O > ~{prefix}.gd_cnv_calls.tsv.gz
        gzip -t ~{prefix}.gd_cnv_calls.tsv.gz
    >>>

    output {
        File gd_calls = "~{prefix}.gd_cnv_calls.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(gd_output_tarball, "GB")) + 10,
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

task CombineGDCalls {
    input {
        Array[File] gd_calls
        Array[File] ploidy_tables
        Array[String] contigs
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        first_line() {
            python3 -c 'import gzip, sys; print(gzip.open(sys.argv[1], "rt").readline(), end="")' "$1"
        }

        first_line ~{gd_calls[0]} > header.tsv
        for f in ~{sep=" " gd_calls}; do
            if [ "$(first_line "$f")" != "$(cat header.tsv)" ]; then
                echo "ERROR: $f header differs from the first GD calls table" >&2
                exit 1
            fi
        done
        { cat header.tsv; for f in ~{sep=" " gd_calls}; do zcat "$f" | tail -n +2; done; } \
            | gzip > ~{prefix}.gd_cnv_calls.tsv.gz

        rm -f contig_files.txt
        for c in ~{sep=" " contigs}; do
            { cat header.tsv; zcat ~{prefix}.gd_cnv_calls.tsv.gz | awk -F'\t' -v c="$c" 'NR > 1 && $4 == c'; } \
                | gzip > ~{prefix}.$c.gd_cnv_calls.tsv.gz
            echo "~{prefix}.$c.gd_cnv_calls.tsv.gz" >> contig_files.txt
        done

        head -n 1 ~{ploidy_tables[0]} > ~{prefix}.ploidy.tsv
        for p in ~{sep=" " ploidy_tables}; do
            if [ "$(head -n 1 "$p")" != "$(head -n 1 ~{prefix}.ploidy.tsv)" ]; then
                echo "ERROR: $p header differs from the first ploidy table" >&2
                exit 1
            fi
            tail -n +2 "$p" >> ~{prefix}.ploidy.tsv
        done
    >>>

    output {
        File combined_gd_calls = "~{prefix}.gd_cnv_calls.tsv.gz"
        Array[File] contig_gd_calls = read_lines("contig_files.txt")
        File ploidy_table = "~{prefix}.ploidy.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 3 * ceil(size(gd_calls, "GB")) + 20,
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

task WriteGDVcf {
    input {
        File header_vcf
        File header_vcf_idx
        File gd_calls
        File ploidy_table
        String contig
        File gd_table
        File par_bed
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        awk -F'\t' -v c="~{contig}" 'NR == 1 || $1 == c' ~{gd_table} > gd_table.~{contig}.tsv

        gatk-sv-gd integrate \
            --vcf ~{header_vcf} \
            --gd-calls ~{gd_calls} \
            --gd-table gd_table.~{contig}.tsv \
            --par-bed ~{par_bed} \
            --ploidy-table ~{ploidy_table} \
            --out-vcf ~{prefix}.~{contig}.gd.vcf.gz \
            --temp-dir "$(pwd)"
    >>>

    output {
        File gd_vcf = "~{prefix}.~{contig}.gd.vcf.gz"
        File gd_vcf_idx = "~{prefix}.~{contig}.gd.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 16,
        disk_gb: 4 * ceil(size(gd_calls, "GB")) + 20,
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
