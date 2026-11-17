version 1.0

# Depth edge scan for CNV discovery (edge_scan_cnv.py): sharp drops / rises of cohort-relative depth found with a step
# scan, paired into DEL / DUP candidates, breakpoints refined by shrinking each edge's window. One task per contig; the
# RD matrix is read on the VM over gs:// with pysam (OAuth token from the metadata server), never localized.

import "Structs.wdl"

workflow EdgeScanCnv {
    input {
        String rd_matrix
        Array[String] contigs
        String samples
        String prefix

        String extra_args = ""

        File median_coverage
        File ref_fai
        File scan_script

        String docker

        RuntimeAttr? runtime_attr_edge_scan
    }

    parameter_meta {
        rd_matrix: "Cohort RD matrix ('#Chr Start End <samples>', bgzipped + tabix-indexed), gs:// path, read in place"
        samples: "Comma-separated samples to scan"
        median_coverage: "Sample median coverage, 2 rows (sample IDs / medians)"
        ref_fai: "Reference FASTA index, for contig lengths"
        scan_script: "edge_scan_cnv.py"
        extra_args: "Extra arguments for scan_script, e.g. '--scale 5000 --min-z 4'"
    }

    scatter (contig in contigs) {
        call EdgeScan {
            input:
                rd_matrix = rd_matrix,
                contig = contig,
                samples = samples,
                prefix = prefix + "." + contig,
                extra_args = extra_args,
                median_coverage = median_coverage,
                ref_fai = ref_fai,
                scan_script = scan_script,
                docker = docker,
                runtime_attr_override = runtime_attr_edge_scan
        }
    }

    output {
        Array[File] candidates = EdgeScan.edges
    }
}

task EdgeScan {
    input {
        String rd_matrix
        String contig
        String samples
        String prefix
        String extra_args
        File median_coverage
        File ref_fai
        File scan_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        len=$(awk -v c=~{contig} '$1 == c {print $2}' ~{ref_fai})

        python3 ~{scan_script} \
            --rd ~{rd_matrix} \
            --medians ~{median_coverage} \
            --contig ~{contig} \
            --length ${len} \
            --samples ~{samples} \
            ~{extra_args} \
            --out ~{prefix}.edges.tsv.gz
    >>>

    output {
        File edges = "~{prefix}.edges.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: 20,
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
