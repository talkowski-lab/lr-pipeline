version 1.0

# Window-level depth for CNV calling (bin_depth_windows.py): the cohort RD matrix (100 bp bins) aggregated into fixed
# windows per contig, for chosen samples, as copy number relative to each sample's genome-wide median and relative to
# the cohort median of each window. One task per contig; the matrix is read on the VM over gs:// with pysam (OAuth
# token from the metadata server), never localized.

import "../utils/Structs.wdl"

workflow BinDepthWindows {
    input {
        String rd_matrix
        Array[String] contigs
        String samples
        String prefix

        Int window = 5000
        Int? step

        File median_coverage
        File ref_fai
        File window_script

        String docker

        RuntimeAttr? runtime_attr_bin_depth
    }

    parameter_meta {
        rd_matrix: "Cohort RD matrix ('#Chr Start End <samples>', bgzipped + tabix-indexed), gs:// path, read in place"
        samples: "Comma-separated samples to report"
        median_coverage: "Sample median coverage, 2 rows (sample IDs / medians)"
        ref_fai: "Reference FASTA index, for contig lengths"
        window_script: "bin_depth_windows.py"
    }

    scatter (contig in contigs) {
        call BinDepth {
            input:
                rd_matrix = rd_matrix,
                contig = contig,
                samples = samples,
                prefix = prefix + "." + contig,
                window = window,
                step = step,
                median_coverage = median_coverage,
                ref_fai = ref_fai,
                window_script = window_script,
                docker = docker,
                runtime_attr_override = runtime_attr_bin_depth
        }
    }

    output {
        Array[File] window_tables = BinDepth.windows
    }
}

task BinDepth {
    input {
        String rd_matrix
        String contig
        String samples
        String prefix
        Int window
        Int? step
        File median_coverage
        File ref_fai
        File window_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        len=$(awk -v c=~{contig} '$1 == c {print $2}' ~{ref_fai})

        python3 ~{window_script} \
            --rd ~{rd_matrix} \
            --medians ~{median_coverage} \
            --contig ~{contig} \
            --length ${len} \
            --samples ~{samples} \
            --window ~{window} \
            ~{"--step " + step} \
            --out ~{prefix}.windows.tsv.gz
    >>>

    output {
        File windows = "~{prefix}.windows.tsv.gz"
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
