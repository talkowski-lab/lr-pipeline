version 1.0

# Clean up overlapping depth CNV records that describe the same CNV (cluster_depth_cnvs.py): cluster overlapping records
# that share most carriers, score every carrier against every record of its cluster from depth (one-copy step at both
# breakpoints), and assign each sample to its best-scoring record. Depth is read on the VM from the cohort RD matrix over
# gs:// with pysam (OAuth token from the metadata server), so the matrix is never localized.

import "../utils/Structs.wdl"

workflow ClusterDepthCnvs {
    input {
        File vcf
        File vcf_idx
        String rd_matrix
        String prefix

        String? region
        Float min_shared = 0.5
        Int flank = 2000
        Int edge = 500

        File median_coverage
        File cluster_script

        String docker

        RuntimeAttr? runtime_attr_cluster_cnvs
    }

    parameter_meta {
        vcf: "Genotyped depth CNV VCF"
        rd_matrix: "Cohort RD matrix ('#Chr Start End <samples>', bgzipped + tabix-indexed), gs:// path, read in place"
        region: "chrom or chrom:start-end to process (default whole VCF)"
        median_coverage: "Sample median coverage, 2 rows (sample IDs / medians)"
        cluster_script: "cluster_depth_cnvs.py"
    }

    call ClusterCnvs {
        input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            rd_matrix = rd_matrix,
            prefix = prefix,
            region = region,
            min_shared = min_shared,
            flank = flank,
            edge = edge,
            median_coverage = median_coverage,
            cluster_script = cluster_script,
            docker = docker,
            runtime_attr_override = runtime_attr_cluster_cnvs
    }

    output {
        File clusters = ClusterCnvs.clusters
        File scores = ClusterCnvs.scores
        File calls = ClusterCnvs.calls
    }
}

task ClusterCnvs {
    input {
        File vcf
        File vcf_idx
        String rd_matrix
        String prefix
        String? region
        Float min_shared
        Int flank
        Int edge
        File median_coverage
        File cluster_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        python3 ~{cluster_script} \
            --vcf ~{vcf} \
            --rd ~{rd_matrix} \
            --medians ~{median_coverage} \
            ~{"--region " + region} \
            --min-shared ~{min_shared} \
            --flank ~{flank} \
            --edge ~{edge} \
            --out-prefix ~{prefix}
    >>>

    output {
        File clusters = "~{prefix}.clusters.tsv"
        File scores = "~{prefix}.scores.tsv"
        File calls = "~{prefix}.calls.tsv.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 8,
        disk_gb: ceil(size(vcf, "GB")) + 20,
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
