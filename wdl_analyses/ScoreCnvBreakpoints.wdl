version 1.0

# Breakpoint check v2 for depth-genotyped CNV carriers (score_cnv_breakpoints.py): SD-aware, size-scaled breakpoint
# windows, the sharpness model plus a minimum copy-number step at both ends. One task per region, scattered. Depth is
# read on the VM from the cohort RD matrix over gs:// with pysam (OAuth token from the metadata server), so the matrix
# is never localized.

import "Structs.wdl"

workflow ScoreCnvBreakpoints {
    input {
        File vcf
        File vcf_idx
        String rd_matrix
        Array[String] regions
        String prefix

        String extra_args = ""

        File median_coverage
        File segdup_bed
        File sharpness_model
        File score_script

        String docker

        RuntimeAttr? runtime_attr_score_breakpoints
    }

    parameter_meta {
        vcf: "Genotyped depth CNV VCF"
        rd_matrix: "Cohort RD matrix ('#Chr Start End <samples>', bgzipped + tabix-indexed), gs:// path, read in place"
        regions: "Regions to score, one task each (chrom or chrom:start-end)"
        extra_args: "Extra arguments for score_script, e.g. '--min-step 0.4'"
        median_coverage: "Sample median coverage, 2 rows (sample IDs / medians)"
        sharpness_model: "Sharpness model JSON (score_sharpness_gbm.py format)"
        score_script: "score_cnv_breakpoints.py"
    }

    scatter (region in regions) {
        call ScoreBreakpoints {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                rd_matrix = rd_matrix,
                region = region,
                prefix = prefix + "." + sub(region, "[:-]", "_"),
                extra_args = extra_args,
                median_coverage = median_coverage,
                segdup_bed = segdup_bed,
                sharpness_model = sharpness_model,
                score_script = score_script,
                docker = docker,
                runtime_attr_override = runtime_attr_score_breakpoints
        }
    }

    output {
        Array[File] scores = ScoreBreakpoints.scores
    }
}

task ScoreBreakpoints {
    input {
        File vcf
        File vcf_idx
        String rd_matrix
        String region
        String prefix
        String extra_args
        File median_coverage
        File segdup_bed
        File sharpness_model
        File score_script
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        export GCS_OAUTH_TOKEN=$(curl -sS -H "Metadata-Flavor: Google" \
            http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token \
            | python3 -c 'import json,sys; print(json.load(sys.stdin)["access_token"])')

        python3 ~{score_script} \
            --vcf ~{vcf} \
            --rd ~{rd_matrix} \
            --medians ~{median_coverage} \
            --segdup ~{segdup_bed} \
            --model ~{sharpness_model} \
            --region ~{region} \
            ~{extra_args} \
            --out ~{prefix}.bp_check_v2.tsv.gz
    >>>

    output {
        File scores = "~{prefix}.bp_check_v2.tsv.gz"
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
