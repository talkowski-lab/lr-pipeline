version 1.0

## PerSampleVariantBurden
##
## Given a set of per-chromosome VCFs from the same cohort (e.g.
## chr1.chr1.annotated.vcf.gz .. chr22.chr22.annotated.vcf.gz), for each VCF:
##   1. count samples in the VCF header
##   2. count non-ref genotypes per sample, overall and split into
##      {snv, del_1_49, ins_1_49, del_50_499, ins_50_499, del_gt499, ins_gt499}
##      using the VCF's own INFO/allele_type + INFO/allele_length annotations
##      (see per_sample_type_counts.sh for the exact bcftools-stats formula)
##
## Then, across all chromosome shards:
##   3. sum each sample's counts genome-wide, restrict to samples present in
##      a 2-column sample/population list (SAS/EAS/EUR/AMR/AFR), sort samples
##      ascending by total non-ref count, and emit a cumulative table where
##      row i = sample i, count of samples so far, and cumulative totals
##      (overall + per type/size bucket) for that sample and all samples
##      ranked below it.
##
## Per-VCF sample counts (n_samples_per_vcf, output below) are expected to
## match across shards from the same cohort -- not enforced here, just
## surfaced for the caller to sanity-check.

workflow PerSampleVariantBurden {

    input {
        Array[File]  vcfs
        Array[File]? vcf_indices               # optional, one per vcf; not required by this
                                                 # workflow (no region queries), accepted for
                                                 # consistency with the usual vcf+index convention
        File         sample_population_list     # 2 cols: sample_id, population (SAS/EAS/EUR/AMR/AFR)
        String       output_basename
        File         per_sample_type_counts_script  # per_sample_type_counts.sh
        File         combine_and_rank_script         # combine_and_rank_samples.py
        String       bcftools_docker = "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
        String       python_docker   = "python:3.11-slim"
        Int          mem_gb      = 8
        Int          disk_gb     = 50
        Int          preemptible = 1
    }

    scatter (vcf in vcfs) {
        call PerSampleTypeCounts {
            input:
                vcf         = vcf,
                script       = per_sample_type_counts_script,
                docker       = bcftools_docker,
                mem_gb       = mem_gb,
                disk_gb      = disk_gb,
                preemptible  = preemptible
        }
    }

    call CombineAndRankSamples {
        input:
            per_sample_tsvs         = PerSampleTypeCounts.per_sample_counts,
            sample_population_list = sample_population_list,
            output_basename         = output_basename,
            script                  = combine_and_rank_script,
            docker                  = python_docker,
            mem_gb                  = mem_gb,
            disk_gb                 = disk_gb,
            preemptible             = preemptible
    }

    output {
        Array[Int]  n_samples_per_vcf        = PerSampleTypeCounts.n_samples
        Array[File] per_sample_counts_per_vcf = PerSampleTypeCounts.per_sample_counts
        File        ranked_sample_burden_table = CombineAndRankSamples.ranked_table
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Per-sample non-ref variant burden (by type/size) across chromosome-shard VCFs, ranked and cumulated ascending by total burden, restricted to a given sample/population list."
    }
}

task PerSampleTypeCounts {
    input {
        File    vcf
        File    script
        String  docker
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    String out_prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        bash ~{script} ~{vcf} ~{out_prefix}
    >>>

    output {
        File per_sample_counts = out_prefix + ".per_sample_counts.tsv"
        Int  n_samples         = read_int(out_prefix + ".n_samples.txt")
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task CombineAndRankSamples {
    input {
        Array[File] per_sample_tsvs
        File        sample_population_list
        String      output_basename
        File        script
        String      docker
        Int         mem_gb
        Int         disk_gb
        Int         preemptible
    }

    String out_name = output_basename + ".ranked_sample_burden.tsv"

    command <<<
        set -euo pipefail
        python3 ~{script} ~{sample_population_list} ~{out_name} ~{sep=" " per_sample_tsvs}
    >>>

    output {
        File ranked_table = out_name
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}
