version 1.0

## PerSampleVariantBurden
##
## Given a set of per-chromosome VCFs from the same cohort (e.g.
## chr1.chr1.annotated.vcf.gz .. chr22.chr22.annotated.vcf.gz):
##   1. (skipped if sample_count_table is provided) for each VCF, count samples in
##      the header and count FILTER=PASS non-ref genotypes per sample, overall and
##      split into {snv, del_1_49, ins_1_49, del_50_499, ins_50_499, del_gt499,
##      ins_gt499}, each further split by INFO/REGION (US/RM/SD/SR), using the
##      VCF's own INFO/allele_type + INFO/allele_length + INFO/REGION annotations
##      (see per_sample_type_counts.sh for the exact bcftools-stats formula); then
##      sum each sample's counts genome-wide into a flat per-sample variant count
##      table (BuildSampleVariantCountTable).
##   2. rank samples ascending by their genome-wide PASS non-ref total, restricted
##      to samples in the sample/population list (RankSamples). The totals come
##      from step 1, or from sample_count_table when provided (a precomputed
##      sample vs. count TSV, e.g. the sample_variant_count_table of an earlier
##      run; header with "sample" or "sample_id" and "n_total_nonref" or "count").
##   3. cumulative union table: row i = the i-th ranked sample, i, and the number
##      of distinct PASS variant SITES where at least one of the top i
##      (smallest-burden) samples is non-ref, a set UNION across samples, split by
##      type/size. Computed in ONE pass per chromosome (UnionMinRankHistogram,
##      scattered over the VCFs): bcftools view -f PASS -S <ranked samples> -c 1
##      subsets each VCF to the ranked samples, each site is assigned the rank of
##      its best-ranked non-ref carrier, and sites are histogrammed by that rank.
##      A site is in the top-i union exactly when that rank is <= i, so summing the
##      per-chromosome histograms and taking the cumulative sum over rank gives
##      the table (CumulativeUnionCounts). This replaces the earlier one
##      bcftools pass per (chromosome x i), i.e. N x n_vcfs passes, and gives
##      identical counts (checked on a chr22 slice at i = 1..292).
##
## Per-VCF sample counts (n_samples_per_vcf) are expected to match across shards
## from the same cohort; not enforced, just surfaced for the caller to check.

workflow PerSampleVariantBurden {

    input {
        Array[File]  vcfs
        Array[File]? vcf_indices               # optional, one per vcf; accepted for the usual
                                                 # vcf+index convention (no region queries here)
        File         sample_population_list     # sample -> population (SAS/EAS/EUR/AMR/AFR)
        File?        sample_count_table         # optional sample vs. count TSV; if given, ranks
                                                 # samples from it and skips PerSampleTypeCounts
                                                 # and BuildSampleVariantCountTable
        String       output_basename
        File?        per_sample_type_counts_script            # per_sample_type_counts.sh; required
                                                               # unless sample_count_table is given
        File?        build_sample_variant_count_table_script  # build_sample_variant_count_table.py;
                                                               # required unless sample_count_table is given
        File         rank_samples_script                      # rank_samples.py
        File         union_min_rank_histogram_script          # union_min_rank_histogram.sh
        File         cumulative_union_from_min_rank_script    # cumulative_union_from_min_rank.py
        String       bcftools_docker = "quay.io/biocontainers/bcftools@sha256:badc3a0c7af72a83e5761ab0e881aa84204694bdead003b47552cb283958f78d"
        String       python_docker   = "python:3.11-slim"
        Int          mem_gb      = 8
        Int          disk_gb     = 50
        Int          preemptible = 1
    }

    if (!defined(sample_count_table)) {
        scatter (vcf in vcfs) {
            call PerSampleTypeCounts {
                input:
                    vcf         = vcf,
                    script      = select_first([per_sample_type_counts_script]),
                    docker      = bcftools_docker,
                    mem_gb      = mem_gb,
                    disk_gb     = disk_gb,
                    preemptible = preemptible
            }
        }

        call BuildSampleVariantCountTable {
            input:
                per_sample_tsvs = PerSampleTypeCounts.per_sample_counts,
                output_basename = output_basename,
                script          = select_first([build_sample_variant_count_table_script]),
                docker          = python_docker,
                mem_gb          = mem_gb,
                disk_gb         = disk_gb,
                preemptible     = preemptible
        }
    }

    Array[File] rank_inputs = if defined(sample_count_table)
                              then [select_first([sample_count_table])]
                              else select_first([PerSampleTypeCounts.per_sample_counts])

    call RankSamples {
        input:
            per_sample_tsvs        = rank_inputs,
            sample_population_list = sample_population_list,
            output_basename        = output_basename,
            script                 = rank_samples_script,
            docker                 = python_docker,
            mem_gb                 = mem_gb,
            disk_gb                = disk_gb,
            preemptible            = preemptible
    }

    # One task per chromosome VCF: a single bcftools pass covering every i.
    scatter (vcf in vcfs) {
        call UnionMinRankHistogram {
            input:
                vcf             = vcf,
                ordered_samples = RankSamples.ordered_samples,
                script          = union_min_rank_histogram_script,
                docker          = bcftools_docker,
                mem_gb          = mem_gb,
                disk_gb         = disk_gb,
                preemptible     = preemptible
        }
    }

    call CumulativeUnionCounts {
        input:
            ordered_samples = RankSamples.ordered_samples,
            histograms      = UnionMinRankHistogram.histogram,
            output_basename = output_basename,
            script          = cumulative_union_from_min_rank_script,
            docker          = python_docker,
            mem_gb          = mem_gb,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    output {
        Array[Int]?  n_samples_per_vcf             = PerSampleTypeCounts.n_samples
        Array[File]? per_sample_counts_per_vcf     = PerSampleTypeCounts.per_sample_counts
        File?        sample_variant_count_table    = BuildSampleVariantCountTable.count_table
        File         sample_rank_table             = RankSamples.sample_rank
        Array[File]  min_rank_histogram_per_vcf    = UnionMinRankHistogram.histogram
        File         ranked_sample_burden_table    = CumulativeUnionCounts.out_table
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Per-sample PASS-only non-ref variant burden across chromosome-shard VCFs: a ranked, cumulative-union-site (not summed) population-list-restricted burden table computed in one pass per chromosome, and (unless a precomputed sample vs. count table is supplied) a flat unranked genome-wide per-sample count table."
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

task RankSamples {
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

    command <<<
        set -euo pipefail
        python3 ~{script} ~{sample_population_list} ~{output_basename} ~{sep=" " per_sample_tsvs}
    >>>

    output {
        File ordered_samples = output_basename + ".ordered_samples.txt"
        File sample_rank     = output_basename + ".sample_rank.tsv"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task UnionMinRankHistogram {
    input {
        File    vcf
        File    ordered_samples
        File    script
        String  docker
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    String out_prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        bash ~{script} ~{vcf} ~{ordered_samples} ~{out_prefix}
    >>>

    output {
        File histogram = out_prefix + ".min_rank_histogram.tsv"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task CumulativeUnionCounts {
    input {
        File        ordered_samples
        Array[File] histograms
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
        python3 ~{script} --ordered-samples ~{ordered_samples} --histograms ~{sep=" " histograms} --out ~{out_name}
    >>>

    output {
        File out_table = out_name
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task BuildSampleVariantCountTable {
    input {
        Array[File] per_sample_tsvs
        String      output_basename
        File        script
        String      docker
        Int         mem_gb
        Int         disk_gb
        Int         preemptible
    }

    String out_name = output_basename + ".sample_variant_count_table.tsv"

    command <<<
        set -euo pipefail
        python3 ~{script} --shards ~{sep=" " per_sample_tsvs} --out ~{out_name}
    >>>

    output {
        File count_table = out_name
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}
