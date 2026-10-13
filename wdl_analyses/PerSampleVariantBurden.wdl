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
##   3. cumulative union table at checkpoint sample sizes (MakeCutoffs): 1..10 by
##      1, 12..20 by 2, 25..100 by 5, 110..1,000 by 10, 1,100..10,000 by 100, then
##      by 1,000, plus the last rank N (e.g. ..., 1,200, 1,300, 1,319). One task per
##      cutoff i (CumulativeUnionCounts, scattered over the cutoffs) counts the
##      distinct PASS variant SITES where at least one of the top i (smallest-
##      burden) samples is non-ref, a set UNION across samples, split by type/size:
##      bcftools view -f PASS -S <top-i samples> -c 1 on every chromosome VCF, with
##      the VCFs processed in parallel across the task's CPUs (cumulative_cpu).
##      ConcatCumulativeCounts stacks the rows in cutoff order.
##
## Only samples in sample_population_list are analyzed: PerSampleTypeCounts reads
## only them (bcftools view -S), so the per-sample count table covers only them;
## RankSamples ranks only those also present in the VCFs / sample_count_table; and
## the cumulative union is over that ranked set. To run on a subset of the VCF's
## samples, pass a population list containing just that subset.
##
## Per-VCF sample counts (n_samples_per_vcf, the VCF header's sample count) are
## expected to match across shards from the same cohort; not enforced, just
## surfaced for the caller to check.

workflow PerSampleVariantBurden {

    input {
        Array[File]  vcfs
        Array[File]? vcf_indices               # optional, one per vcf; accepted for the usual
                                                 # vcf+index convention (no region queries here)
        File         sample_population_list     # sample -> population (SAS/EAS/EUR/AMR/AFR); only
                                                 # these samples are analyzed
        File?        sample_count_table         # optional sample vs. count TSV; if given, ranks
                                                 # samples from it and skips PerSampleTypeCounts
                                                 # and BuildSampleVariantCountTable
        String       output_basename
        File?        per_sample_type_counts_script            # per_sample_type_counts.sh; required
                                                               # unless sample_count_table is given
        File?        build_sample_variant_count_table_script  # build_sample_variant_count_table.py;
                                                               # required unless sample_count_table is given
        File         rank_samples_script                      # rank_samples.py
        File         make_cutoffs_script                      # make_cutoffs.py
        File         cumulative_union_counts_script           # cumulative_union_counts.sh
        File         concat_cumulative_counts_script          # concat_cumulative_counts.py
        String       bcftools_docker = "quay.io/biocontainers/bcftools@sha256:badc3a0c7af72a83e5761ab0e881aa84204694bdead003b47552cb283958f78d"
        String       python_docker   = "python:3.11-slim"
        Int          mem_gb      = 8
        Int          disk_gb     = 50
        Int          preemptible = 1
        Int          cumulative_cpu = 4     # per-cutoff task: chromosome VCFs counted in parallel
    }

    if (!defined(sample_count_table)) {
        scatter (vcf in vcfs) {
            call PerSampleTypeCounts {
                input:
                    vcf         = vcf,
                    sample_population_list = sample_population_list,
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

    call MakeCutoffs {
        input:
            ordered_samples = RankSamples.ordered_samples,
            script          = make_cutoffs_script,
            docker          = python_docker,
            preemptible     = preemptible
    }

    # One task per cutoff i: union of non-ref PASS sites over the top-i ranked samples.
    scatter (i in MakeCutoffs.cutoffs) {
        call CumulativeUnionCounts {
            input:
                ordered_samples = RankSamples.ordered_samples,
                i               = i,
                vcfs            = vcfs,
                script          = cumulative_union_counts_script,
                docker          = bcftools_docker,
                cpu             = cumulative_cpu,
                mem_gb          = mem_gb,
                preemptible     = preemptible
        }
    }

    call ConcatCumulativeCounts {
        input:
            cumulative_rows = CumulativeUnionCounts.row,
            output_basename = output_basename,
            script          = concat_cumulative_counts_script,
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
        Array[Int]   cumulative_cutoffs            = MakeCutoffs.cutoffs
        File         ranked_sample_burden_table    = ConcatCumulativeCounts.out_table
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Per-sample PASS-only non-ref variant burden across chromosome-shard VCFs: a ranked, cumulative-union-site (not summed) population-list-restricted burden table at checkpoint sample sizes (one task per cutoff), and (unless a precomputed sample vs. count table is supplied) a flat unranked genome-wide per-sample count table, all restricted to the samples in sample_population_list."
    }
}

task PerSampleTypeCounts {
    input {
        File    vcf
        File    sample_population_list
        File    script
        String  docker
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    String out_prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        bash ~{script} ~{vcf} ~{out_prefix} ~{sample_population_list}
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

task MakeCutoffs {
    input {
        File    ordered_samples
        File    script
        String  docker
        Int     preemptible
    }

    command <<<
        set -euo pipefail
        python3 ~{script} ~{ordered_samples} cutoffs.txt
    >>>

    output {
        Array[Int] cutoffs = read_lines("cutoffs.txt")
    }

    runtime {
        docker:      docker
        memory:      "2 GB"
        cpu:         1
        disks:       "local-disk 10 HDD"
        preemptible: preemptible
    }
}

task CumulativeUnionCounts {
    input {
        File        ordered_samples
        Int         i
        Array[File] vcfs
        File        script
        String      docker
        Int         cpu
        Int         mem_gb
        Int         preemptible
    }

    String out_prefix = "cumulative." + i
    Int disk_gb = ceil(size(vcfs, "GB")) + 20

    command <<<
        set -euo pipefail
        N_PARALLEL=~{cpu} bash ~{script} ~{ordered_samples} ~{i} ~{out_prefix} ~{sep=" " vcfs}
    >>>

    output {
        File row = out_prefix + ".tsv"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         cpu
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task ConcatCumulativeCounts {
    input {
        Array[File] cumulative_rows
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
        python3 ~{script} ~{out_name} ~{sep=" " cumulative_rows}
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
