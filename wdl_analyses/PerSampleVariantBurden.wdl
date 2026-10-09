version 1.0

## PerSampleVariantBurden
##
## Given a set of per-chromosome VCFs from the same cohort (e.g.
## chr1.chr1.annotated.vcf.gz .. chr22.chr22.annotated.vcf.gz), for each VCF:
##   1. count samples in the VCF header
##   2. count FILTER=PASS non-ref genotypes per sample, overall and split
##      into {snv, del_1_49, ins_1_49, del_50_499, ins_50_499, del_gt499,
##      ins_gt499}, each of those 7 further split by INFO/REGION (US/RM/SD/SR
##      genomic context), using the VCF's own INFO/allele_type +
##      INFO/allele_length + INFO/REGION annotations (see
##      per_sample_type_counts.sh for the exact bcftools-stats formula)
##
## Then, across all chromosome shards:
##   3. sum each sample's individual PASS non-ref site count genome-wide,
##      restrict to samples present in a sample/population list
##      (SAS/EAS/EUR/AMR/AFR), and rank samples ascending by that individual
##      total (RankSamples). Using that rank order, for i = 1..N compute the
##      cumulative table where row i = the i-th ranked sample, i, and the
##      number of distinct PASS variant SITES where at least one of the top
##      i (smallest-burden) samples is non-ref -- a set UNION across samples,
##      not a sum, so a site shared by multiple of the i samples is counted
##      once (CumulativeUnionCounts, one bcftools `-S <sample_list> -c 1`
##      pass per chromosome shard per i; this recomputes INFO/AC restricted
##      to the listed samples and keeps sites where that's >=1).
##   4. separately, also sum each sample's counts genome-wide into a flat
##      (non-cumulative, non-ranked, unfiltered by population) per-sample
##      variant count table covering all type/size x REGION categories.
##
## Per-VCF sample counts (n_samples_per_vcf, output below) are expected to
## match across shards from the same cohort -- not enforced here, just
## surfaced for the caller to sanity-check.
##
## Cost note: step 3's cumulative union counting does one bcftools pass per
## (chromosome shard x i), i.e. N x (number of vcfs) bcftools invocations,
## since bcftools -S needs to re-scan the file to recompute the subset AC
## regardless of subset size -- expensive if scaled to many chromosomes at
## once, though embarrassingly parallel across i.

workflow PerSampleVariantBurden {

    input {
        Array[File]  vcfs
        Array[File]? vcf_indices               # optional, one per vcf; not required by this
                                                 # workflow (no region queries), accepted for
                                                 # consistency with the usual vcf+index convention
        File         sample_population_list     # 2 cols: sample_id, population (SAS/EAS/EUR/AMR/AFR)
        String       output_basename
        File         per_sample_type_counts_script     # per_sample_type_counts.sh
        File         rank_samples_script               # rank_samples.py
        File         cumulative_union_counts_script    # cumulative_union_counts.sh
        File         concat_cumulative_counts_script   # concat_cumulative_counts.py
        File         build_sample_variant_count_table_script  # build_sample_variant_count_table.py
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

    call RankSamples {
        input:
            per_sample_tsvs         = PerSampleTypeCounts.per_sample_counts,
            sample_population_list = sample_population_list,
            output_basename         = output_basename,
            script                  = rank_samples_script,
            docker                  = python_docker,
            mem_gb                  = mem_gb,
            disk_gb                 = disk_gb,
            preemptible             = preemptible
    }

    Int n_ranked_samples = length(read_lines(RankSamples.ordered_samples))

    scatter (i in range(n_ranked_samples)) {
        call CumulativeUnionCounts {
            input:
                ordered_samples = RankSamples.ordered_samples,
                i                = i + 1,
                vcfs             = vcfs,
                script           = cumulative_union_counts_script,
                docker           = bcftools_docker,
                mem_gb           = mem_gb,
                disk_gb          = disk_gb,
                preemptible      = preemptible
        }
    }

    call ConcatCumulativeCounts {
        input:
            cumulative_rows  = CumulativeUnionCounts.row,
            output_basename  = output_basename,
            script           = concat_cumulative_counts_script,
            docker           = python_docker,
            mem_gb           = mem_gb,
            disk_gb          = disk_gb,
            preemptible      = preemptible
    }

    call BuildSampleVariantCountTable {
        input:
            per_sample_tsvs = PerSampleTypeCounts.per_sample_counts,
            output_basename = output_basename,
            script           = build_sample_variant_count_table_script,
            docker           = python_docker,
            mem_gb           = mem_gb,
            disk_gb          = disk_gb,
            preemptible      = preemptible
    }

    output {
        Array[Int]  n_samples_per_vcf         = PerSampleTypeCounts.n_samples
        Array[File] per_sample_counts_per_vcf = PerSampleTypeCounts.per_sample_counts
        File        sample_rank_table         = RankSamples.sample_rank
        File        ranked_sample_burden_table = ConcatCumulativeCounts.out_table
        File        sample_variant_count_table = BuildSampleVariantCountTable.count_table
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Per-sample PASS-only non-ref variant burden across chromosome-shard VCFs: a ranked, cumulative-union-site (not summed) population-list-restricted burden table, and a flat unranked genome-wide per-sample count table."
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

task CumulativeUnionCounts {
    input {
        File        ordered_samples
        Int         i
        Array[File] vcfs
        File        script
        String      docker
        Int         mem_gb
        Int         disk_gb
        Int         preemptible
    }

    String out_prefix = "cumulative." + i

    command <<<
        set -euo pipefail
        bash ~{script} ~{ordered_samples} ~{i} ~{out_prefix} ~{sep=" " vcfs}
    >>>

    output {
        File row = out_prefix + ".tsv"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
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
