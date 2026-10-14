version 1.0

## PerSampleVariantBurden
##
## Given a set of per-chromosome VCFs from the same cohort (e.g.
## chr1.chr1.annotated.vcf.gz .. chr22.chr22.annotated.vcf.gz):
##   1. (always run, in parallel with 2-3) for each VCF, count samples in
##      the header and count FILTER=PASS non-ref genotypes per sample, overall and
##      split into {snv, del_1_49, ins_1_49, del_50_499, ins_50_499, del_gt499,
##      ins_gt499}, each further split by INFO/REGION (US/RM/SD/SR), using the
##      VCF's own INFO/allele_type + INFO/allele_length + INFO/REGION annotations
##      (see per_sample_type_counts.sh for the exact bcftools-stats formula); then
##      sum each sample's counts genome-wide into a flat per-sample variant count
##      table (BuildSampleVariantCountTable).
##   2. sample order for the cumulative union. If ranked_samples is given (one
##      sample ID per line, rank 1 first), that order is used as is and no ranking
##      is computed. Otherwise samples are ranked ascending by their PASS non-ref
##      total, restricted to samples in the sample/population list (RankSamples),
##      ONCE, from a single small contig: chr22 if one of the VCFs is chr22 (by
##      filename), else the shortest autosome provided (GRCh38 lengths), counted
##      with a quick total-only pass (RankingContigCounts). The same order is used
##      for all contigs. The per-sample counts of step 1 do not feed the ranking,
##      so they run in parallel with the ranking and the cumulative union.
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
## split_by_allele_type (default true) splits every count by type/size using
## INFO/allele_type + INFO/allele_length (and by INFO/REGION when that is defined).
## Both tasks check the VCF header first and fail with a clear message if those
## tags are missing. Set it to false for VCFs without them: only total PASS non-ref
## counts are produced (per sample, and for the cumulative union).
## allele_type_tag (default "allele_type") names the INFO tag used as the variant
## type in place of INFO/allele_type; its values must still be snv/del/ins.
##
## Only samples in sample_population_list are analyzed: PerSampleTypeCounts reads
## only them (bcftools view -S), so the per-sample count table covers only them;
## RankSamples ranks only those also present in the VCFs; and
## the cumulative union is over that ranked set. To run on a subset of the VCF's
## samples, pass a population list containing just that subset.
## A ranked_samples list is used as given (not filtered by the population list);
## every sample in it must be present in each VCF.
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
        File?        ranked_samples             # optional ranked sample list (one ID per line, rank 1
                                                 # first); if given, CumulativeUnionCounts uses this
                                                 # order and the chr22 ranking is skipped
        String       output_basename
        # Scripts: give scripts_dir (a folder, e.g. gs://bucket/path, holding all six
        # scripts under their standard names below), and/or individual script inputs,
        # which override scripts_dir for that script.
        String?      scripts_dir
        File?        per_sample_type_counts_script            # per_sample_type_counts.sh
        File?        build_sample_variant_count_table_script  # build_sample_variant_count_table.py
        File?        rank_samples_script                      # rank_samples.py
        File?        make_cutoffs_script                      # make_cutoffs.py
        File?        cumulative_union_counts_script           # cumulative_union_counts.sh
        File?        concat_cumulative_counts_script          # concat_cumulative_counts.py
        String       bcftools_docker = "quay.io/biocontainers/bcftools@sha256:badc3a0c7af72a83e5761ab0e881aa84204694bdead003b47552cb283958f78d"
        String       python_docker   = "python:3.11-slim"
        Int          mem_gb      = 8
        Int          disk_gb     = 50
        Int          preemptible = 1
        Int          cumulative_cpu = 4     # per-cutoff task: chromosome VCFs counted in parallel
        Boolean      split_by_allele_type = true   # false if the VCFs lack INFO/allele_type + allele_length
        String       allele_type_tag = "allele_type"   # INFO tag holding the type (values snv/del/ins)
    }

    # Resolve each script: explicit input if given, else <scripts_dir>/<standard name>.
    String scripts_base = sub(select_first([scripts_dir, ""]), "/+$", "")
    File per_sample_type_counts_sh = if defined(per_sample_type_counts_script)
        then select_first([per_sample_type_counts_script])
        else scripts_base + "/per_sample_type_counts.sh"
    File build_sample_variant_count_table_py = if defined(build_sample_variant_count_table_script)
        then select_first([build_sample_variant_count_table_script])
        else scripts_base + "/build_sample_variant_count_table.py"
    File rank_samples_py = if defined(rank_samples_script)
        then select_first([rank_samples_script])
        else scripts_base + "/rank_samples.py"
    File make_cutoffs_py = if defined(make_cutoffs_script)
        then select_first([make_cutoffs_script])
        else scripts_base + "/make_cutoffs.py"
    File cumulative_union_counts_sh = if defined(cumulative_union_counts_script)
        then select_first([cumulative_union_counts_script])
        else scripts_base + "/cumulative_union_counts.sh"
    File concat_cumulative_counts_py = if defined(concat_cumulative_counts_script)
        then select_first([concat_cumulative_counts_script])
        else scripts_base + "/concat_cumulative_counts.py"

    # 1. Sample order for the cumulative union. If ranked_samples is given, it is used
    #    as is. Otherwise samples are ranked once, from total PASS non-ref counts on a
    #    single small contig (chr22 if provided, else the shortest autosome provided);
    #    the same ranking is used for every contig in step 2a.
    if (!defined(ranked_samples)) {
        Array[String] vcf_paths = vcfs
        call PickRankingVcf {
            input:
                vcf_paths   = vcf_paths,
                docker      = python_docker,
                preemptible = preemptible
        }

        call PerSampleTypeCounts as RankingContigCounts {
            input:
                vcf                    = vcfs[PickRankingVcf.index],
                sample_population_list = sample_population_list,
                split_by_allele_type   = false,
                allele_type_tag        = allele_type_tag,
                script                 = per_sample_type_counts_sh,
                docker                 = bcftools_docker,
                mem_gb                 = mem_gb,
                disk_gb                = disk_gb,
                preemptible            = preemptible
        }

        call RankSamples {
            input:
                per_sample_tsvs        = [RankingContigCounts.per_sample_counts],
                sample_population_list = sample_population_list,
                output_basename        = output_basename,
                script                 = rank_samples_py,
                docker                 = python_docker,
                mem_gb                 = mem_gb,
                disk_gb                = disk_gb,
                preemptible            = preemptible
        }
    }

    File ordered_samples = select_first([ranked_samples, RankSamples.ordered_samples])

    call MakeCutoffs {
        input:
            ordered_samples = ordered_samples,
            script          = make_cutoffs_py,
            docker          = python_docker,
            preemptible     = preemptible
    }

    # 2a. Cumulative union burden: one task per cutoff i, over all contigs.
    scatter (i in MakeCutoffs.cutoffs) {
        call CumulativeUnionCounts {
            input:
                ordered_samples = ordered_samples,
                i               = i,
                vcfs            = vcfs,
                split_by_allele_type = split_by_allele_type,
                allele_type_tag = allele_type_tag,
                script          = cumulative_union_counts_sh,
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
            script          = concat_cumulative_counts_py,
            docker          = python_docker,
            mem_gb          = mem_gb,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    # 2b. Per-sample variant counts on every contig (always run, also when
    #     ranked_samples is given). Independent of the ranking, so it runs in
    #     parallel with steps 1 and 2a.
    scatter (vcf in vcfs) {
        call PerSampleTypeCounts {
            input:
                vcf         = vcf,
                sample_population_list = sample_population_list,
                split_by_allele_type = split_by_allele_type,
                allele_type_tag = allele_type_tag,
                script      = per_sample_type_counts_sh,
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
            script          = build_sample_variant_count_table_py,
            docker          = python_docker,
            mem_gb          = mem_gb,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    output {
        Array[Int]   n_samples_per_vcf             = PerSampleTypeCounts.n_samples
        Array[File]  per_sample_counts_per_vcf     = PerSampleTypeCounts.per_sample_counts
        File         sample_variant_count_table    = BuildSampleVariantCountTable.count_table
        String?      ranking_contig                = PickRankingVcf.contig
        File?        sample_rank_table             = RankSamples.sample_rank
        File         cumulative_sample_order       = ordered_samples
        Array[Int]   cumulative_cutoffs            = MakeCutoffs.cutoffs
        File         ranked_sample_burden_table    = ConcatCumulativeCounts.out_table
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Per-sample PASS-only non-ref variant burden across chromosome-shard VCFs: a ranked, cumulative-union-site (not summed) population-list-restricted burden table at checkpoint sample sizes (one task per cutoff), and a flat unranked genome-wide per-sample count table, all restricted to the samples in sample_population_list."
    }
}

# Picks the VCF used for ranking from the chrN token in each filename (no VCF is
# localized): chr22 if present, else the shortest GRCh38 autosome among those found.
task PickRankingVcf {
    input {
        Array[String] vcf_paths
        String        docker
        Int           preemptible
    }

    command <<<
        set -euo pipefail
        python3 <<'EOF'
import os
import re
import sys

LENGTHS = {1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
           6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
           11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
           16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
           21: 46709983, 22: 50818468}
paths = [l.rstrip("\n") for l in open("~{write_lines(vcf_paths)}") if l.strip()]
autosome = {}
for i, p in enumerate(paths):
    found = {int(n) for n in re.findall(r"chr([1-9]|1[0-9]|2[0-2])(?![0-9])", os.path.basename(p))}
    if len(found) == 1:
        autosome[i] = found.pop()
if not autosome:
    sys.exit("No VCF filename names a single autosome (chr1..chr22); cannot pick a ranking contig. "
             "Provide ranked_samples, or name per-contig VCFs with their chrN.")
chr22 = [i for i, c in autosome.items() if c == 22]
index = chr22[0] if chr22 else min(autosome, key=lambda i: LENGTHS[autosome[i]])
print(f"Ranking contig: chr{autosome[index]} ({paths[index]})")
open("index.txt", "w").write(f"{index}\n")
open("contig.txt", "w").write(f"chr{autosome[index]}\n")
EOF
    >>>

    output {
        Int    index  = read_int("index.txt")
        String contig = read_string("contig.txt")
    }

    runtime {
        docker:      docker
        memory:      "1 GB"
        cpu:         1
        disks:       "local-disk 10 HDD"
        preemptible: preemptible
    }
}

task PerSampleTypeCounts {
    input {
        File    vcf
        File    sample_population_list
        Boolean split_by_allele_type
        String  allele_type_tag
        File    script
        String  docker
        Int     mem_gb
        Int     disk_gb
        Int     preemptible
    }

    String out_prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        BY_ALLELE_TYPE=~{split_by_allele_type} ALLELE_TYPE_TAG=~{allele_type_tag} bash ~{script} ~{vcf} ~{out_prefix} ~{sample_population_list}
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
        Boolean     split_by_allele_type
        String      allele_type_tag
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
        N_PARALLEL=~{cpu} BY_ALLELE_TYPE=~{split_by_allele_type} ALLELE_TYPE_TAG=~{allele_type_tag} bash ~{script} ~{ordered_samples} ~{i} ~{out_prefix} ~{sep=" " vcfs}
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
