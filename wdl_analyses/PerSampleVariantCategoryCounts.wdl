version 1.0

## PerSampleVariantCategoryCounts
##
## Given a set of per-chromosome, jointly-genotyped VCFs from the same cohort
## (e.g. chr1.chr1.annotated.vcf.gz .. chr22.chr22.annotated.vcf.gz, as
## produced by this cohort's AnnotateVcf pipeline):
##   1. extract the sample list once, from the first VCF (bcftools query -l;
##      all contig VCFs are assumed to share the same sample set, matching
##      PerSampleVariantBurden.wdl's convention -- not enforced here, just
##      surfaced for the caller to sanity-check).
##   2-2c. for each contig VCF, restricted to FILTER=PASS: for each sample,
##       count pLoF SNVs, pLoF indels (DEL/INS, <50bp), pLoF SVs (DEL/INS,
##       >=50bp), missense, synonymous, intronic and intergenic non-ref
##       genotypes (see per_sample_category_counts.py for the exact
##       INFO/vep consequence-term matching used for each category -- these
##       are independent substring/term checks, not mutually exclusive, e.g.
##       a variant can be both pLoF on one overlapping transcript and
##       missense on another). Also collects, per sample, the set of gene
##       symbols disrupted by pLoF SNVs/indels/SVs (genes are only credited
##       from transcript annotations that themselves carry a HIGH-impact
##       consequence term).
##   3. sum each sample's counts and union each sample's pLoF gene sets
##      across all contig shards into one genome-wide table.
##
## pLoF is defined by VEP's documented HIGH-impact consequence terms
## (transcript_ablation, splice_acceptor_variant, splice_donor_variant,
## stop_gained, frameshift_variant, stop_lost, start_lost,
## transcript_amplification, feature_elongation, feature_truncation), matched
## against INFO/vep (format: Allele|Consequence|IMPACT|SYMBOL|Gene|...,
## comma-separated across overlapping transcripts), since IMPACT is not
## independently exposed as its own INFO field in this dataset.

workflow PerSampleVariantCategoryCounts {

    input {
        Array[File] vcfs
        String      output_basename
        File        per_sample_category_counts_script
        File        concat_sample_category_counts_script
        String      bcftools_docker = "quay.io/biocontainers/bcftools:1.20--h8b25389_0"
        String      python_docker   = "python:3.11-slim"
        Int         mem_gb      = 8
        Int         disk_gb     = 50
        Int         preemptible = 1
    }

    call ExtractSampleIds {
        input:
            vcf         = vcfs[0],
            docker      = bcftools_docker,
            mem_gb      = mem_gb,
            disk_gb     = disk_gb,
            preemptible = preemptible
    }

    scatter (vcf in vcfs) {
        call PerSampleCategoryCounts {
            input:
                vcf         = vcf,
                script      = per_sample_category_counts_script,
                docker      = python_docker,
                mem_gb      = mem_gb,
                disk_gb     = disk_gb,
                preemptible = preemptible
        }
    }

    call ConcatAcrossContigs {
        input:
            per_contig_tsvs = PerSampleCategoryCounts.category_counts,
            output_basename = output_basename,
            script          = concat_sample_category_counts_script,
            docker          = python_docker,
            mem_gb          = mem_gb,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    output {
        File        sample_ids                    = ExtractSampleIds.samples
        Array[File] per_contig_category_counts     = PerSampleCategoryCounts.category_counts
        File        sample_variant_category_table  = ConcatAcrossContigs.out_table
    }

    meta {
        author: "gnomAD LR analysis"
        description: "Per-sample pLoF (SNV/indel/SV)/missense/synonymous/intronic/intergenic non-ref variant counts across chromosome-shard VCFs, plus per-sample pLoF gene attribution, genome-wide."
    }
}

task ExtractSampleIds {
    input {
        File   vcf
        String docker
        Int    mem_gb
        Int    disk_gb
        Int    preemptible
    }

    command <<<
        set -euo pipefail
        bcftools query -l ~{vcf} > samples.txt
    >>>

    output {
        File samples = "samples.txt"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         1
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task PerSampleCategoryCounts {
    input {
        File   vcf
        File   script
        String docker
        Int    mem_gb
        Int    disk_gb
        Int    preemptible
    }

    String out_prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        export DEBIAN_FRONTEND=noninteractive
        apt-get update -qq && apt-get install -y -qq bcftools > /dev/null
        python3 ~{script} ~{vcf} ~{out_prefix}
    >>>

    output {
        File category_counts = out_prefix + ".category_counts.tsv"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         2
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task ConcatAcrossContigs {
    input {
        Array[File] per_contig_tsvs
        String      output_basename
        File        script
        String      docker
        Int         mem_gb
        Int         disk_gb
        Int         preemptible
    }

    String out_name = output_basename + ".sample_variant_category_table.tsv"

    command <<<
        set -euo pipefail
        python3 ~{script} ~{out_name} ~{sep=" " per_contig_tsvs}
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
