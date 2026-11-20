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
##      Each category (plof, missense, synonymous, intronic, intergenic) runs
##      as its own independent task per contig VCF, all in parallel; each
##      task fills only its own category's columns (others are 0 / empty).
##   3. ConcatCategoryCounts: sum each sample's counts and union each
##      sample's pLoF gene sets across all (contig, category) task outputs
##      into one genome-wide table with every category's columns.
##   4. Steps 2-3 are run a second time additionally restricted to rare sites
##      (MAX(INFO/<af_field>) < max_af, default cohort INFO/AF < 0.01; a
##      multiallelic site is kept only if all its ALT alleles are rare),
##      producing sample_variant_category_table_rare.
##
## pLoF is defined by VEP's documented HIGH-impact consequence terms
## (transcript_ablation, splice_acceptor_variant, splice_donor_variant,
## stop_gained, frameshift_variant, stop_lost, start_lost,
## transcript_amplification, feature_elongation, feature_truncation), matched
## against INFO/vep (format: Allele|Consequence|IMPACT|SYMBOL|Gene|...,
## comma-separated across overlapping transcripts), since IMPACT is not
## independently exposed as its own INFO field in this dataset.
##
## Step 2's per-contig counting itself runs in parallel across
## n_chunks_per_contig genomic-position chunks (per_sample_category_counts_
## parallel.sh): each chunk is sliced out with `bcftools view -r`, which is
## tabix-index-seekable, so N concurrent chunk workers divide the contig's
## total work by N rather than each rescanning the whole file -- true data
## parallelism. (An earlier design instead ran each of the 9 count
## categories + 3 gene-list categories as its own concurrent `bcftools view
## -i <category filter>` pass over the full file; measured SLOWER in
## practice, since an arbitrary INFO-field `-i` filter can't use the index
## and still requires a full linear scan, so parallelizing across
## categories just multiplied total I/O instead of dividing it.)

workflow PerSampleVariantCategoryCounts {

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String      output_basename
        File        per_sample_category_counts_script
        File        per_sample_category_counts_parallel_script
        File        concat_sample_category_counts_script
        String      bcftools_docker = "quay.io/biocontainers/bcftools@sha256:badc3a0c7af72a83e5761ab0e881aa84204694bdead003b47552cb283958f78d"
        String      python_docker   = "python:3.11-slim"
        Array[String] categories = ["plof", "missense", "synonymous", "intronic", "intergenic"]
        String      af_field = "AF"
        String      max_af   = "0.01"
        Int         n_chunks_per_contig = 6
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

    # One independent task per (contig VCF, category): each data collection
    # runs as its own job, in parallel. Flattened into a single scatter via
    # cross() rather than nesting scatters.
    scatter (job in cross(range(length(vcfs)), categories)) {
        call PerSampleCategoryCounts {
            input:
                vcf             = vcfs[job.left],
                vcf_idx         = vcf_idxs[job.left],
                category        = job.right,
                per_sample_script = per_sample_category_counts_script,
                parallel_script   = per_sample_category_counts_parallel_script,
                concat_script     = concat_sample_category_counts_script,
                n_chunks        = n_chunks_per_contig,
                docker          = python_docker,
                mem_gb          = mem_gb,
                disk_gb         = disk_gb,
                preemptible     = preemptible
        }
    }

    # Same per-(contig, category) tasks, additionally restricted to rare
    # sites: MAX(INFO/<af_field>) < max_af.
    scatter (job in cross(range(length(vcfs)), categories)) {
        call PerSampleCategoryCounts as PerSampleCategoryCountsRare {
            input:
                vcf             = vcfs[job.left],
                vcf_idx         = vcf_idxs[job.left],
                category        = job.right,
                af_field        = af_field,
                max_af          = max_af,
                per_sample_script = per_sample_category_counts_script,
                parallel_script   = per_sample_category_counts_parallel_script,
                concat_script     = concat_sample_category_counts_script,
                n_chunks        = n_chunks_per_contig,
                docker          = python_docker,
                mem_gb          = mem_gb,
                disk_gb         = disk_gb,
                preemptible     = preemptible
        }
    }

    call ConcatCategoryCounts {
        input:
            category_count_tsvs = PerSampleCategoryCounts.category_counts,
            output_basename = output_basename,
            script          = concat_sample_category_counts_script,
            docker          = python_docker,
            mem_gb          = mem_gb,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    call ConcatCategoryCounts as ConcatCategoryCountsRare {
        input:
            category_count_tsvs = PerSampleCategoryCountsRare.category_counts,
            output_basename = output_basename + ".PASS_" + af_field + "_lt_" + max_af,
            script          = concat_sample_category_counts_script,
            docker          = python_docker,
            mem_gb          = mem_gb,
            disk_gb         = disk_gb,
            preemptible     = preemptible
    }

    output {
        File        sample_ids                    = ExtractSampleIds.samples
        Array[File] per_contig_per_category_counts = PerSampleCategoryCounts.category_counts
        Array[File] per_contig_per_category_counts_rare = PerSampleCategoryCountsRare.category_counts
        File        sample_variant_category_table_rare = ConcatCategoryCountsRare.out_table
        File        sample_variant_category_table  = ConcatCategoryCounts.out_table
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
        File   vcf_idx
        String category
        String af_field = "none"
        String max_af   = "none"
        File   per_sample_script
        File   parallel_script
        File   concat_script
        Int    n_chunks
        String docker
        Int    mem_gb
        Int    disk_gb
        Int    preemptible
    }

    String af_suffix  = if max_af == "none" then "" else "." + af_field + "_lt_" + max_af
    String out_prefix = basename(vcf, ".vcf.gz") + "." + category + af_suffix

    command <<<
        set -euo pipefail
        export DEBIAN_FRONTEND=noninteractive
        apt-get update -qq && apt-get install -y -qq bcftools tabix > /dev/null
        # Fail fast on stale staged scripts that predate per-category runs:
        # they would ignore the category arg and count every category,
        # inflating the summed counts in ConcatCategoryCounts.
        # Likewise for the AF restriction: a stale script would silently
        # ignore it and the "rare" table would equal the unrestricted one.
        if ! grep -q 'MAX_AF=' ~{parallel_script} || ! grep -q 'af_field, max_af' ~{per_sample_script}; then
            echo "ERROR: per_sample_category_counts scripts do not support per-category / AF-restricted runs; re-upload them" >&2
            exit 1
        fi
        # tabix -l / bcftools view -r need the index next to the VCF; the
        # index may have been localized to a different directory.
        ln -s ~{vcf} input.vcf.gz
        ln -s ~{vcf_idx} input.vcf.gz.tbi
        bash ~{parallel_script} input.vcf.gz ~{out_prefix} ~{n_chunks} ~{per_sample_script} ~{concat_script} ~{category} ~{af_field} ~{max_af}
    >>>

    output {
        File category_counts = out_prefix + ".category_counts.tsv"
    }

    runtime {
        docker:      docker
        memory:      mem_gb + " GB"
        cpu:         n_chunks
        disks:       "local-disk " + disk_gb + " HDD"
        preemptible: preemptible
    }
}

task ConcatCategoryCounts {
    input {
        Array[File] category_count_tsvs
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
        python3 ~{script} ~{out_name} ~{sep=" " category_count_tsvs}
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
