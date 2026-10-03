version 1.0

## nonref_plof_summary.wdl
##
## Reads a list of VCFs (one per contig, all containing the same sample set),
## extracts per-sample non-reference calls, parses VEP (CSQ/ANN) annotations
## and the SV-style INFO/PREDICTED_LOF field, and summarizes pLoF / missense
## burden per sample.
##
## ASSUMPTIONS (adjust workflow inputs if these don't hold for your data):
##   1. All contig VCFs are bgzipped (.vcf.gz) and contain the identical set
##      of samples (in any order). Each VCF's tabix/CSI index must be
##      supplied explicitly via `contig_vcf_indices`, in the same order as
##      `contig_vcfs` (this workflow does not (re)generate indices itself).
##   2. Variants are VEP-annotated, with the annotation stored in a single
##      INFO field (default "CSQ"; set `vep_info_field = "ANN"` if needed),
##      and a gene-symbol subfield (default "SYMBOL"; set `vep_gene_field`
##      to "Gene" if you want Ensembl gene IDs instead).
##   3. If `canonical_only = true` (default), CSQ/ANN must contain a
##      "CANONICAL" subfield (i.e. VEP was run with --canonical). Set to
##      false if that subfield is absent, or all transcripts will be
##      filtered out.
##   4. INFO/PREDICTED_LOF (default field name, override via
##      `predicted_lof_info_field`) is a comma- or "&"-separated list of
##      gene symbols at sites where an SV is predicted to cause LoF (the
##      convention used by GATK-SV / gnomAD-SV SVAnnotate output). No
##      per-allele parsing is attempted for this field, since it is already
##      site-level.
##   5. bcftools (>=1.12, with plugins, for `+split-vep`) is available in
##      the docker image used for bcftools tasks.
##
## Variant identity is CHROM:POS:REF:ALT throughout, which is what's used
## to de-duplicate and union variant lists across annotation sources.

workflow NonRefVariantPLoFSummary {
  input {
    # One VCF per contig; same sample set assumed across all of them.
    Array[File] contig_vcfs

    # Required: index file (.tbi or .csi) for each VCF in `contig_vcfs`,
    # in the same order.
    Array[File] contig_vcf_indices

    # Optional: sample IDs to process. If omitted, the sample list is
    # extracted from the first VCF in `contig_vcfs`.
    Array[String]? samples

    # VEP annotation settings
    String vep_info_field = "CSQ"
    String vep_gene_field = "SYMBOL"
    Boolean canonical_only = true
    Array[String] lof_consequences = [
      "transcript_ablation",
      "splice_acceptor_variant",
      "splice_donor_variant",
      "stop_gained",
      "frameshift_variant",
      "stop_lost",
      "start_lost"
    ]
    String missense_consequence = "missense_variant"

    # SV-style INFO field with predicted LoF gene list
    String predicted_lof_info_field = "PREDICTED_LOF"

    # Docker images
    String bcftools_docker = "staphb/bcftools:1.19"
  }

  # ---------------------------------------------------------------------
  # 1. Extract sample list from the first VCF, unless provided.
  # ---------------------------------------------------------------------
  if (!defined(samples)) {
    call GetSampleList {
      input:
        vcf       = contig_vcfs[0],
        vcf_index = contig_vcf_indices[0],
        docker    = bcftools_docker
    }
  }
  Array[String] sample_ids = select_first([samples, GetSampleList.samples])

  # ---------------------------------------------------------------------
  # 2-5. Per-sample: pull non-ref calls across all contigs, then extract
  #      VEP pLoF / missense variants+genes and PREDICTED_LOF variants+genes.
  # ---------------------------------------------------------------------
  scatter (sample_id in sample_ids) {

    call ExtractSampleNonRef {
      input:
        sample_id   = sample_id,
        vcfs        = contig_vcfs,
        vcf_indices = contig_vcf_indices,
        docker      = bcftools_docker
    }

    call ExtractVepConsequenceVariants as ExtractPlofVep {
      input:
        sample_id        = sample_id,
        vcf              = ExtractSampleNonRef.sample_nonref_vcf,
        vcf_index        = ExtractSampleNonRef.sample_nonref_vcf_index,
        vep_info_field   = vep_info_field,
        vep_gene_field   = vep_gene_field,
        consequence_terms = lof_consequences,
        canonical_only   = canonical_only,
        label            = "plof",
        docker           = bcftools_docker
    }

    call ExtractVepConsequenceVariants as ExtractMissenseVep {
      input:
        sample_id        = sample_id,
        vcf              = ExtractSampleNonRef.sample_nonref_vcf,
        vcf_index        = ExtractSampleNonRef.sample_nonref_vcf_index,
        vep_info_field   = vep_info_field,
        vep_gene_field   = vep_gene_field,
        consequence_terms = [missense_consequence],
        canonical_only   = canonical_only,
        label            = "missense",
        docker           = bcftools_docker
    }

    call ExtractPredictedLofField {
      input:
        sample_id  = sample_id,
        vcf        = ExtractSampleNonRef.sample_nonref_vcf,
        vcf_index  = ExtractSampleNonRef.sample_nonref_vcf_index,
        info_field = predicted_lof_info_field,
        docker     = bcftools_docker
    }

    # -------------------------------------------------------------------
    # 6. Per-sample summary row (includes the union columns 6 & 7).
    # -------------------------------------------------------------------
    call CombineSampleSummary {
      input:
        sample_id                   = sample_id,
        n_plof_vep_variants         = ExtractPlofVep.n_variants,
        n_plof_vep_genes            = ExtractPlofVep.n_genes,
        n_predicted_lof_variants    = ExtractPredictedLofField.n_variants,
        n_predicted_lof_genes       = ExtractPredictedLofField.n_genes,
        n_missense_vep_variants     = ExtractMissenseVep.n_variants,
        n_missense_vep_genes        = ExtractMissenseVep.n_genes,
        plof_vep_variants_list      = ExtractPlofVep.variants_list,
        predicted_lof_variants_list = ExtractPredictedLofField.variants_list,
        plof_vep_genes_list         = ExtractPlofVep.genes_list,
        predicted_lof_genes_list    = ExtractPredictedLofField.genes_list,
        docker                      = bcftools_docker
    }
  }

  call MergeSummaryTables {
    input:
      summary_rows = CombineSampleSummary.summary_row,
      docker       = bcftools_docker
  }

  output {
    File         summary_table              = MergeSummaryTables.final_summary_table
    Array[File]  sample_nonref_vcfs          = ExtractSampleNonRef.sample_nonref_vcf
    Array[File]  sample_plof_vep_tables      = ExtractPlofVep.variant_gene_table
    Array[File]  sample_missense_vep_tables  = ExtractMissenseVep.variant_gene_table
    Array[File]  sample_predicted_lof_tables = ExtractPredictedLofField.variant_gene_table
  }
}

# ===========================================================================
# TASKS
# ===========================================================================

task GetSampleList {
  input {
    File vcf
    File vcf_index
    String docker
  }

  command <<<
    set -euo pipefail
    ln -s ~{vcf} in.vcf.gz
    ln -s ~{vcf_index} in.vcf.gz.tbi
    bcftools query -l in.vcf.gz > samples.txt
  >>>

  output {
    Array[String] samples = read_lines("samples.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 20 HDD"
  }
}

# Subsets each contig VCF down to one sample's non-ref genotypes
# (`bcftools view -s <sample> -c 1`, the standard idiom: -c/--min-ac 1 is
# evaluated *after* subsetting, so it keeps only sites where that sample
# is het or hom-alt for at least one ALT allele), then concatenates the
# per-contig results (in the same order as `contig_vcfs` was supplied) into
# a single per-sample VCF.
task ExtractSampleNonRef {
  input {
    String sample_id
    Array[File] vcfs
    Array[File] vcf_indices
    String docker
  }
  Int disk_gb = ceil(size(vcfs, "GB") * 2) + 20

  command <<<
    set -euo pipefail
    paste ~{write_lines(vcfs)} ~{write_lines(vcf_indices)} > manifest.tsv

    mkdir -p persample
    : > filelist.txt
    i=0
    while IFS=$'\t' read -r vcf_path idx_path; do
      i=$((i+1))
      ln -sf "${vcf_path}" "contig_${i}.vcf.gz"
      ln -sf "${idx_path}" "contig_${i}.vcf.gz.tbi"
      bcftools view -s "~{sample_id}" --force-samples -c 1 \
        -Oz -o "persample/nonref_${i}.vcf.gz" "contig_${i}.vcf.gz"
      tabix -p vcf "persample/nonref_${i}.vcf.gz"
      echo "persample/nonref_${i}.vcf.gz" >> filelist.txt
    done < manifest.tsv

    bcftools concat -a -f filelist.txt -Oz -o "~{sample_id}.nonref.vcf.gz"
    tabix -p vcf "~{sample_id}.nonref.vcf.gz"
  >>>

  output {
    File sample_nonref_vcf       = "~{sample_id}.nonref.vcf.gz"
    File sample_nonref_vcf_index = "~{sample_id}.nonref.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
}

# Uses `bcftools +split-vep` to explode the VEP CSQ/ANN annotation into one
# row per transcript, filters to rows whose Consequence matches any of
# `consequence_terms` (substring match, so it also catches compound
# annotations like "missense_variant&splice_region_variant"), optionally
# restricted to the canonical transcript, and reports the deduplicated
# variant list and gene list.
task ExtractVepConsequenceVariants {
  input {
    String sample_id
    File vcf
    File vcf_index
    String vep_info_field
    String vep_gene_field
    Array[String] consequence_terms
    Boolean canonical_only
    String label
    String docker
  }
  Int disk_gb = ceil(size(vcf, "GB") * 3) + 10

  command <<<
    set -euo pipefail
    ln -s ~{vcf} in.vcf.gz
    ln -s ~{vcf_index} in.vcf.gz.tbi

    TERMS=(~{sep=" " consequence_terms})
    EXPR=""
    for t in "${TERMS[@]}"; do
      if [ -z "${EXPR}" ]; then
        EXPR="Consequence~\"${t}\""
      else
        EXPR="${EXPR} || Consequence~\"${t}\""
      fi
    done

    CANON_FILTER=""
    if [ "~{canonical_only}" = "true" ]; then
      CANON_FILTER=" && CANONICAL=\"YES\""
    fi

    FULL_EXPR="(${EXPR})${CANON_FILTER}"

    bcftools +split-vep -a ~{vep_info_field} \
      -f "%CHROM:%POS:%REF:%ALT\t%~{vep_gene_field}\n" \
      -d \
      -i "${FULL_EXPR}" \
      in.vcf.gz > raw_hits.tsv || true

    sort -u raw_hits.tsv > "~{sample_id}.~{label}.variant_gene.tsv"
    cut -f1 "~{sample_id}.~{label}.variant_gene.tsv" | sort -u \
      > "~{sample_id}.~{label}.variants.txt"
    cut -f2 "~{sample_id}.~{label}.variant_gene.tsv" | sort -u \
      | sed '/^\.$/d;/^$/d' > "~{sample_id}.~{label}.genes.txt"

    wc -l < "~{sample_id}.~{label}.variants.txt" | tr -d ' ' > n_variants.txt
    wc -l < "~{sample_id}.~{label}.genes.txt" | tr -d ' ' > n_genes.txt
  >>>

  output {
    File variant_gene_table = "~{sample_id}.~{label}.variant_gene.tsv"
    File variants_list      = "~{sample_id}.~{label}.variants.txt"
    File genes_list         = "~{sample_id}.~{label}.genes.txt"
    Int  n_variants         = read_int("n_variants.txt")
    Int  n_genes            = read_int("n_genes.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
}

# Parses INFO/PREDICTED_LOF (or whichever field name is supplied), which is
# expected to already be a site-level, comma/"&"-separated list of gene
# symbols (no per-allele CSQ-style parsing needed/possible for this field).
task ExtractPredictedLofField {
  input {
    String sample_id
    File vcf
    File vcf_index
    String info_field
    String docker
  }
  Int disk_gb = ceil(size(vcf, "GB") * 3) + 10

  command <<<
    set -euo pipefail
    ln -s ~{vcf} in.vcf.gz
    ln -s ~{vcf_index} in.vcf.gz.tbi

    bcftools query -f "%CHROM:%POS:%REF:%ALT\t%INFO/~{info_field}\n" in.vcf.gz \
      | awk -F'\t' '$2!="." && $2!=""' > raw.tsv

    awk -F'\t' '{
      n = split($2, genes, /[,&]/);
      for (i = 1; i <= n; i++) {
        if (genes[i] != "" && genes[i] != ".") print $1"\t"genes[i];
      }
    }' raw.tsv | sort -u > "~{sample_id}.predicted_lof.variant_gene.tsv"

    cut -f1 "~{sample_id}.predicted_lof.variant_gene.tsv" | sort -u \
      > "~{sample_id}.predicted_lof.variants.txt"
    cut -f2 "~{sample_id}.predicted_lof.variant_gene.tsv" | sort -u \
      > "~{sample_id}.predicted_lof.genes.txt"

    wc -l < "~{sample_id}.predicted_lof.variants.txt" | tr -d ' ' > n_variants.txt
    wc -l < "~{sample_id}.predicted_lof.genes.txt" | tr -d ' ' > n_genes.txt
  >>>

  output {
    File variant_gene_table = "~{sample_id}.predicted_lof.variant_gene.tsv"
    File variants_list      = "~{sample_id}.predicted_lof.variants.txt"
    File genes_list         = "~{sample_id}.predicted_lof.genes.txt"
    Int  n_variants         = read_int("n_variants.txt")
    Int  n_genes            = read_int("n_genes.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
}

# Builds one summary row per sample: columns 1-5 and 8-9 are passed straight
# through as counts; columns 6-7 are the union of the VEP-pLoF and
# PREDICTED_LOF variant/gene lists, computed here.
task CombineSampleSummary {
  input {
    String sample_id
    Int n_plof_vep_variants
    Int n_plof_vep_genes
    Int n_predicted_lof_variants
    Int n_predicted_lof_genes
    Int n_missense_vep_variants
    Int n_missense_vep_genes
    File plof_vep_variants_list
    File predicted_lof_variants_list
    File plof_vep_genes_list
    File predicted_lof_genes_list
    String docker
  }

  command <<<
    set -euo pipefail
    UNION_VAR=$(cat ~{plof_vep_variants_list} ~{predicted_lof_variants_list} | sort -u | wc -l | tr -d ' ')
    UNION_GENE=$(cat ~{plof_vep_genes_list} ~{predicted_lof_genes_list} | sort -u | wc -l | tr -d ' ')

    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" \
      "~{sample_id}" \
      "~{n_plof_vep_variants}" \
      "~{n_plof_vep_genes}" \
      "~{n_predicted_lof_variants}" \
      "~{n_predicted_lof_genes}" \
      "${UNION_VAR}" \
      "${UNION_GENE}" \
      "~{n_missense_vep_variants}" \
      "~{n_missense_vep_genes}" \
      > "~{sample_id}.summary.tsv"
  >>>

  output {
    File summary_row = "~{sample_id}.summary.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
}

task MergeSummaryTables {
  input {
    Array[File] summary_rows
    String docker
  }

  command <<<
    set -euo pipefail
    printf "sample_id\tn_plof_variants_vep\tn_plof_genes_vep\tn_lof_variants_predicted_lof\tn_genes_predicted_lof\tn_variants_union_plof_predictedlof\tn_genes_union_plof_predictedlof\tn_missense_variants_vep\tn_missense_genes_vep\n" > final_summary.tsv
    cat ~{sep=" " summary_rows} >> final_summary.tsv
  >>>

  output {
    File final_summary_table = "final_summary.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
}
