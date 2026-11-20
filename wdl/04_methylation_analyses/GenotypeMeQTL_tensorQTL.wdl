version 1.0

# Cis-meQTL scan: diploid genotype dosage vs. per-sample methylation level,
# via tensorQTL (permutation-based cis mapping), wrapping the core task from
# AoU-Multiomics-Analysis/tensorQTL_cis_permutations. That upstream workflow
# takes pre-made plink2 files + a phenotype bed + covariates directly; this
# one instead takes a VCF and a wide-format methylation table per contig and
# builds those inputs itself (bcftools norm, plink2 --make-pgen, call-rate
# filtering + mean imputation into a phenotype bed, and a covariates file).
#
# Per contig (scattered): split multiallelic sites to biallelic, convert to
# plink2 pgen/pvar/psam, filter+reshape the methylation table into a
# phenotype bed, build a covariates file, then run one tensorQTL cis
# permutation pass over the whole contig (no per-site chunking - tensorQTL
# is vectorized/GPU-based, unlike SAIGE's per-site null model). Per-contig
# results are gathered into one genome-wide table.
#
# `vcfs`, `methylation_files`, and `contigs` are parallel arrays: index i
# in all three must refer to the same contig.

import "MeQTLTensorQTLTasks.wdl"

workflow GenotypeMeQTL_tensorQTL {
    input {
        Array[File] vcfs
        Array[File] methylation_files
        Array[String] contigs
        String prefix

        Float min_call_rate = 0.9
        Int cis_window = 1000000
        String vcf_half_call = "missing"
        File? covariates_file
        File? phenotype_groups
        Float? fdr
        Float? qvalue_lambda
        Float? pval_threshold
        Int? seed
        String? flags

        String plink2_docker = "quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0"
        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
        String python_docker = "wzhou88/saige:1.3.6"
        String tensorqtl_docker = "gcr.io/broad-cga-francois-gtex/tensorqtl:latest"

        RuntimeAttr? runtime_attr_normalize_vcf
        RuntimeAttr? runtime_attr_convert_pgen
        RuntimeAttr? runtime_attr_build_phenotype_bed
        RuntimeAttr? runtime_attr_bgzip_tabix
        RuntimeAttr? runtime_attr_build_covariates
        RuntimeAttr? runtime_attr_run_tensorqtl
        RuntimeAttr? runtime_attr_concat_contigs
    }

    Int n_contigs = length(contigs)

    scatter (i in range(n_contigs)) {
        File methylation_file = methylation_files[i]
        String contig = contigs[i]
        String contig_prefix = "~{prefix}.~{contig}"

        call MeQTLTensorQTLTasks.NormalizeVcf {
            input:
                vcf = vcfs[i],
                prefix = contig_prefix,
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_normalize_vcf
        }

        call MeQTLTensorQTLTasks.ConvertVcfToPgen {
            input:
                vcf = NormalizeVcf.normalized_vcf,
                vcf_half_call = vcf_half_call,
                prefix = contig_prefix,
                docker = plink2_docker,
                runtime_attr_override = runtime_attr_convert_pgen
        }

        call MeQTLTensorQTLTasks.BuildPhenotypeBed {
            input:
                methylation_bed = methylation_file,
                call_rate_threshold = min_call_rate,
                prefix = contig_prefix,
                docker = python_docker,
                runtime_attr_override = runtime_attr_build_phenotype_bed
        }

        call MeQTLTensorQTLTasks.BgzipTabixBed {
            input:
                bed = BuildPhenotypeBed.phenotype_bed_plain,
                prefix = contig_prefix,
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_bgzip_tabix
        }

        call MeQTLTensorQTLTasks.BuildCovariates {
            input:
                phenotype_bed_plain = BuildPhenotypeBed.phenotype_bed_plain,
                covariates_file = covariates_file,
                prefix = contig_prefix,
                docker = python_docker,
                runtime_attr_override = runtime_attr_build_covariates
        }

        call MeQTLTensorQTLTasks.TensorQTLCisPermutations {
            input:
                plink_pgen = ConvertVcfToPgen.pgen,
                plink_pvar = ConvertVcfToPgen.pvar,
                plink_psam = ConvertVcfToPgen.psam,
                phenotype_bed = BgzipTabixBed.bed_gz,
                phenotype_bed_index = BgzipTabixBed.bed_gz_index,
                covariates = BuildCovariates.covariates,
                cis_window = cis_window,
                phenotype_groups = phenotype_groups,
                fdr = fdr,
                qvalue_lambda = qvalue_lambda,
                pval_threshold = pval_threshold,
                seed = seed,
                flags = flags,
                prefix = contig_prefix,
                docker = tensorqtl_docker,
                runtime_attr_override = runtime_attr_run_tensorqtl
        }
    }

    call MeQTLTensorQTLTasks.ConcatenateGzippedTsvs {
        input:
            gz_tsvs = TensorQTLCisPermutations.cis_qtl,
            prefix = "~{prefix}.cis_meQTL.tensorqtl.genotype.all_contigs",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_contigs
    }

    output {
        Array[File] per_contig_cis_qtl = TensorQTLCisPermutations.cis_qtl
        Array[File] per_contig_log = TensorQTLCisPermutations.log
        File combined_cis_qtl = ConcatenateGzippedTsvs.merged
        Array[File] per_contig_phenotype_bed = BgzipTabixBed.bed_gz
        Array[File] per_contig_covariates = BuildCovariates.covariates
        Array[Int] per_contig_n_sites = BuildPhenotypeBed.n_sites
    }
}
