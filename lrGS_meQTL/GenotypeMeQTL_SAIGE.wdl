version 1.0

# Cis-meQTL scan: diploid genotype dosage vs. per-sample methylation level,
# via plink (LD-pruning, sparse-GRM marker set) + SAIGE (sparse-GRM null
# model + SPA test).
#
# Per contig (scattered): LD-prune + extract a plink bfile from the contig
# VCF, build a sparse GRM from that bfile, filter the contig's methylation
# sites to those with >= min_call_rate non-missing samples, then run SAIGE
# step1_fitNULLGLMM.R (one null model per site) + step2_SPAtests.R
# (restricted to a +/- cis_window region around the site) for every
# qualifying site via RunCisMeQTLContig - one task per contig, sites
# parallelized across CPU cores within it (see MeQTLTasks.wdl for why: a
# second, nested WDL scatter over site chunks hit a Cromwell bug resolving
# optional inputs two scatter levels deep). Per-contig results are gathered
# across contigs into one final association table.
#
# `vcfs`, `methylation_files`, and `contigs` are parallel arrays: index i
# in all three must refer to the same contig.

import "MeQTLTasks.wdl"

workflow GenotypeMeQTL_SAIGE {
    input {
        Array[File] vcfs
        Array[File] methylation_files
        Array[String] contigs
        String prefix

        Float min_call_rate = 0.9
        Int cis_window = 2000000
        Int ld_prune_window_kb = 50
        Int ld_prune_step = 5
        Float ld_prune_r2 = 0.2
        String vcf_half_call = "missing"
        Int num_random_markers_for_grm = 2000
        Float relatedness_cutoff = 0.125
        Float min_maf_for_grm = 0.01
        Float max_missing_rate_for_grm = 0.15
        File? covariates_file
        String covar_col_list = ""
        String qcovar_col_list = ""
        Int min_samples_per_site = 20
        Boolean inv_normalize = true
        Int n_parallel_workers = 8

        String plink_docker = "quay.io/biocontainers/plink:1.90b6.21--h031d066_5"
        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
        String saige_docker = "wzhou88/saige:1.3.6"

        RuntimeAttr? runtime_attr_index_vcf
        RuntimeAttr? runtime_attr_ld_prune
        RuntimeAttr? runtime_attr_create_grm
        RuntimeAttr? runtime_attr_filter_sites
        RuntimeAttr? runtime_attr_run_contig
        RuntimeAttr? runtime_attr_concat_contigs
    }

    Int n_contigs = length(contigs)

    scatter (i in range(n_contigs)) {
        File vcf = vcfs[i]
        File methylation_file = methylation_files[i]
        String contig = contigs[i]
        String contig_prefix = "~{prefix}.~{contig}"

        call MeQTLTasks.IndexVcf {
            input:
                vcf = vcf,
                prefix = contig_prefix,
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_index_vcf
        }

        call MeQTLTasks.LdPruneAndExtract {
            input:
                vcf = vcf,
                prefix = contig_prefix,
                window_size_kb = ld_prune_window_kb,
                step_size = ld_prune_step,
                r2_threshold = ld_prune_r2,
                vcf_half_call = vcf_half_call,
                docker = plink_docker,
                runtime_attr_override = runtime_attr_ld_prune
        }

        call MeQTLTasks.CreateSparseGRM {
            input:
                bed = LdPruneAndExtract.bed,
                bim = LdPruneAndExtract.bim,
                fam = LdPruneAndExtract.fam,
                prefix = contig_prefix,
                num_random_markers = num_random_markers_for_grm,
                relatedness_cutoff = relatedness_cutoff,
                min_maf_for_grm = min_maf_for_grm,
                max_missing_rate_for_grm = max_missing_rate_for_grm,
                n_threads = 4,
                docker = saige_docker,
                runtime_attr_override = runtime_attr_create_grm
        }

        call MeQTLTasks.FilterMethylationSites {
            input:
                methylation_bed = methylation_file,
                call_rate_threshold = min_call_rate,
                prefix = contig_prefix,
                docker = saige_docker,
                runtime_attr_override = runtime_attr_filter_sites
        }

        call MeQTLTasks.RunCisMeQTLContig {
            input:
                filtered_sites = FilterMethylationSites.filtered_sites,
                contig = contig,
                vcf = vcf,
                vcf_csi = IndexVcf.vcf_csi,
                pruned_bed = LdPruneAndExtract.bed,
                pruned_bim = LdPruneAndExtract.bim,
                pruned_fam = LdPruneAndExtract.fam,
                sparse_grm = CreateSparseGRM.sparse_grm,
                sparse_grm_samples = CreateSparseGRM.sparse_grm_samples,
                relatedness_cutoff = relatedness_cutoff,
                cis_window = cis_window,
                covariates_file = covariates_file,
                covar_col_list = covar_col_list,
                qcovar_col_list = qcovar_col_list,
                min_samples_per_site = min_samples_per_site,
                vcf_field = "GT",
                inv_normalize = inv_normalize,
                n_parallel_workers = n_parallel_workers,
                prefix = "~{contig_prefix}.cis_meQTL",
                docker = saige_docker,
                runtime_attr_override = runtime_attr_run_contig
        }
    }

    call MeQTLTasks.ConcatenateTsvs as ConcatenateAllContigs {
        input:
            tsvs = RunCisMeQTLContig.contig_assoc,
            prefix = "~{prefix}.cis_meQTL.genotype.all_contigs",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_contigs
    }

    output {
        Array[File] per_contig_assoc = RunCisMeQTLContig.contig_assoc
        File combined_assoc = ConcatenateAllContigs.merged
        Array[File] per_contig_pruned_bed = LdPruneAndExtract.bed
        Array[File] per_contig_sparse_grm = CreateSparseGRM.sparse_grm
        Array[File] per_contig_skipped_logs = RunCisMeQTLContig.skipped_sites_log
    }
}
