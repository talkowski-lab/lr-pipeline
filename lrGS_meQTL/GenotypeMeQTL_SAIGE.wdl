version 1.0

# Cis-meQTL scan: diploid genotype dosage vs. per-sample methylation level.
#
# Per contig (scattered): LD-prune + extract a plink bfile from the contig
# VCF, build a sparse GRM from that bfile, filter the contig's methylation
# sites to those with >= min_call_rate non-missing samples, then chunk the
# qualifying sites and - per chunk - run SAIGE step1_fitNULLGLMM.R (one
# null model per site) + step2_SPAtests.R (restricted to a +/- cis_window
# region around the site) via RunCisMeQTLChunk. Chunk results are gathered
# per contig, then across contigs into one final association table.
#
# `vcfs`, `methylation_files`, and `contigs` are parallel arrays: index i
# in all three must refer to the same contig.

import "Tasks.meQTL.wdl" as Tasks

workflow GenotypeMeQTL_SAIGE {
    input {
        Array[File] vcfs
        Array[File] methylation_files
        Array[String] contigs

        Float min_call_rate = 0.9
        Int cis_window = 2000000
        Int sites_per_shard = 200

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

        String plink_docker = "quay.io/biocontainers/plink:1.90b6.21--h031d066_5"
        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
        String saige_docker = "wzhou88/saige:1.3.6"
    }

    Int n_contigs = length(contigs)

    scatter (i in range(n_contigs)) {
        File vcf = vcfs[i]
        File methylation_file = methylation_files[i]
        String contig = contigs[i]

        call Tasks.IndexVcf {
            input:
                vcf = vcf,
                docker = bcftools_docker
        }

        call Tasks.LdPruneAndExtract {
            input:
                vcf = vcf,
                contig = contig,
                window_size_kb = ld_prune_window_kb,
                step_size = ld_prune_step,
                r2_threshold = ld_prune_r2,
                vcf_half_call = vcf_half_call,
                docker = plink_docker
        }

        call Tasks.CreateSparseGRM {
            input:
                bed = LdPruneAndExtract.bed,
                bim = LdPruneAndExtract.bim,
                fam = LdPruneAndExtract.fam,
                num_random_markers = num_random_markers_for_grm,
                relatedness_cutoff = relatedness_cutoff,
                min_maf_for_grm = min_maf_for_grm,
                max_missing_rate_for_grm = max_missing_rate_for_grm,
                docker = saige_docker
        }

        call Tasks.FilterMethylationSites {
            input:
                methylation_bed = methylation_file,
                call_rate_threshold = min_call_rate,
                docker = saige_docker
        }

        Int n_chunks = (FilterMethylationSites.n_sites + sites_per_shard - 1) / sites_per_shard

        scatter (chunk_idx in range(n_chunks)) {
            Int row_start = chunk_idx * sites_per_shard
            Int row_end = if (row_start + sites_per_shard) < FilterMethylationSites.n_sites
                          then row_start + sites_per_shard
                          else FilterMethylationSites.n_sites

            call Tasks.RunCisMeQTLChunk {
                input:
                    filtered_sites = FilterMethylationSites.filtered_sites,
                    row_start = row_start,
                    row_end = row_end,
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
                    docker = saige_docker
            }
        }

        call Tasks.ConcatenateTsvs as ConcatenateContigChunks {
            input:
                tsvs = RunCisMeQTLChunk.chunk_assoc,
                out_prefix = contig + ".cis_meQTL",
                docker = bcftools_docker
        }
    }

    call Tasks.ConcatenateTsvs as ConcatenateAllContigs {
        input:
            tsvs = ConcatenateContigChunks.merged,
            out_prefix = "cis_meQTL.genotype.all_contigs",
            docker = bcftools_docker
    }

    output {
        Array[File] per_contig_assoc = ConcatenateContigChunks.merged
        File combined_assoc = ConcatenateAllContigs.merged
        Array[File] per_contig_pruned_bed = LdPruneAndExtract.bed
        Array[File] per_contig_sparse_grm = CreateSparseGRM.sparse_grm
        Array[Array[File]] per_chunk_skipped_logs = RunCisMeQTLChunk.skipped_sites_log
    }
}
