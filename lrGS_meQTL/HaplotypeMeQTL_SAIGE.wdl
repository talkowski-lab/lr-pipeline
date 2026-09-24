version 1.0

# Cis-meQTL scan: per-haplotype genotype vs. per-haplotype methylation level.
#
# Identical to GenotypeMeQTL.wdl except each phased contig VCF is first
# rewritten (Tasks.SplitPhasedVcfToHaplotypes) into a pseudo-haploid VCF
# with 2x the samples ("S_hap1", "S_hap2" per original sample "S", each a
# homozygous pseudo-diploid encoding of that one haplotype's allele). LD
# pruning, the sparse GRM, and the per-locus SAIGE step1/step2 calls then
# run on that expanded sample set unchanged, matched against a per-haplotype
# methylation file (e.g. hprc_methylated.chr22.haplotype.bed.gz) whose
# sample columns already use the same "*_hap1" / "*_hap2" naming.
#
# `vcfs`, `methylation_files`, and `contigs` are parallel arrays: index i
# in all three must refer to the same contig. `vcfs` must be phased.

import "Tasks.meQTL.wdl" as Tasks

workflow HaplotypeMeQTL_SAIGE {
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
        File methylation_file = methylation_files[i]
        String contig = contigs[i]

        call Tasks.SplitPhasedVcfToHaplotypes {
            input:
                vcf = vcfs[i],
                contig = contig,
                docker = bcftools_docker
        }

        call Tasks.LdPruneAndExtract {
            input:
                vcf = SplitPhasedVcfToHaplotypes.haplotype_vcf,
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
                    vcf = SplitPhasedVcfToHaplotypes.haplotype_vcf,
                    vcf_csi = SplitPhasedVcfToHaplotypes.haplotype_vcf_csi,
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
                out_prefix = contig + ".cis_meQTL.haplotype",
                docker = bcftools_docker
        }
    }

    call Tasks.ConcatenateTsvs as ConcatenateAllContigs {
        input:
            tsvs = ConcatenateContigChunks.merged,
            out_prefix = "cis_meQTL.haplotype.all_contigs",
            docker = bcftools_docker
    }

    output {
        Array[File] per_contig_assoc = ConcatenateContigChunks.merged
        File combined_assoc = ConcatenateAllContigs.merged
        Array[File] per_contig_haplotype_vcf = SplitPhasedVcfToHaplotypes.haplotype_vcf
        Array[File] per_contig_pruned_bed = LdPruneAndExtract.bed
        Array[File] per_contig_sparse_grm = CreateSparseGRM.sparse_grm
        Array[Array[File]] per_chunk_skipped_logs = RunCisMeQTLChunk.skipped_sites_log
    }
}
