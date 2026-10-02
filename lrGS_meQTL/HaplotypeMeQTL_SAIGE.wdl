version 1.0

# Cis-meQTL scan: per-haplotype genotype vs. per-haplotype methylation level,
# via plink (LD-pruning, sparse-GRM marker set) + SAIGE (sparse-GRM null
# model + SPA test).
#
# Identical to GenotypeMeQTL_SAIGE.wdl except each phased contig VCF is
# first rewritten (MeQTLTasks.SplitPhasedVcfToHaplotypes) into a
# pseudo-haploid VCF with 2x the samples ("S_hap1", "S_hap2" per original
# sample "S", each a homozygous pseudo-diploid encoding of that one
# haplotype's allele). LD pruning, the sparse GRM, and the per-locus SAIGE
# step1/step2 calls then run on that expanded sample set unchanged, matched
# against a per-haplotype methylation file (e.g.
# hprc_methylated.chr22.haplotype.bed.gz) whose sample columns already use
# the same "*_hap1" / "*_hap2" naming.
#
# As in GenotypeMeQTL_SAIGE.wdl: per-contig setup (haplotype-splitting, LD
# pruning, sparse GRM, site filtering, splitting into shards) happens in
# one scatter over contigs; per-shard SAIGE runs (RunCisMeQTLChunk, further
# parallelized across CPU cores within each shard) happen in a separate,
# single flat scatter spanning every contig's shards combined - see
# MeQTLTasks.wdl's RunCisMeQTLChunk comment for why it isn't nested
# directly inside the per-contig scatter.
#
# `vcfs`, `methylation_files`, and `contigs` are parallel arrays: index i
# in all three must refer to the same contig. `vcfs` must be phased.

import "MeQTLTasks.wdl"

workflow HaplotypeMeQTL_SAIGE {
    input {
        Array[File] vcfs
        Array[File] methylation_files
        Array[String] contigs
        String prefix

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
        Int n_parallel_workers = 4
        Int max_sites_without_override = 5000
        Boolean allow_large_scan = false

        String plink_docker = "quay.io/biocontainers/plink:1.90b6.21--h031d066_5"
        String bcftools_docker = "quay.io/biocontainers/bcftools:1.19--h8b25389_1"
        String saige_docker = "wzhou88/saige:1.3.6"

        RuntimeAttr? runtime_attr_split_haplotypes
        RuntimeAttr? runtime_attr_ld_prune
        RuntimeAttr? runtime_attr_create_grm
        RuntimeAttr? runtime_attr_filter_sites
        RuntimeAttr? runtime_attr_split_sites
        RuntimeAttr? runtime_attr_run_chunk
        RuntimeAttr? runtime_attr_gather
        RuntimeAttr? runtime_attr_concat_contigs
    }

    Int n_contigs = length(contigs)

    scatter (i in range(n_contigs)) {
        File methylation_file = methylation_files[i]
        String contig = contigs[i]
        String contig_prefix = "~{prefix}.~{contig}"

        call MeQTLTasks.SplitPhasedVcfToHaplotypes {
            input:
                vcf = vcfs[i],
                prefix = contig_prefix,
                docker = bcftools_docker,
                runtime_attr_override = runtime_attr_split_haplotypes
        }

        call MeQTLTasks.LdPruneAndExtract {
            input:
                vcf = SplitPhasedVcfToHaplotypes.haplotype_vcf,
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

        call MeQTLTasks.SplitFilteredSites {
            input:
                filtered_sites = FilterMethylationSites.filtered_sites,
                sites_per_shard = sites_per_shard,
                prefix = contig_prefix,
                docker = saige_docker,
                runtime_attr_override = runtime_attr_split_sites
        }

        Int n_shards_this_contig = length(SplitFilteredSites.shards)

        # Plain value repetition only (no task calls) - broadcasts this
        # contig's per-contig context onto each of its shards, so it can be
        # flattened alongside the shards themselves into flat, per-shard
        # arrays below. See RunCisMeQTLChunk's comment for why the actual
        # per-shard task call must not be nested in this scatter.
        scatter (j in range(n_shards_this_contig)) {
            File broadcast_vcf = SplitPhasedVcfToHaplotypes.haplotype_vcf
            File broadcast_vcf_csi = SplitPhasedVcfToHaplotypes.haplotype_vcf_csi
            File broadcast_pruned_bed = LdPruneAndExtract.bed
            File broadcast_pruned_bim = LdPruneAndExtract.bim
            File broadcast_pruned_fam = LdPruneAndExtract.fam
            File broadcast_sparse_grm = CreateSparseGRM.sparse_grm
            File broadcast_sparse_grm_samples = CreateSparseGRM.sparse_grm_samples
            String broadcast_contig = contig
            String broadcast_prefix = "~{contig_prefix}.shard~{j}"
        }
    }

    Array[File] flat_shard_sites = flatten(SplitFilteredSites.shards)
    Array[File] flat_vcf = flatten(broadcast_vcf)
    Array[File] flat_vcf_csi = flatten(broadcast_vcf_csi)
    Array[File] flat_pruned_bed = flatten(broadcast_pruned_bed)
    Array[File] flat_pruned_bim = flatten(broadcast_pruned_bim)
    Array[File] flat_pruned_fam = flatten(broadcast_pruned_fam)
    Array[File] flat_sparse_grm = flatten(broadcast_sparse_grm)
    Array[File] flat_sparse_grm_samples = flatten(broadcast_sparse_grm_samples)
    Array[String] flat_contig = flatten(broadcast_contig)
    Array[String] flat_prefix = flatten(broadcast_prefix)
    Int n_total_shards = length(flat_shard_sites)

    scatter (k in range(n_total_shards)) {
        call MeQTLTasks.RunCisMeQTLChunk {
            input:
                filtered_sites = flat_shard_sites[k],
                contig = flat_contig[k],
                vcf = flat_vcf[k],
                vcf_csi = flat_vcf_csi[k],
                pruned_bed = flat_pruned_bed[k],
                pruned_bim = flat_pruned_bim[k],
                pruned_fam = flat_pruned_fam[k],
                sparse_grm = flat_sparse_grm[k],
                sparse_grm_samples = flat_sparse_grm_samples[k],
                relatedness_cutoff = relatedness_cutoff,
                cis_window = cis_window,
                covariates_file = covariates_file,
                covar_col_list = covar_col_list,
                qcovar_col_list = qcovar_col_list,
                min_samples_per_site = min_samples_per_site,
                vcf_field = "GT",
                inv_normalize = inv_normalize,
                n_parallel_workers = n_parallel_workers,
                max_sites_without_override = max_sites_without_override,
                allow_large_scan = allow_large_scan,
                prefix = flat_prefix[k],
                docker = saige_docker,
                runtime_attr_override = runtime_attr_run_chunk
        }
    }

    call MeQTLTasks.GatherChunksByContig {
        input:
            chunk_files = RunCisMeQTLChunk.chunk_assoc,
            contigs = flat_contig,
            prefix = "~{prefix}.cis_meQTL.haplotype",
            docker = saige_docker,
            runtime_attr_override = runtime_attr_gather
    }

    call MeQTLTasks.ConcatenateTsvs as ConcatenateAllContigs {
        input:
            tsvs = RunCisMeQTLChunk.chunk_assoc,
            prefix = "~{prefix}.cis_meQTL.haplotype.all_contigs",
            docker = bcftools_docker,
            runtime_attr_override = runtime_attr_concat_contigs
    }

    output {
        Array[File] per_contig_assoc = GatherChunksByContig.per_contig_assoc
        File combined_assoc = ConcatenateAllContigs.merged
        Array[File] per_contig_haplotype_vcf = SplitPhasedVcfToHaplotypes.haplotype_vcf
        Array[File] per_contig_pruned_bed = LdPruneAndExtract.bed
        Array[File] per_contig_sparse_grm = CreateSparseGRM.sparse_grm
        Array[File] per_shard_skipped_logs = RunCisMeQTLChunk.skipped_sites_log
    }
}
