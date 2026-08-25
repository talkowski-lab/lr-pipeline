version 1.0

import "LRCNVs.wdl"
import "DepthPreprocessing.wdl"
import "DepthClustering.wdl"
import "GenotypeDepth.wdl"

workflow LongReadCNVs {
    meta {
        description: "Workflow to run GATK-gCNV on long-read samples, cluster the calls and then genotype them in all samples."
    }

    input {
        File intervals
        Array[String]+ sample_ids
        Array[String]+ depth_profiles
        String batch_id
        File contig_ploidy_priors
        Int? num_intervals_per_scatter
        File merged_bincov
        File ref_fa
        File ref_fai
        File ref_dict

        File pedigree
        File primary_contigs_list
        File training_intervals
        File median_coverage
        Int gcnv_qs_cutoff = 30

        String output_prefix
        String variant_prefix

        String? chr_x
        String? chr_y

        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
    }

    call LRCNVs.LRCNVs {
        input:
            intervals = intervals,
            sample_ids = sample_ids,
            depth_profiles = depth_profiles,
            cohort_id = batch_id,
            contig_ploidy_priors = contig_ploidy_priors,
            num_intervals_per_scatter = num_intervals_per_scatter,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            gatk_docker = gatk_docker
    }

    call DepthPreprocessing.DepthPreprocessing {
        input:
            sample_ids = sample_ids,
            genotyped_segments_vcfs = LRCNVs.genotyped_segments_vcfs,
            contig_ploidy_calls_tar = LRCNVs.contig_ploidy_calls_tar,
            primary_contigs_list = primary_contigs_list,
            ref_fai = ref_fai,
            pedigree = pedigree,
            batch_id = batch_id,
            chr_x = chr_x,
            chr_y = chr_y,
            gcnv_qs_cutoff = gcnv_qs_cutoff,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_docker = sv_pipeline_docker
    }

    call DepthClustering.DepthClustering {
        input:
            depth_vcf = DepthPreprocessing.merged_vcf,
            output_prefix = output_prefix,
            variant_prefix = variant_prefix,
            pedigree = pedigree,
            contig_list = primary_contigs_list,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            chr_x = chr_x,
            chr_y = chr_y,
            gatk_docker = gatk_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_docker = sv_pipeline_docker
    }

    call GenotypeDepth.GenotypeDepth {
        input:
           batch_id = batch_id,
           vcf = DepthClustering.clustered_vcf,
           training_intervals = training_intervals,
           median_coverage = median_coverage,
           rd_file = merged_bincov,
           ref_dict = ref_dict,
           ploidy_table = DepthClustering.ploidy_table,
           contig_list = primary_contigs_list,
           chr_x = chr_x,
           chr_y = chr_y,
           gatk_docker = gatk_docker,
           sv_base_mini_docker = sv_base_mini_docker
    }

    output {
        File merged_cnvs_vcf = DepthPreprocessing.merged_vcf
        File merged_cnvs_vcf_index = DepthPreprocessing.merged_vcf_index
        File ploidy_table = DepthPreprocessing.ploidy_table
        File genotyped_depth_vcf = GenotypeDepth.genotyped_depth_vcf
        File genotyped_depth_vcf_index = GenotypeDepth.genotyped_depth_vcf_index
        File genotyping_rd_table = GenotypeDepth.genotyping_rd_table
    }
}
