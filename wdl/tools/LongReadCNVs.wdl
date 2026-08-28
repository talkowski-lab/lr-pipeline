version 1.0

import "../utils/LRCNVs.wdl"
import "../utils/DepthPreprocessing.wdl"
import "../utils/DepthClustering.wdl"
import "../utils/GenotypeDepth.wdl"
import "../utils/Structs.wdl"

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
        File? contig_subset_list
        File training_intervals
        File median_coverage
        Int gcnv_qs_cutoff = 30

        String prefix
        String variant_prefix

        String? chr_x
        String? chr_y

        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_annotate_intervals
        RuntimeAttr? runtime_attr_filter_intervals
        RuntimeAttr? runtime_attr_scatter_intervals
        RuntimeAttr? runtime_attr_determine_contig_ploidy
        RuntimeAttr? runtime_attr_germline_cnv_caller
        RuntimeAttr? runtime_attr_postprocess_germline_cnv_calls
        RuntimeAttr? runtime_attr_collect_sample_quality_metrics
        RuntimeAttr? runtime_attr_collect_model_quality_metrics
        RuntimeAttr? runtime_attr_gcnv_vcf_to_bed
        RuntimeAttr? runtime_attr_merge_sample
        RuntimeAttr? runtime_attr_merge_set
        RuntimeAttr? runtime_attr_make_ploidy_table
        RuntimeAttr? runtime_attr_cnv_bed_to_vcf
        RuntimeAttr? runtime_attr_concat_preprocessed_vcfs
        RuntimeAttr? runtime_attr_create_ploidy_table
        RuntimeAttr? runtime_attr_sv_cluster
        RuntimeAttr? runtime_attr_exclude_intervals
        RuntimeAttr? runtime_attr_gatk_to_svtk_vcf
        RuntimeAttr? runtime_attr_concat_clustered_vcfs
        RuntimeAttr? runtime_attr_train_sv_genotyping
        RuntimeAttr? runtime_attr_genotype_svs
        RuntimeAttr? runtime_attr_concat_genotyped_vcfs
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
            gatk_docker = gatk_docker,
            runtime_attr_annotate_intervals = runtime_attr_annotate_intervals,
            runtime_attr_filter_intervals = runtime_attr_filter_intervals,
            runtime_attr_scatter_intervals = runtime_attr_scatter_intervals,
            runtime_attr_determine_contig_ploidy = runtime_attr_determine_contig_ploidy,
            runtime_attr_germline_cnv_caller = runtime_attr_germline_cnv_caller,
            runtime_attr_postprocess_germline_cnv_calls = runtime_attr_postprocess_germline_cnv_calls,
            runtime_attr_collect_sample_quality_metrics = runtime_attr_collect_sample_quality_metrics,
            runtime_attr_collect_model_quality_metrics = runtime_attr_collect_model_quality_metrics
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
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_gcnv_vcf_to_bed = runtime_attr_gcnv_vcf_to_bed,
            runtime_attr_merge_sample = runtime_attr_merge_sample,
            runtime_attr_merge_set = runtime_attr_merge_set,
            runtime_attr_make_ploidy_table = runtime_attr_make_ploidy_table,
            runtime_attr_cnv_bed_to_vcf = runtime_attr_cnv_bed_to_vcf,
            runtime_attr_concat_vcfs = runtime_attr_concat_preprocessed_vcfs
    }

    call DepthClustering.DepthClustering {
        input:
            depth_vcf = DepthPreprocessing.merged_vcf,
            prefix = prefix,
            variant_prefix = variant_prefix,
            pedigree = pedigree,
            contig_list = primary_contigs_list,
            contig_subset_list = contig_subset_list,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            chr_x = chr_x,
            chr_y = chr_y,
            gatk_docker = gatk_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_create_ploidy_table = runtime_attr_create_ploidy_table,
            runtime_attr_sv_cluster = runtime_attr_sv_cluster,
            runtime_attr_exclude_intervals = runtime_attr_exclude_intervals,
            runtime_attr_gatk_to_svtk_vcf = runtime_attr_gatk_to_svtk_vcf,
            runtime_attr_concat_vcfs = runtime_attr_concat_clustered_vcfs
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
            contig_subset_list = contig_subset_list,
            chr_x = chr_x,
            chr_y = chr_y,
            gatk_docker = gatk_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_train_sv_genotyping = runtime_attr_train_sv_genotyping,
            runtime_attr_genotype_svs = runtime_attr_genotype_svs,
            runtime_attr_concat_vcfs = runtime_attr_concat_genotyped_vcfs
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
