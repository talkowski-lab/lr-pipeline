version 1.0

import "../utils/LRCNVs.wdl"
import "../utils/DepthPreprocessing.wdl"
import "../utils/DepthClustering.wdl"
import "../utils/GenotypeDepth.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow LongReadCNVs {
    meta {
        description: "Workflow to run GATK-gCNV on long-read samples, cluster the calls and then genotype them in all samples."
    }

    input {
        File intervals
        Array[String]+ sample_ids
        Array[File]+ depth_profiles
        Boolean sort_depth_profiles
        String batch_id
        File contig_ploidy_priors
        Int num_intervals_per_scatter = 10000
        File merged_bincov
        File merged_bincov_idx
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

        String chr_x = "chrX"
        String chr_y = "chrY"

        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        File? gatk4_jar_override
        File? mappability_track_bed
        File? mappability_track_bed_idx
        File? segmental_duplication_track_bed
        File? segmental_duplication_track_bed_idx
        Int? feature_query_lookahead
        File? blacklist_intervals
        Int? low_count_filter_count_threshold
        Float? low_count_filter_percentage_of_samples
        Float? extreme_count_filter_minimum_percentile
        Float? extreme_count_filter_maximum_percentile
        Float? extreme_count_filter_percentage_of_samples
        Float? ploidy_mean_bias_standard_deviation
        Float? ploidy_mapping_error_rate
        Float? ploidy_global_psi_scale
        Float? ploidy_sample_psi_scale
        Float? gcnv_p_alt
        Float? gcnv_p_active
        Float? gcnv_cnv_coherence_length
        Float? gcnv_class_coherence_length
        Int? gcnv_max_copy_number
        Int? gcnv_max_bias_factors
        Float? gcnv_mapping_error_rate
        Float? gcnv_interval_psi_scale
        Float? gcnv_sample_psi_scale
        Float? gcnv_depth_correction_tau
        Float? gcnv_log_mean_bias_standard_deviation
        Float? gcnv_init_ard_rel_unexplained_variance
        Int? gcnv_num_gc_bins
        Float? gcnv_gc_curve_standard_deviation
        String? gcnv_copy_number_posterior_expectation_mode
        Boolean? gcnv_enable_bias_factors
        Int? gcnv_active_class_padding_hybrid_mode
        Float? gcnv_learning_rate
        Float? gcnv_adamax_beta_1
        Float? gcnv_adamax_beta_2
        Int? gcnv_log_emission_samples_per_round
        Float? gcnv_log_emission_sampling_median_rel_error
        Int? gcnv_log_emission_sampling_rounds
        Int? gcnv_max_advi_iter_first_epoch
        Int? gcnv_max_advi_iter_subsequent_epochs
        Int? gcnv_min_training_epochs
        Int? gcnv_max_training_epochs
        Float? gcnv_initial_temperature
        Int? gcnv_num_thermal_advi_iters
        Int? gcnv_convergence_snr_averaging_window
        Float? gcnv_convergence_snr_trigger_threshold
        Int? gcnv_convergence_snr_countdown_window
        Int? gcnv_max_calling_iters
        Float? gcnv_caller_update_convergence_threshold
        Float? gcnv_caller_internal_admixing_rate
        Float? gcnv_caller_external_admixing_rate
        Boolean? gcnv_disable_annealing
        Int ref_copy_number_autosomal_contigs = 2
        Array[String]? allosomal_contigs
        Int maximum_number_events_per_sample = 1000
        Float? defragment_max_dist
        Boolean fast_mode = true
        String clustering_algorithm = "SINGLE_LINKAGE"
        Boolean? enable_cnv
        Boolean? default_no_call
        Boolean? omit_members
        String? breakpoint_summary_strategy
        Float? defrag_padding_fraction
        Float? defrag_sample_overlap
        Float depth_sample_overlap = 0
        Float depth_interval_overlap = 0.8
        Float? depth_size_similarity
        Int depth_breakend_window = 10000000
        File? exclude_intervals
        Float exclude_overlap_fraction = 0.5
        File? gatk_to_svtk_script
        Boolean svtk_set_pass = false

        RuntimeAttr? runtime_attr_sort_depth_profiles
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
        RuntimeAttr? runtime_attr_sv_cluster
        RuntimeAttr? runtime_attr_exclude_intervals
        RuntimeAttr? runtime_attr_gatk_to_svtk_vcf
        RuntimeAttr? runtime_attr_concat_clustered_vcfs
        RuntimeAttr? runtime_attr_train_sv_genotyping
        RuntimeAttr? runtime_attr_genotype_svs
        RuntimeAttr? runtime_attr_concat_genotyped_vcfs
    }

    if (sort_depth_profiles) {
        scatter (i in range(length(depth_profiles))) {
            call Helpers.SortReadCounts {
                input:
                    read_counts = depth_profiles[i],
                    prefix = sample_ids[i] + ".sorted_counts",
                    docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_sort_depth_profiles
            }
        }
    }

    Array[File]+ depth_profiles_ = select_first([SortReadCounts.sorted_read_counts, depth_profiles])

    call LRCNVs.LRCNVs {
        input:
            intervals = intervals,
            sample_ids = sample_ids,
            depth_profiles = depth_profiles_,
            prefix = prefix,
            cohort_id = batch_id,
            contig_ploidy_priors = contig_ploidy_priors,
            num_intervals_per_scatter = num_intervals_per_scatter,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            gatk_docker = gatk_docker,
            gatk4_jar_override = gatk4_jar_override,
            mappability_track_bed = mappability_track_bed,
            mappability_track_bed_idx = mappability_track_bed_idx,
            segmental_duplication_track_bed = segmental_duplication_track_bed,
            segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
            feature_query_lookahead = feature_query_lookahead,
            blacklist_intervals = blacklist_intervals,
            low_count_filter_count_threshold = low_count_filter_count_threshold,
            low_count_filter_percentage_of_samples = low_count_filter_percentage_of_samples,
            extreme_count_filter_minimum_percentile = extreme_count_filter_minimum_percentile,
            extreme_count_filter_maximum_percentile = extreme_count_filter_maximum_percentile,
            extreme_count_filter_percentage_of_samples = extreme_count_filter_percentage_of_samples,
            ploidy_mean_bias_standard_deviation = ploidy_mean_bias_standard_deviation,
            ploidy_mapping_error_rate = ploidy_mapping_error_rate,
            ploidy_global_psi_scale = ploidy_global_psi_scale,
            ploidy_sample_psi_scale = ploidy_sample_psi_scale,
            gcnv_p_alt = gcnv_p_alt,
            gcnv_p_active = gcnv_p_active,
            gcnv_cnv_coherence_length = gcnv_cnv_coherence_length,
            gcnv_class_coherence_length = gcnv_class_coherence_length,
            gcnv_max_copy_number = gcnv_max_copy_number,
            gcnv_max_bias_factors = gcnv_max_bias_factors,
            gcnv_mapping_error_rate = gcnv_mapping_error_rate,
            gcnv_interval_psi_scale = gcnv_interval_psi_scale,
            gcnv_sample_psi_scale = gcnv_sample_psi_scale,
            gcnv_depth_correction_tau = gcnv_depth_correction_tau,
            gcnv_log_mean_bias_standard_deviation = gcnv_log_mean_bias_standard_deviation,
            gcnv_init_ard_rel_unexplained_variance = gcnv_init_ard_rel_unexplained_variance,
            gcnv_num_gc_bins = gcnv_num_gc_bins,
            gcnv_gc_curve_standard_deviation = gcnv_gc_curve_standard_deviation,
            gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
            gcnv_enable_bias_factors = gcnv_enable_bias_factors,
            gcnv_active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
            gcnv_learning_rate = gcnv_learning_rate,
            gcnv_adamax_beta_1 = gcnv_adamax_beta_1,
            gcnv_adamax_beta_2 = gcnv_adamax_beta_2,
            gcnv_log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
            gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
            gcnv_log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
            gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
            gcnv_max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
            gcnv_min_training_epochs = gcnv_min_training_epochs,
            gcnv_max_training_epochs = gcnv_max_training_epochs,
            gcnv_initial_temperature = gcnv_initial_temperature,
            gcnv_num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
            gcnv_convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
            gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
            gcnv_convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
            gcnv_max_calling_iters = gcnv_max_calling_iters,
            gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
            gcnv_caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
            gcnv_caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
            gcnv_disable_annealing = gcnv_disable_annealing,
            ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
            allosomal_contigs = allosomal_contigs,
            maximum_number_events_per_sample = maximum_number_events_per_sample,
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
            genotyped_segments_vcf_idxs = LRCNVs.genotyped_segments_vcf_idxs,
            contig_ploidy_calls_tar = LRCNVs.contig_ploidy_calls_tar,
            primary_contigs_list = primary_contigs_list,
            ref_fai = ref_fai,
            pedigree = pedigree,
            batch_id = batch_id,
            prefix = prefix,
            chr_x = chr_x,
            chr_y = chr_y,
            gcnv_qs_cutoff = gcnv_qs_cutoff,
            defragment_max_dist = defragment_max_dist,
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
            depth_vcf_idx = DepthPreprocessing.merged_vcf_idx,
            ploidy_table = DepthPreprocessing.ploidy_table,
            prefix = prefix,
            variant_prefix = variant_prefix,
            contig_list = primary_contigs_list,
            contig_subset_list = contig_subset_list,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            gatk_docker = gatk_docker,
            sv_base_mini_docker = sv_base_mini_docker,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_sv_cluster = runtime_attr_sv_cluster,
            runtime_attr_exclude_intervals = runtime_attr_exclude_intervals,
            runtime_attr_gatk_to_svtk_vcf = runtime_attr_gatk_to_svtk_vcf,
            runtime_attr_concat_vcfs = runtime_attr_concat_clustered_vcfs,
            fast_mode = fast_mode,
            clustering_algorithm = clustering_algorithm,
            enable_cnv = enable_cnv,
            default_no_call = default_no_call,
            omit_members = omit_members,
            breakpoint_summary_strategy = breakpoint_summary_strategy,
            defrag_padding_fraction = defrag_padding_fraction,
            defrag_sample_overlap = defrag_sample_overlap,
            depth_sample_overlap = depth_sample_overlap,
            depth_interval_overlap = depth_interval_overlap,
            depth_size_similarity = depth_size_similarity,
            depth_breakend_window = depth_breakend_window,
            exclude_intervals = exclude_intervals,
            exclude_overlap_fraction = exclude_overlap_fraction,
            gatk_to_svtk_script = gatk_to_svtk_script,
            svtk_set_pass = svtk_set_pass
    }

    call GenotypeDepth.GenotypeDepth {
        input:
            prefix = prefix,
            vcf = DepthClustering.clustered_vcf,
            vcf_idx = DepthClustering.clustered_vcf_idx,
            training_intervals = training_intervals,
            median_coverage = median_coverage,
            rd_file = merged_bincov,
            rd_file_idx = merged_bincov_idx,
            ref_dict = ref_dict,
            ploidy_table = DepthPreprocessing.ploidy_table,
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
        File merged_cnvs_vcf_idx = DepthPreprocessing.merged_vcf_idx
        File ploidy_table = DepthPreprocessing.ploidy_table
        File genotyped_depth_vcf = GenotypeDepth.genotyped_depth_vcf
        File genotyped_depth_vcf_idx = GenotypeDepth.genotyped_depth_vcf_idx
        File genotyping_rd_table = GenotypeDepth.genotyping_rd_table
    }
}
