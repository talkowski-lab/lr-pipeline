# BSD 3-Clause License
#
# Copyright (c) 2019, Broad Institute
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

version 1.0

import "Helpers.wdl"
import "Structs.wdl"

workflow LRCNVs {
    meta {
        description: [
            "This component calls copy-number variants across a cohort using GATK germline CNV (gCNV) cohort mode. From per-sample depth profiles over a shared interval list it annotates and filters intervals, determines contig ploidy, fits gCNV across scattered interval shards, post-processes per-sample calls into genotyped interval and segment VCFs, and collects sample- and model-level QC.",
            "Setting `num_training_samples` below the cohort size switches the component to a hybrid case-cohort mode: that many samples are drawn at random and used to fit the contig-ploidy and gCNV models in cohort mode, and every remaining sample is then called against those models in case mode. Interval annotation, interval filtering and the fitted models therefore derive from the training samples alone, while the per-sample genotyped VCFs, denoised copy ratios, sample QC and contig-ploidy calls still cover the whole cohort in `sample_ids` order."
        ]
    }

    parameter_meta {
        intervals: "Interval list over which CNVs are called."
        sample_ids: "Sample IDs in the cohort."
        depth_profiles: "Per-sample read-depth profiles, aligned to `sample_ids`."
        cohort_id: "Identifier for the cohort."
        contig_ploidy_priors: "Contig ploidy priors used to determine per-sample contig ploidy."
        ref_fa: "From references."
        ref_fai: "From references."
        ref_dict: "From references."
        num_intervals_per_scatter: "Number of intervals processed per gCNV scatter shard. GermlineCNVCaller memory grows with samples times intervals per shard, so raising this above the default needs more memory in `runtime_attr_germline_cnv_caller`."
        num_training_samples: "Number of samples drawn at random to fit the contig-ploidy and gCNV models in cohort mode, with every remaining sample called against those models in case mode. Left unset, or set to at least the cohort size, every sample is called in cohort mode. Interval filtering percentages then apply over the training samples only, so a training set of fewer than a few dozen samples degrades the fitted models."
        subsample_seed: "Random seed used to draw the training samples."
        gatk4_jar_override: "Override GATK4 jar."
        mappability_track_bed: "Mappability track used to annotate intervals."
        mappability_track_bed_idx: "Index for `mappability_track_bed`."
        segmental_duplication_track_bed: "Segmental-duplication track used to annotate intervals."
        segmental_duplication_track_bed_idx: "Index for `segmental_duplication_track_bed`."
        feature_query_lookahead: "Base pairs to look ahead when querying interval-annotation feature tracks."
        blacklist_intervals: "Intervals to exclude from calling."
        low_count_filter_count_threshold: "Minimum read count for an interval to be considered well-covered in a sample."
        low_count_filter_percentage_of_samples: "Minimum percentage of samples that must meet `low_count_filter_count_threshold` for an interval to pass."
        extreme_count_filter_minimum_percentile: "Lower count percentile below which an interval is considered an outlier."
        extreme_count_filter_maximum_percentile: "Upper count percentile above which an interval is considered an outlier."
        extreme_count_filter_percentage_of_samples: "Minimum percentage of samples that must pass the extreme-count percentile bounds for an interval to pass."
        ploidy_mean_bias_standard_deviation: "DetermineGermlineContigPloidy --mean-bias-standard-deviation."
        ploidy_mapping_error_rate: "DetermineGermlineContigPloidy --mapping-error-rate."
        ploidy_global_psi_scale: "DetermineGermlineContigPloidy --global-psi-scale."
        ploidy_sample_psi_scale: "DetermineGermlineContigPloidy --sample-psi-scale."
        gcnv_p_alt: "GermlineCNVCaller --p-alt."
        gcnv_p_active: "GermlineCNVCaller --p-active."
        gcnv_cnv_coherence_length: "GermlineCNVCaller --cnv-coherence-length."
        gcnv_class_coherence_length: "GermlineCNVCaller --class-coherence-length."
        gcnv_max_copy_number: "GermlineCNVCaller --max-copy-number."
        gcnv_max_bias_factors: "GermlineCNVCaller --max-bias-factors."
        gcnv_mapping_error_rate: "GermlineCNVCaller --mapping-error-rate."
        gcnv_interval_psi_scale: "GermlineCNVCaller --interval-psi-scale."
        gcnv_sample_psi_scale: "GermlineCNVCaller --sample-psi-scale."
        gcnv_depth_correction_tau: "GermlineCNVCaller --depth-correction-tau."
        gcnv_log_mean_bias_standard_deviation: "GermlineCNVCaller --log-mean-bias-standard-deviation."
        gcnv_init_ard_rel_unexplained_variance: "GermlineCNVCaller --init-ard-rel-unexplained-variance."
        gcnv_num_gc_bins: "GermlineCNVCaller --num-gc-bins."
        gcnv_gc_curve_standard_deviation: "GermlineCNVCaller --gc-curve-standard-deviation."
        gcnv_copy_number_posterior_expectation_mode: "GermlineCNVCaller --copy-number-posterior-expectation-mode."
        gcnv_enable_bias_factors: "GermlineCNVCaller --enable-bias-factors."
        gcnv_active_class_padding_hybrid_mode: "GermlineCNVCaller --active-class-padding-hybrid-mode."
        gcnv_learning_rate: "GermlineCNVCaller --learning-rate."
        gcnv_adamax_beta_1: "GermlineCNVCaller --adamax-beta-1."
        gcnv_adamax_beta_2: "GermlineCNVCaller --adamax-beta-2."
        gcnv_log_emission_samples_per_round: "GermlineCNVCaller --log-emission-samples-per-round."
        gcnv_log_emission_sampling_median_rel_error: "GermlineCNVCaller --log-emission-sampling-median-rel-error."
        gcnv_log_emission_sampling_rounds: "GermlineCNVCaller --log-emission-sampling-rounds."
        gcnv_max_advi_iter_first_epoch: "GermlineCNVCaller --max-advi-iter-first-epoch."
        gcnv_max_advi_iter_subsequent_epochs: "GermlineCNVCaller --max-advi-iter-subsequent-epochs."
        gcnv_min_training_epochs: "GermlineCNVCaller --min-training-epochs."
        gcnv_max_training_epochs: "GermlineCNVCaller --max-training-epochs."
        gcnv_initial_temperature: "GermlineCNVCaller --initial-temperature."
        gcnv_num_thermal_advi_iters: "GermlineCNVCaller --num-thermal-advi-iters."
        gcnv_convergence_snr_averaging_window: "GermlineCNVCaller --convergence-snr-averaging-window."
        gcnv_convergence_snr_trigger_threshold: "GermlineCNVCaller --convergence-snr-trigger-threshold."
        gcnv_convergence_snr_countdown_window: "GermlineCNVCaller --convergence-snr-countdown-window."
        gcnv_max_calling_iters: "GermlineCNVCaller --max-calling-iters."
        gcnv_caller_update_convergence_threshold: "GermlineCNVCaller --caller-update-convergence-threshold."
        gcnv_caller_internal_admixing_rate: "GermlineCNVCaller --caller-internal-admixing-rate."
        gcnv_caller_external_admixing_rate: "GermlineCNVCaller --caller-external-admixing-rate."
        gcnv_disable_annealing: "GermlineCNVCaller --disable-annealing."
        ref_copy_number_autosomal_contigs: "Reference copy number for autosomes."
        allosomal_contigs: "Contigs treated as allosomal."
        maximum_number_events_per_sample: "Maximum number of events permitted per sample."
        annotated_intervals: "Intervals annotated with GC content and tracks."
        filtered_intervals: "Intervals retained after filtering the training samples."
        contig_ploidy_model_tar: "Contig-ploidy model fitted on the training samples."
        contig_ploidy_calls_tar: "Per-sample contig-ploidy calls for every sample, ordered as `sample_ids`."
        gcnv_model_tars: "gCNV models fitted on the training samples, one per scatter shard."
        gcnv_calls_tars: "Per-shard gCNV calls for the training samples."
        gcnv_tracking_tars: "Per-shard model-fitting tracking files for the training samples."
        genotyped_intervals_vcfs: "Per-sample genotyped interval VCFs."
        genotyped_intervals_vcf_idxs: "Indexes for `genotyped_intervals_vcfs`."
        genotyped_segments_vcfs: "Per-sample genotyped segment VCFs."
        genotyped_segments_vcf_idxs: "Indexes for `genotyped_segments_vcfs`."
        sample_qc_status_files: "Per-sample QC status files."
        sample_qc_status_strings: "Per-sample QC status strings."
        model_qc_status_file: "Model-level QC status file for the models fitted on the training samples."
        model_qc_string: "Model-level QC status string for the models fitted on the training samples."
        denoised_copy_ratios: "Per-sample denoised copy ratios."
    }

    input {
        File intervals
        Array[String]+ sample_ids
        Array[File]+ depth_profiles
        String prefix
        String cohort_id

        File contig_ploidy_priors
        File ref_fa
        File ref_fai
        File ref_dict
        String gatk_docker
        String sv_pipeline_docker

        Int num_intervals_per_scatter
        Int? num_training_samples
        Int subsample_seed = 42

        File? gatk4_jar_override

        # AnnotateIntervals
        File? mappability_track_bed
        File? mappability_track_bed_idx
        File? segmental_duplication_track_bed
        File? segmental_duplication_track_bed_idx
        Int? feature_query_lookahead

        # FilterIntervals
        File? blacklist_intervals
        Int? low_count_filter_count_threshold
        Float? low_count_filter_percentage_of_samples
        Float? extreme_count_filter_minimum_percentile
        Float? extreme_count_filter_maximum_percentile
        Float? extreme_count_filter_percentage_of_samples

        # DeterminGermlineContigPloidyCohortMode
        Float? ploidy_mean_bias_standard_deviation
        Float? ploidy_mapping_error_rate
        Float? ploidy_global_psi_scale
        Float? ploidy_sample_psi_scale

        # GermlineCNVCallerCohortMode
        Float? gcnv_p_alt
        Float? gcnv_p_active
        Float? gcnv_cnv_coherence_length
        Float? gcnv_class_coherence_length
        Int? gcnv_max_copy_number

        # GermlineCNVCallerCohortMode - germline CNV denoising model
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

        # GermlineCNVCallerCohortMode - Hybrid ADVI
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

        # PostprocessGermlineCNVCalls
        Int ref_copy_number_autosomal_contigs = 2
        Array[String]? allosomal_contigs

        # CollectSampleQualityMetrics
        Int maximum_number_events_per_sample = 1000

        RuntimeAttr? runtime_attr_subsample_indices
        RuntimeAttr? runtime_attr_annotate_intervals
        RuntimeAttr? runtime_attr_filter_intervals
        RuntimeAttr? runtime_attr_scatter_intervals
        RuntimeAttr? runtime_attr_determine_contig_ploidy
        RuntimeAttr? runtime_attr_determine_contig_ploidy_case
        RuntimeAttr? runtime_attr_germline_cnv_caller
        RuntimeAttr? runtime_attr_germline_cnv_caller_case
        RuntimeAttr? runtime_attr_postprocess_germline_cnv_calls
        RuntimeAttr? runtime_attr_collect_sample_quality_metrics
        RuntimeAttr? runtime_attr_collect_model_quality_metrics
        RuntimeAttr? runtime_attr_merge_contig_ploidy_calls
    }

    Int num_training_samples_ = select_first([num_training_samples, length(sample_ids)])
    Boolean run_case_mode = num_training_samples_ < length(sample_ids)

    if (run_case_mode) {
        call Helpers.SubsampleIndices {
            input:
                num_items = length(sample_ids),
                num_subsampled = num_training_samples_,
                seed = subsample_seed,
                prefix = prefix,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_subsample_indices
        }
    }

    Array[Int] training_indices = select_first([SubsampleIndices.subsampled_indices, range(length(sample_ids))])

    scatter (training_index in training_indices) {
        String training_sample_ids = sample_ids[training_index]
        File training_depth_profiles = depth_profiles[training_index]
    }

    call AnnotateIntervals {
        input:
            intervals = intervals,
            prefix = prefix,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            mappability_track_bed = mappability_track_bed,
            mappability_track_bed_idx = mappability_track_bed_idx,
            segmental_duplication_track_bed = segmental_duplication_track_bed,
            segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
            feature_query_lookahead = feature_query_lookahead,
            gatk4_jar_override = gatk4_jar_override,
            docker = gatk_docker,
            runtime_attr_override = runtime_attr_annotate_intervals
    }

    call FilterIntervals {
        input:
            intervals = intervals,
            prefix = prefix,
            annotated_intervals = AnnotateIntervals.annotated_intervals,
            blacklist_intervals = blacklist_intervals,
            read_count_files = training_depth_profiles,
            low_count_filter_count_threshold = low_count_filter_count_threshold,
            low_count_filter_percentage_of_samples = low_count_filter_percentage_of_samples,
            extreme_count_filter_minimum_percentile = extreme_count_filter_minimum_percentile,
            extreme_count_filter_maximum_percentile = extreme_count_filter_maximum_percentile,
            extreme_count_filter_percentage_of_samples = extreme_count_filter_percentage_of_samples,
            gatk4_jar_override = gatk4_jar_override,
            docker = gatk_docker,
            runtime_attr_override = runtime_attr_filter_intervals
    }

    call DetermineGermlineContigPloidyCohortMode {
        input:
            cohort_id = cohort_id,
            prefix = prefix,
            intervals = FilterIntervals.filtered_intervals,
            read_count_files = training_depth_profiles,
            contig_ploidy_priors = contig_ploidy_priors,
            gatk4_jar_override = gatk4_jar_override,
            docker = gatk_docker,
            mean_bias_standard_deviation = ploidy_mean_bias_standard_deviation,
            mapping_error_rate = ploidy_mapping_error_rate,
            global_psi_scale = ploidy_global_psi_scale,
            sample_psi_scale = ploidy_sample_psi_scale,
            runtime_attr_override = runtime_attr_determine_contig_ploidy
    }

    call ScatterIntervals {
        input:
            interval_list = FilterIntervals.filtered_intervals,
            prefix = prefix,
            num_intervals_per_scatter = num_intervals_per_scatter,
            docker = gatk_docker,
            runtime_attr_override = runtime_attr_scatter_intervals
    }

    scatter (scatter_index in range(length(ScatterIntervals.scattered_interval_lists))) {
        call GermlineCNVCallerCohortMode {
            input:
                scatter_index = scatter_index,
                cohort_id = cohort_id,
                prefix = prefix,
                read_count_files = training_depth_profiles,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar,
                intervals = ScatterIntervals.scattered_interval_lists[scatter_index],
                annotated_intervals = AnnotateIntervals.annotated_intervals,
                gatk4_jar_override = gatk4_jar_override,
                docker = gatk_docker,
                p_alt = gcnv_p_alt,
                p_active = gcnv_p_active,
                cnv_coherence_length = gcnv_cnv_coherence_length,
                class_coherence_length = gcnv_class_coherence_length,
                max_copy_number = gcnv_max_copy_number,
                max_bias_factors = gcnv_max_bias_factors,
                mapping_error_rate = gcnv_mapping_error_rate,
                interval_psi_scale = gcnv_interval_psi_scale,
                sample_psi_scale = gcnv_sample_psi_scale,
                depth_correction_tau = gcnv_depth_correction_tau,
                log_mean_bias_standard_deviation = gcnv_log_mean_bias_standard_deviation,
                init_ard_rel_unexplained_variance = gcnv_init_ard_rel_unexplained_variance,
                num_gc_bins = gcnv_num_gc_bins,
                gc_curve_standard_deviation = gcnv_gc_curve_standard_deviation,
                copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                enable_bias_factors = gcnv_enable_bias_factors,
                active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                learning_rate = gcnv_learning_rate,
                adamax_beta_1 = gcnv_adamax_beta_1,
                adamax_beta_2 = gcnv_adamax_beta_2,
                log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                min_training_epochs = gcnv_min_training_epochs,
                max_training_epochs = gcnv_max_training_epochs,
                initial_temperature = gcnv_initial_temperature,
                num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                max_calling_iters = gcnv_max_calling_iters,
                caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                disable_annealing = gcnv_disable_annealing,
                runtime_attr_override = runtime_attr_germline_cnv_caller
        }
    }

    Array[Array[File]] call_tars_sample_by_shard = transpose(GermlineCNVCallerCohortMode.gcnv_call_tars)

    scatter (sample_index in range(length(training_indices))) {
        call PostprocessGermlineCNVCalls {
            input:
                prefix = prefix + "." + training_sample_ids[sample_index],
                gcnv_calls_tars = call_tars_sample_by_shard[sample_index],
                gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar,
                calling_configs = GermlineCNVCallerCohortMode.calling_config_json,
                denoising_configs = GermlineCNVCallerCohortMode.denoising_config_json,
                gcnvkernel_version = GermlineCNVCallerCohortMode.gcnvkernel_version_json,
                sharded_interval_lists = GermlineCNVCallerCohortMode.sharded_interval_list,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar,
                allosomal_contigs = allosomal_contigs,
                ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                sample_index = sample_index,
                gatk4_jar_override = gatk4_jar_override,
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_postprocess_germline_cnv_calls
        }

        call CollectSampleQualityMetrics {
            input:
                genotyped_segments_vcf = PostprocessGermlineCNVCalls.genotyped_segments_vcf,
                prefix = prefix + "." + training_sample_ids[sample_index],
                maximum_number_events = maximum_number_events_per_sample,
                docker = gatk_docker,
                runtime_attr_override = runtime_attr_collect_sample_quality_metrics
        }
    }

    call CollectModelQualityMetrics {
        input:
            gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar,
            prefix = prefix,
            docker = gatk_docker,
            runtime_attr_override = runtime_attr_collect_model_quality_metrics
    }

    if (run_case_mode) {
        Array[Int] case_indices = select_first([SubsampleIndices.remaining_indices])

        scatter (case_index in case_indices) {
            String case_sample_ids = sample_ids[case_index]
            File case_depth_profiles = depth_profiles[case_index]
        }

        call DetermineGermlineContigPloidyCaseMode {
            input:
                prefix = prefix,
                read_count_files = case_depth_profiles,
                contig_ploidy_model_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_model_tar,
                gatk4_jar_override = gatk4_jar_override,
                docker = gatk_docker,
                mapping_error_rate = ploidy_mapping_error_rate,
                sample_psi_scale = ploidy_sample_psi_scale,
                runtime_attr_override = runtime_attr_determine_contig_ploidy_case
        }

        scatter (scatter_index in range(length(GermlineCNVCallerCohortMode.gcnv_model_tar))) {
            call GermlineCNVCallerCaseMode {
                input:
                    scatter_index = scatter_index,
                    prefix = prefix,
                    sample_ids = case_sample_ids,
                    read_count_files = case_depth_profiles,
                    contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                    gcnv_model_tar = GermlineCNVCallerCohortMode.gcnv_model_tar[scatter_index],
                    gatk4_jar_override = gatk4_jar_override,
                    docker = gatk_docker,
                    p_alt = gcnv_p_alt,
                    cnv_coherence_length = gcnv_cnv_coherence_length,
                    max_copy_number = gcnv_max_copy_number,
                    mapping_error_rate = gcnv_mapping_error_rate,
                    sample_psi_scale = gcnv_sample_psi_scale,
                    depth_correction_tau = gcnv_depth_correction_tau,
                    copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                    active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                    learning_rate = gcnv_learning_rate,
                    adamax_beta_1 = gcnv_adamax_beta_1,
                    adamax_beta_2 = gcnv_adamax_beta_2,
                    log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                    log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                    log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                    max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                    max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                    min_training_epochs = gcnv_min_training_epochs,
                    max_training_epochs = gcnv_max_training_epochs,
                    initial_temperature = gcnv_initial_temperature,
                    num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                    convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                    convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                    convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                    max_calling_iters = gcnv_max_calling_iters,
                    caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                    caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                    caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                    disable_annealing = gcnv_disable_annealing,
                    runtime_attr_override = runtime_attr_germline_cnv_caller_case
            }
        }

        Array[Array[File]] case_call_tars_sample_by_shard = transpose(GermlineCNVCallerCaseMode.gcnv_call_tars)

        scatter (case_sample_index in range(length(case_indices))) {
            call PostprocessGermlineCNVCalls as PostprocessGermlineCNVCallsCase {
                input:
                    prefix = prefix + "." + case_sample_ids[case_sample_index],
                    gcnv_calls_tars = case_call_tars_sample_by_shard[case_sample_index],
                    gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar,
                    calling_configs = GermlineCNVCallerCaseMode.calling_config_json,
                    denoising_configs = GermlineCNVCallerCaseMode.denoising_config_json,
                    gcnvkernel_version = GermlineCNVCallerCaseMode.gcnvkernel_version_json,
                    sharded_interval_lists = GermlineCNVCallerCaseMode.sharded_interval_list,
                    contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                    allosomal_contigs = allosomal_contigs,
                    ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                    sample_index = case_sample_index,
                    gatk4_jar_override = gatk4_jar_override,
                    docker = gatk_docker,
                    runtime_attr_override = runtime_attr_postprocess_germline_cnv_calls
            }

            call CollectSampleQualityMetrics as CollectSampleQualityMetricsCase {
                input:
                    genotyped_segments_vcf = PostprocessGermlineCNVCallsCase.genotyped_segments_vcf,
                    prefix = prefix + "." + case_sample_ids[case_sample_index],
                    maximum_number_events = maximum_number_events_per_sample,
                    docker = gatk_docker,
                    runtime_attr_override = runtime_attr_collect_sample_quality_metrics
            }
        }

        call MergeContigPloidyCalls {
            input:
                prefix = prefix,
                training_contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar,
                case_contig_ploidy_calls_tar = DetermineGermlineContigPloidyCaseMode.contig_ploidy_calls_tar,
                training_indices = training_indices,
                case_indices = case_indices,
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_merge_contig_ploidy_calls
        }
    }

    Array[Int] concatenated_positions = select_first([SubsampleIndices.concatenated_positions, range(length(sample_ids))])
    Array[File] concatenated_genotyped_intervals_vcfs = flatten([PostprocessGermlineCNVCalls.genotyped_intervals_vcf, select_first([PostprocessGermlineCNVCallsCase.genotyped_intervals_vcf, []])])
    Array[File] concatenated_genotyped_intervals_vcf_idxs = flatten([PostprocessGermlineCNVCalls.genotyped_intervals_vcf_idx, select_first([PostprocessGermlineCNVCallsCase.genotyped_intervals_vcf_idx, []])])
    Array[File] concatenated_genotyped_segments_vcfs = flatten([PostprocessGermlineCNVCalls.genotyped_segments_vcf, select_first([PostprocessGermlineCNVCallsCase.genotyped_segments_vcf, []])])
    Array[File] concatenated_genotyped_segments_vcf_idxs = flatten([PostprocessGermlineCNVCalls.genotyped_segments_vcf_idx, select_first([PostprocessGermlineCNVCallsCase.genotyped_segments_vcf_idx, []])])
    Array[File] concatenated_denoised_copy_ratios = flatten([PostprocessGermlineCNVCalls.denoised_copy_ratios, select_first([PostprocessGermlineCNVCallsCase.denoised_copy_ratios, []])])
    Array[File] concatenated_qc_status_files = flatten([CollectSampleQualityMetrics.qc_status_file, select_first([CollectSampleQualityMetricsCase.qc_status_file, []])])
    Array[String] concatenated_qc_status_strings = flatten([CollectSampleQualityMetrics.qc_status_string, select_first([CollectSampleQualityMetricsCase.qc_status_string, []])])

    scatter (output_index in range(length(sample_ids))) {
        Int concatenated_position = concatenated_positions[output_index]
        File ordered_genotyped_intervals_vcf = concatenated_genotyped_intervals_vcfs[concatenated_position]
        File ordered_genotyped_intervals_vcf_idx = concatenated_genotyped_intervals_vcf_idxs[concatenated_position]
        File ordered_genotyped_segments_vcf = concatenated_genotyped_segments_vcfs[concatenated_position]
        File ordered_genotyped_segments_vcf_idx = concatenated_genotyped_segments_vcf_idxs[concatenated_position]
        File ordered_denoised_copy_ratios = concatenated_denoised_copy_ratios[concatenated_position]
        File ordered_qc_status_file = concatenated_qc_status_files[concatenated_position]
        String ordered_qc_status_string = concatenated_qc_status_strings[concatenated_position]
    }

    output {
        File annotated_intervals = AnnotateIntervals.annotated_intervals
        File filtered_intervals = FilterIntervals.filtered_intervals
        File contig_ploidy_model_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_model_tar
        File contig_ploidy_calls_tar = select_first([MergeContigPloidyCalls.contig_ploidy_calls_tar, DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar])
        Array[File] gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar
        Array[Array[File]] gcnv_calls_tars = GermlineCNVCallerCohortMode.gcnv_call_tars
        Array[File] gcnv_tracking_tars = GermlineCNVCallerCohortMode.gcnv_tracking_tar
        Array[File] genotyped_intervals_vcfs = ordered_genotyped_intervals_vcf
        Array[File] genotyped_intervals_vcf_idxs = ordered_genotyped_intervals_vcf_idx
        Array[File] genotyped_segments_vcfs = ordered_genotyped_segments_vcf
        Array[File] genotyped_segments_vcf_idxs = ordered_genotyped_segments_vcf_idx
        Array[File] sample_qc_status_files = ordered_qc_status_file
        Array[String] sample_qc_status_strings = ordered_qc_status_string
        File model_qc_status_file = CollectModelQualityMetrics.qc_status_file
        String model_qc_string = CollectModelQualityMetrics.qc_status_string
        Array[File] denoised_copy_ratios = ordered_denoised_copy_ratios
    }
}

task AnnotateIntervals {
    input {
        File intervals
        String prefix
        File ref_fa
        File ref_fai
        File ref_dict
        File? mappability_track_bed
        File? mappability_track_bed_idx
        File? segmental_duplication_track_bed
        File? segmental_duplication_track_bed_idx
        Int? feature_query_lookahead
        File? gatk4_jar_override
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" AnnotateIntervals \
            -L ~{intervals} \
            --reference ~{ref_fa} \
            --sequence-dictionary ~{ref_dict} \
            ~{"--mappability-track " + mappability_track_bed} \
            ~{"--segmental-duplication-track " + segmental_duplication_track_bed} \
            --feature-query-lookahead ~{default=1000000 feature_query_lookahead} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{prefix}.annotated_intervals.tsv
    >>>

    output {
        File annotated_intervals = "~{prefix}.annotated_intervals.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size([intervals, ref_fa], "GB")) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task FilterIntervals {
    input {
        File intervals
        String prefix
        File annotated_intervals
        File? blacklist_intervals
        Array[File]? read_count_files
        Int? low_count_filter_count_threshold
        Float? low_count_filter_percentage_of_samples
        Float? extreme_count_filter_minimum_percentile
        Float? extreme_count_filter_maximum_percentile
        Float? extreme_count_filter_percentage_of_samples
        File? gatk4_jar_override
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        gatk --java-options "-Xmx~{command_mem_mb}m" FilterIntervals \
            -L ~{intervals} \
            ~{"-XL " + blacklist_intervals} \
            ~{if defined(read_count_files) then "--input " else ""} ~{sep=" --input " read_count_files} \
            ~{"--annotated-intervals " + annotated_intervals} \
            --low-count-filter-count-threshold ~{default="5" low_count_filter_count_threshold} \
            --low-count-filter-percentage-of-samples ~{default="90.0" low_count_filter_percentage_of_samples} \
            --extreme-count-filter-minimum-percentile ~{default="1.0" extreme_count_filter_minimum_percentile} \
            --extreme-count-filter-maximum-percentile ~{default="99.0" extreme_count_filter_maximum_percentile} \
            --extreme-count-filter-percentage-of-samples ~{default="90.0" extreme_count_filter_percentage_of_samples} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{prefix}.filtered_intervals.interval_list
    >>>

    output {
        File filtered_intervals = "~{prefix}.filtered_intervals.interval_list"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 7,
        disk_gb: ceil(size([intervals, annotated_intervals], "GB") * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ScatterIntervals {
    input {
        File interval_list
        String prefix
        Int num_intervals_per_scatter
        String? output_dir
        File? gatk4_jar_override
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    # Default the output directory to the task prefix
    String output_dir_ = select_first([output_dir, prefix + ".scattered_intervals"])

    command <<<
        set -euo pipefail

        # Create the output directory because IntervalListTools fails if it does not exist
        mkdir ~{output_dir_}
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        # Integer division gives the shard count because IntervalListTools puts remainder intervals in the last shard
        NUM_INTERVALS=$(grep -v '@' ~{interval_list} | wc -l)
        NUM_SCATTERS=$(echo $((NUM_INTERVALS / ~{num_intervals_per_scatter})))

        if [ $NUM_SCATTERS -le 1 ]; then
            # Copy the original interval list when only a single shard is required
            >&2 echo "Not running IntervalListTools because only a single shard is required. Copying original interval list..."
            cp ~{interval_list} ~{output_dir_}/~{prefix}.scattered.0001.interval_list
        else
            gatk --java-options "-Xmx~{command_mem_mb}m" IntervalListTools \
                --INPUT ~{interval_list} \
                --SUBDIVISION_MODE INTERVAL_COUNT \
                --SCATTER_CONTENT ~{num_intervals_per_scatter} \
                --OUTPUT ~{output_dir_}

            # Rename the per-shard interval lists from the IntervalListTools temp directory layout to the task prefix
            ls -v ~{output_dir_}/*/scattered.interval_list | \
                cat -n | \
                while read n filename; do mv $filename ~{output_dir_}/~{prefix}.scattered.$(printf "%04d" $n).interval_list; done
            rm -rf ~{output_dir_}/temp_*_of_*
        fi
    >>>

    output {
        Array[File] scattered_interval_lists = glob("~{output_dir_}/~{prefix}.scattered.*.interval_list")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size(interval_list, "GB") * 2) + 40,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task PostprocessGermlineCNVCalls {
    input {
        String prefix
        Array[File] gcnv_calls_tars
        Array[File] gcnv_model_tars
        Array[File] calling_configs
        Array[File] denoising_configs
        Array[File] gcnvkernel_version
        Array[File] sharded_interval_lists
        File contig_ploidy_calls_tar
        Array[String]? allosomal_contigs
        Int ref_copy_number_autosomal_contigs
        Int sample_index
        File? gatk4_jar_override
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    String genotyped_intervals_vcf_filename = "~{prefix}.genotyped_intervals.vcf.gz"
    String genotyped_segments_vcf_filename = "~{prefix}.genotyped_segments.vcf.gz"
    String denoised_copy_ratios_filename = "~{prefix}.denoised_copy_ratios.tsv"

    Array[String] allosomal_contigs_args = if defined(allosomal_contigs) then prefix("--allosomal-contig ", select_first([allosomal_contigs])) else []

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        sharded_interval_lists_array=(~{sep=" " sharded_interval_lists})

        # Untar calls into CALLS_0, CALLS_1, etc. directories with their shard config and interval files
        gcnv_calls_tar_array=(~{sep=" " gcnv_calls_tars})
        calling_configs_array=(~{sep=" " calling_configs})
        denoising_configs_array=(~{sep=" " denoising_configs})
        gcnvkernel_version_array=(~{sep=" " gcnvkernel_version})
        sharded_interval_lists_array=(~{sep=" " sharded_interval_lists})
        calls_args=""
        for index in ${!gcnv_calls_tar_array[@]}; do
            gcnv_calls_tar=${gcnv_calls_tar_array[$index]}
            mkdir -p CALLS_$index/SAMPLE_~{sample_index}
            tar xzf $gcnv_calls_tar -C CALLS_$index/SAMPLE_~{sample_index}
            cp ${calling_configs_array[$index]} CALLS_$index/
            cp ${denoising_configs_array[$index]} CALLS_$index/
            cp ${gcnvkernel_version_array[$index]} CALLS_$index/
            cp ${sharded_interval_lists_array[$index]} CALLS_$index/
            calls_args="$calls_args --calls-shard-path CALLS_$index"
        done

        # Untar models into MODEL_0, MODEL_1, etc. directories and build the command line
        gcnv_model_tar_array=(~{sep=" " gcnv_model_tars})
        model_args=""
        for index in ${!gcnv_model_tar_array[@]}; do
            gcnv_model_tar=${gcnv_model_tar_array[$index]}
            mkdir MODEL_$index
            tar xzf $gcnv_model_tar -C MODEL_$index
            model_args="$model_args --model-shard-path MODEL_$index"
        done

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        gatk --java-options "-Xmx~{command_mem_mb}m" PostprocessGermlineCNVCalls \
            $calls_args \
            $model_args \
            ~{sep=" " allosomal_contigs_args} \
            --autosomal-ref-copy-number ~{ref_copy_number_autosomal_contigs} \
            --contig-ploidy-calls contig-ploidy-calls \
            --sample-index ~{sample_index} \
            --output-genotyped-intervals ~{genotyped_intervals_vcf_filename} \
            --output-genotyped-segments ~{genotyped_segments_vcf_filename} \
            --output-denoised-copy-ratios ~{denoised_copy_ratios_filename}

        rm -rf CALLS_*
        rm -rf MODEL_*
        rm -rf contig-ploidy-calls
    >>>

    output {
        File genotyped_intervals_vcf = genotyped_intervals_vcf_filename
        File genotyped_intervals_vcf_idx = "~{genotyped_intervals_vcf_filename}.tbi"
        File genotyped_segments_vcf = genotyped_segments_vcf_filename
        File genotyped_segments_vcf_idx = "~{genotyped_segments_vcf_filename}.tbi"
        File denoised_copy_ratios = denoised_copy_ratios_filename
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 7,
        disk_gb: ceil((size(gcnv_calls_tars, "GB") + size(gcnv_model_tars, "GB") + size(calling_configs, "GB") + size(denoising_configs, "GB") + size(gcnvkernel_version, "GB") + size(sharded_interval_lists, "GB") + size(contig_ploidy_calls_tar, "GB")) * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CollectSampleQualityMetrics {
    input {
        File genotyped_segments_vcf
        String prefix
        Int maximum_number_events
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        NUM_SEGMENTS=$(gunzip -c ~{genotyped_segments_vcf} | awk '!/^#/ {count++} END {print count + 0}')
        if [ $NUM_SEGMENTS -lt ~{maximum_number_events} ]; then
            echo "PASS" >> ~{prefix}.qc_status.txt
        else
            echo "EXCESSIVE_NUMBER_OF_EVENTS" >> ~{prefix}.qc_status.txt
        fi
    >>>

    output {
        File qc_status_file = "~{prefix}.qc_status.txt"
        String qc_status_string = read_string("~{prefix}.qc_status.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: ceil(size(genotyped_segments_vcf, "GB")) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CollectModelQualityMetrics {
    input {
        Array[File] gcnv_model_tars
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        qc_status="PASS"

        gcnv_model_tar_array=(~{sep=" " gcnv_model_tars})
        for index in ${!gcnv_model_tar_array[@]}; do
            gcnv_model_tar=${gcnv_model_tar_array[$index]}
            mkdir MODEL_$index
            tar xzf $gcnv_model_tar -C MODEL_$index
            ard_file="MODEL_$index/mu_ard_u_interval__.tsv"

            # Check whether all ARD values are less than or equal to one
            NUM_ARD_VALUES_ABOVE_ONE=$(awk '!/^@/ && $1 != "VALUE_0" { ard = 1e10 / (1 + exp(-$1)); if (ard > 1.0) count++ } END { print count + 0 }' "$ard_file")
            if [ $NUM_ARD_VALUES_ABOVE_ONE -eq 0 ]; then
                qc_status="ALL_PRINCIPAL_COMPONENTS_USED"
                break
            fi
        done
        echo $qc_status >> ~{prefix}.model_qc_status.txt
    >>>

    output {
        File qc_status_file = "~{prefix}.model_qc_status.txt"
        String qc_status_string = read_string("~{prefix}.model_qc_status.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: ceil(size(gcnv_model_tars, "GB") * 2) + 40,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task DetermineGermlineContigPloidyCohortMode {
    input {
        String cohort_id
        String prefix
        File? intervals
        Array[File] read_count_files
        File contig_ploidy_priors
        String? output_dir
        File? gatk4_jar_override
        Float? mean_bias_standard_deviation
        Float? mapping_error_rate
        Float? global_psi_scale
        Float? sample_psi_scale
        String docker
        RuntimeAttr? runtime_attr_override
    }

    # Hybrid ADVI parameters are not exposed because the defaults are adequate

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    # Default the output directory to "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}
        export OMP_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}

        gatk --java-options "-Xmx~{command_mem_mb}m"  DetermineGermlineContigPloidy \
            ~{"-L " + intervals} \
            --input ~{sep=" --input " read_count_files} \
            --contig-ploidy-priors ~{contig_ploidy_priors} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{output_dir_} \
            --output-prefix ~{cohort_id} \
            --verbosity DEBUG \
            --mean-bias-standard-deviation ~{default="0.01" mean_bias_standard_deviation} \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --global-psi-scale ~{default="0.001" global_psi_scale} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale}

        tar czf ~{prefix}.contig_ploidy_model.tar.gz -C ~{output_dir_}/~{cohort_id}-model .
        tar czf ~{prefix}.contig_ploidy_calls.tar.gz -C ~{output_dir_}/~{cohort_id}-calls .
    >>>

    output {
        File contig_ploidy_model_tar = "~{prefix}.contig_ploidy_model.tar.gz"
        File contig_ploidy_calls_tar = "~{prefix}.contig_ploidy_calls.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 7,
        disk_gb: ceil(size(read_count_files, "GB") * 2 + size(contig_ploidy_priors, "GB")) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GermlineCNVCallerCohortMode {
    input {
        Int scatter_index
        String cohort_id
        String prefix
        Array[File] read_count_files
        File contig_ploidy_calls_tar
        File intervals
        File? annotated_intervals
        String? output_dir
        File? gatk4_jar_override
        Float? p_alt
        Float? p_active
        Float? cnv_coherence_length
        Float? class_coherence_length
        Int? max_copy_number
        Int? max_bias_factors
        Float? mapping_error_rate
        Float? interval_psi_scale
        Float? sample_psi_scale
        Float? depth_correction_tau
        Float? log_mean_bias_standard_deviation
        Float? init_ard_rel_unexplained_variance
        Int? num_gc_bins
        Float? gc_curve_standard_deviation
        String? copy_number_posterior_expectation_mode
        Boolean? enable_bias_factors
        Int? active_class_padding_hybrid_mode
        Float? learning_rate
        Float? adamax_beta_1
        Float? adamax_beta_2
        Int? log_emission_samples_per_round
        Float? log_emission_sampling_median_rel_error
        Int? log_emission_sampling_rounds
        Int? max_advi_iter_first_epoch
        Int? max_advi_iter_subsequent_epochs
        Int? min_training_epochs
        Int? max_training_epochs
        Float? initial_temperature
        Int? num_thermal_advi_iters
        Int? convergence_snr_averaging_window
        Float? convergence_snr_trigger_threshold
        Int? convergence_snr_countdown_window
        Int? max_calling_iters
        Float? caller_update_convergence_threshold
        Float? caller_internal_admixing_rate
        Float? caller_external_admixing_rate
        Boolean? disable_annealing
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    # Default the output directory to "out"
    String output_dir_ = select_first([output_dir, "out"])
    Int num_samples = length(read_count_files)

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}
        export OMP_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        gatk --java-options "-Xmx~{command_mem_mb}m"  GermlineCNVCaller \
            --run-mode COHORT \
            -L ~{intervals} \
            --input ~{sep=" --input " read_count_files} \
            --contig-ploidy-calls contig-ploidy-calls \
            ~{"--annotated-intervals " + annotated_intervals} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ~{output_dir_} \
            --output-prefix ~{cohort_id} \
            --verbosity DEBUG \
            --p-alt ~{default="1e-6" p_alt} \
            --p-active ~{default="1e-2" p_active} \
            --cnv-coherence-length ~{default="10000.0" cnv_coherence_length} \
            --class-coherence-length ~{default="10000.0" class_coherence_length} \
            --max-copy-number ~{default="5" max_copy_number} \
            --max-bias-factors ~{default="5" max_bias_factors} \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --interval-psi-scale ~{default="0.001" interval_psi_scale} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale} \
            --depth-correction-tau ~{default="10000.0" depth_correction_tau} \
            --log-mean-bias-standard-deviation ~{default="0.1" log_mean_bias_standard_deviation} \
            --init-ard-rel-unexplained-variance ~{default="0.1" init_ard_rel_unexplained_variance} \
            --num-gc-bins ~{default="20" num_gc_bins} \
            --gc-curve-standard-deviation ~{default="1.0" gc_curve_standard_deviation} \
            --copy-number-posterior-expectation-mode ~{default="HYBRID" copy_number_posterior_expectation_mode} \
            --enable-bias-factors ~{default="true" enable_bias_factors} \
            --active-class-padding-hybrid-mode ~{default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ~{default="0.05" learning_rate} \
            --adamax-beta-1 ~{default="0.9" adamax_beta_1} \
            --adamax-beta-2 ~{default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ~{default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ~{default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ~{default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ~{default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ~{default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ~{default="10" min_training_epochs} \
            --max-training-epochs ~{default="100" max_training_epochs} \
            --initial-temperature ~{default="2.0" initial_temperature} \
            --num-thermal-advi-iters ~{default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ~{default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ~{default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ~{default="10" convergence_snr_countdown_window} \
            --max-calling-iters ~{default="10" max_calling_iters} \
            --caller-update-convergence-threshold ~{default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ~{default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ~{default="1.00" caller_external_admixing_rate} \
            --disable-annealing ~{default="false" disable_annealing}

        tar czf ~{prefix}.gcnv_model.shard_~{scatter_index}.tar.gz -C ~{output_dir_}/~{cohort_id}-model .
        tar czf ~{prefix}.gcnv_tracking.shard_~{scatter_index}.tar.gz -C ~{output_dir_}/~{cohort_id}-tracking .

        CURRENT_SAMPLE=0
        NUM_SAMPLES=~{num_samples}
        NUM_DIGITS=${#NUM_SAMPLES}
        while [ $CURRENT_SAMPLE -lt $NUM_SAMPLES ]; do
            CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $CURRENT_SAMPLE)
            tar czf ~{prefix}.gcnv_calls.shard_~{scatter_index}.sample_$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz -C ~{output_dir_}/~{cohort_id}-calls/SAMPLE_$CURRENT_SAMPLE .
            CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
        done

        rm -rf contig-ploidy-calls
    >>>

    output {
        File gcnv_model_tar = "~{prefix}.gcnv_model.shard_~{scatter_index}.tar.gz"
        Array[File] gcnv_call_tars = glob("~{prefix}.gcnv_calls.shard_~{scatter_index}.sample_*.tar.gz")
        File gcnv_tracking_tar = "~{prefix}.gcnv_tracking.shard_~{scatter_index}.tar.gz"
        File calling_config_json = "~{output_dir_}/~{cohort_id}-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/~{cohort_id}-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/~{cohort_id}-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/~{cohort_id}-calls/interval_list.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 10,
        disk_gb: ceil((size(read_count_files, "GB") + size([contig_ploidy_calls_tar, intervals], "GB")) * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task DetermineGermlineContigPloidyCaseMode {
    input {
        String prefix
        Array[File] read_count_files
        File contig_ploidy_model_tar
        String? output_dir
        File? gatk4_jar_override
        Float? mapping_error_rate
        Float? sample_psi_scale
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    # Default the output directory to "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}
        export OMP_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}

        mkdir contig-ploidy-model
        tar xzf ~{contig_ploidy_model_tar} -C contig-ploidy-model

        # Case mode takes its intervals and ploidy priors from the fitted model, so neither may be passed again
        gatk --java-options "-Xmx~{command_mem_mb}m"  DetermineGermlineContigPloidy \
            --input ~{sep=" --input " read_count_files} \
            --model contig-ploidy-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale}

        tar czf ~{prefix}.case_contig_ploidy_calls.tar.gz -C ~{output_dir_}/case-calls .

        rm -rf contig-ploidy-model
    >>>

    output {
        File contig_ploidy_calls_tar = "~{prefix}.case_contig_ploidy_calls.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 7,
        disk_gb: ceil(size(read_count_files, "GB") * 2 + size(contig_ploidy_model_tar, "GB")) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GermlineCNVCallerCaseMode {
    input {
        Int scatter_index
        String prefix
        Array[String] sample_ids
        Array[File] read_count_files
        File contig_ploidy_calls_tar
        File gcnv_model_tar
        String? output_dir
        File? gatk4_jar_override
        Float? p_alt
        Float? cnv_coherence_length
        Int? max_copy_number
        Float? mapping_error_rate
        Float? sample_psi_scale
        Float? depth_correction_tau
        String? copy_number_posterior_expectation_mode
        Int? active_class_padding_hybrid_mode
        Float? learning_rate
        Float? adamax_beta_1
        Float? adamax_beta_2
        Int? log_emission_samples_per_round
        Float? log_emission_sampling_median_rel_error
        Int? log_emission_sampling_rounds
        Int? max_advi_iter_first_epoch
        Int? max_advi_iter_subsequent_epochs
        Int? min_training_epochs
        Int? max_training_epochs
        Float? initial_temperature
        Int? num_thermal_advi_iters
        Int? convergence_snr_averaging_window
        Float? convergence_snr_trigger_threshold
        Int? convergence_snr_countdown_window
        Int? max_calling_iters
        Float? caller_update_convergence_threshold
        Float? caller_internal_admixing_rate
        Float? caller_external_admixing_rate
        Boolean? disable_annealing
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int command_mem_mb = ceil(select_first([runtime_attr.mem_gb, default_attr.mem_gb]) * 0.8 * 1024)

    # Default the output directory to "out"
    String output_dir_ = select_first([output_dir, "out"])
    Int num_samples = length(read_count_files)

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}
        export OMP_NUM_THREADS=~{select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])}

        mkdir contig-ploidy-calls
        tar xzf ~{contig_ploidy_calls_tar} -C contig-ploidy-calls

        mkdir gcnv-model
        tar xzf ~{gcnv_model_tar} -C gcnv-model

        # Case mode takes its intervals and denoising hyperparameters from the fitted model, so neither may be passed again
        gatk --java-options "-Xmx~{command_mem_mb}m"  GermlineCNVCaller \
            --run-mode CASE \
            --input ~{sep=" --input " read_count_files} \
            --contig-ploidy-calls contig-ploidy-calls \
            --model gcnv-model \
            --output ~{output_dir_} \
            --output-prefix case \
            --verbosity DEBUG \
            --p-alt ~{default="1e-6" p_alt} \
            --cnv-coherence-length ~{default="10000.0" cnv_coherence_length} \
            --max-copy-number ~{default="5" max_copy_number} \
            --mapping-error-rate ~{default="0.01" mapping_error_rate} \
            --sample-psi-scale ~{default="0.0001" sample_psi_scale} \
            --depth-correction-tau ~{default="10000.0" depth_correction_tau} \
            --copy-number-posterior-expectation-mode ~{default="HYBRID" copy_number_posterior_expectation_mode} \
            --active-class-padding-hybrid-mode ~{default="50000" active_class_padding_hybrid_mode} \
            --learning-rate ~{default="0.05" learning_rate} \
            --adamax-beta-1 ~{default="0.9" adamax_beta_1} \
            --adamax-beta-2 ~{default="0.99" adamax_beta_2} \
            --log-emission-samples-per-round ~{default="50" log_emission_samples_per_round} \
            --log-emission-sampling-median-rel-error ~{default="0.005" log_emission_sampling_median_rel_error} \
            --log-emission-sampling-rounds ~{default="10" log_emission_sampling_rounds} \
            --max-advi-iter-first-epoch ~{default="5000" max_advi_iter_first_epoch} \
            --max-advi-iter-subsequent-epochs ~{default="100" max_advi_iter_subsequent_epochs} \
            --min-training-epochs ~{default="10" min_training_epochs} \
            --max-training-epochs ~{default="100" max_training_epochs} \
            --initial-temperature ~{default="2.0" initial_temperature} \
            --num-thermal-advi-iters ~{default="2500" num_thermal_advi_iters} \
            --convergence-snr-averaging-window ~{default="500" convergence_snr_averaging_window} \
            --convergence-snr-trigger-threshold ~{default="0.1" convergence_snr_trigger_threshold} \
            --convergence-snr-countdown-window ~{default="10" convergence_snr_countdown_window} \
            --max-calling-iters ~{default="10" max_calling_iters} \
            --caller-update-convergence-threshold ~{default="0.001" caller_update_convergence_threshold} \
            --caller-internal-admixing-rate ~{default="0.75" caller_internal_admixing_rate} \
            --caller-external-admixing-rate ~{default="1.00" caller_external_admixing_rate} \
            --disable-annealing ~{default="false" disable_annealing}

        tar czf ~{prefix}.case_gcnv_tracking.shard_~{scatter_index}.tar.gz -C ~{output_dir_}/case-tracking .

        # Fail loudly if the call directories are not in the order the read counts were given, which downstream indexing assumes
        expected_sample_ids=(~{sep=" " sample_ids})
        for index in ${!expected_sample_ids[@]}; do
            actual_sample_id="$(cat ~{output_dir_}/case-calls/SAMPLE_$index/sample_name.txt)"
            if [[ "${expected_sample_ids[$index]}" != "${actual_sample_id}" ]]; then
                printf 'Expected sample ID does not match actual sample ID for SAMPLE_%s\n' "$index" >&2
                printf 'Expected: %s\n' "${expected_sample_ids[$index]}" >&2
                printf 'Actual: %s\n' "${actual_sample_id}" >&2
                exit 1
            fi
        done

        CURRENT_SAMPLE=0
        NUM_SAMPLES=~{num_samples}
        NUM_DIGITS=${#NUM_SAMPLES}
        while [ $CURRENT_SAMPLE -lt $NUM_SAMPLES ]; do
            CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $CURRENT_SAMPLE)
            tar czf ~{prefix}.case_gcnv_calls.shard_~{scatter_index}.sample_$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz -C ~{output_dir_}/case-calls/SAMPLE_$CURRENT_SAMPLE .
            CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
        done

        rm -rf contig-ploidy-calls
        rm -rf gcnv-model
    >>>

    output {
        Array[File] gcnv_call_tars = glob("~{prefix}.case_gcnv_calls.shard_~{scatter_index}.sample_*.tar.gz")
        File gcnv_tracking_tar = "~{prefix}.case_gcnv_tracking.shard_~{scatter_index}.tar.gz"
        File calling_config_json = "~{output_dir_}/case-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/case-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/case-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/case-calls/interval_list.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 8,
        mem_gb: 10,
        disk_gb: ceil((size(read_count_files, "GB") + size([contig_ploidy_calls_tar, gcnv_model_tar], "GB")) * 2) + 50,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task MergeContigPloidyCalls {
    input {
        String prefix
        File training_contig_ploidy_calls_tar
        File case_contig_ploidy_calls_tar
        Array[Int] training_indices
        Array[Int] case_indices
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir training case merged
        tar xzf ~{training_contig_ploidy_calls_tar} -C training
        tar xzf ~{case_contig_ploidy_calls_tar} -C case

        # Renumber each call directory from its position within its own run to its position in the full sample list
        training_indices=(~{sep=" " training_indices})
        for index in ${!training_indices[@]}; do
            mv training/SAMPLE_$index merged/SAMPLE_${training_indices[$index]}
        done
        case_indices=(~{sep=" " case_indices})
        for index in ${!case_indices[@]}; do
            mv case/SAMPLE_$index merged/SAMPLE_${case_indices[$index]}
        done

        # Carry over the run-level files that sit alongside the call directories
        find training -maxdepth 1 -type f -exec cp {} merged/ \;

        tar czf ~{prefix}.contig_ploidy_calls.tar.gz -C merged .

        rm -rf training case merged
    >>>

    output {
        File contig_ploidy_calls_tar = "~{prefix}.contig_ploidy_calls.tar.gz"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 2,
        disk_gb: ceil(size([training_contig_ploidy_calls_tar, case_contig_ploidy_calls_tar], "GB") * 4) + 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
