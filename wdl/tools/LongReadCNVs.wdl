version 1.0

import "../utils/LRCNVs.wdl"
import "../utils/DepthPreprocessing.wdl"
import "../utils/DepthClustering.wdl"
import "../utils/GenotypeDepth.wdl"
import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow LongReadCNVs {
    meta {
        description: [
            "This workflow calls cohort CNVs from long-read depth profiles with GATK gCNV, then converts, clusters and genotypes the depth calls. It outputs merged CNV calls, ploidy, and genotyped depth VCFs.",
            "By default every sample is called in gCNV cohort mode. Setting `num_training_samples` below the cohort size instead runs a hybrid case-cohort mode, fitting the contig-ploidy and gCNV models on that many randomly drawn samples and calling every remaining sample against those models in case mode. All downstream steps and outputs still cover the whole cohort either way."
        ]
    }

    parameter_meta {
        intervals: "Interval list over which CNVs are called."
        sample_ids: "Sample IDs in the cohort."
        depth_profiles: "Per-sample read-depth profiles, aligned to `sample_ids`."
        sort_depth_profiles: "Whether to sort each depth profile by contig and position first, for profiles that are not already coordinate sorted."
        batch_id: "Identifier for the cohort batch."
        contig_ploidy_priors: "Contig ploidy priors used to determine per-sample contig ploidy."
        merged_bincov: "Merged read-depth evidence and its tabix index for depth genotyping."
        merged_bincov_idx: "Index for merged_bincov."
        ref_fa: "From references."
        ref_fai: "From references."
        ref_dict: "From references."
        pedigree: "Cohort pedigree used by depth genotyping."
        primary_contigs_list: "Primary contigs, one per line in reference dictionary order, used for the ploidy table, VCF headers, clustering and genotyping."
        training_intervals: "Intervals used to train the depth genotyping model."
        median_coverage: "Per-sample median coverage table used by depth genotyping."
        contig_subset_list: "Optional subset of `primary_contigs_list` to restrict depth clustering and genotyping to."
        variant_prefix: "Prefix used for generated variant IDs."
        gcnv_qs_cutoff: "Minimum gCNV quality score for a segment to be kept."
        num_intervals_per_scatter: "Number of intervals processed per gCNV scatter shard. GermlineCNVCaller memory grows with samples times intervals per shard, so raising this above the default needs more memory in `runtime_attr_germline_cnv_caller`."
        num_training_samples: "Number of samples drawn at random to fit the contig-ploidy and gCNV models in cohort mode, with every remaining sample called against those models in case mode. Left unset, or set to at least the cohort size, every sample is called in cohort mode. Interval filtering percentages then apply over the training samples only, so a training set of fewer than a few dozen samples degrades the fitted models."
        subsample_seed: "Random seed used to draw the training samples."
        chr_x: "Name of the X contig in the reference."
        chr_y: "Name of the Y contig in the reference."
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
        defragment_max_dist: "Maximum gap, as a fraction of call length, across which adjacent calls are defragmented."
        fast_mode: "Use SVCluster fast mode."
        clustering_algorithm: "SVCluster algorithm."
        enable_cnv: "SVCluster behavior flags."
        default_no_call: "SVCluster behavior flags."
        omit_members: "SVCluster behavior flags."
        breakpoint_summary_strategy: "SVCluster behavior flags."
        defrag_padding_fraction: "Defragmentation thresholds."
        defrag_sample_overlap: "Defragmentation thresholds."
        depth_sample_overlap: "Required sample overlap for depth clustering."
        depth_interval_overlap: "Required reciprocal interval overlap."
        depth_size_similarity: "Required size similarity."
        depth_breakend_window: "Breakend join window in base pairs."
        exclude_intervals: "Intervals whose overlapping calls are dropped."
        exclude_overlap_fraction: "Overlap fraction at which a call is excluded."
        gatk_to_svtk_script: "Override for the GATK-to-svtk conversion script."
        svtk_set_pass: "Set FILTER to PASS during conversion."
        merged_cnvs_vcf: "Cohort CNV VCF after depth preprocessing."
        merged_cnvs_vcf_idx: "Index for `merged_cnvs_vcf`."
        ploidy_table: "Per-sample ploidy table."
        genotyped_depth_vcf: "Clustered CNV VCF genotyped from read depth."
        genotyped_depth_vcf_idx: "Index for `genotyped_depth_vcf`."
        genotyping_rd_table: "Read-depth evidence used for genotyping."
    }

    input {
        File intervals
        Array[String]+ sample_ids
        Array[File]+ depth_profiles
        Boolean sort_depth_profiles
        String batch_id
        File contig_ploidy_priors
        File merged_bincov
        File merged_bincov_idx
        File ref_fa
        File ref_fai
        File ref_dict

        File pedigree
        File primary_contigs_list
        File training_intervals
        File median_coverage
        File? contig_subset_list

        String prefix
        String variant_prefix

        Int gcnv_qs_cutoff = 30
        Int num_intervals_per_scatter = 1500
        Int? num_training_samples
        Int subsample_seed = 42
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
            num_training_samples = num_training_samples,
            subsample_seed = subsample_seed,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            gatk_docker = gatk_docker,
            sv_pipeline_docker = sv_pipeline_docker,
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
            runtime_attr_subsample_indices = runtime_attr_subsample_indices,
            runtime_attr_annotate_intervals = runtime_attr_annotate_intervals,
            runtime_attr_filter_intervals = runtime_attr_filter_intervals,
            runtime_attr_scatter_intervals = runtime_attr_scatter_intervals,
            runtime_attr_determine_contig_ploidy = runtime_attr_determine_contig_ploidy,
            runtime_attr_determine_contig_ploidy_case = runtime_attr_determine_contig_ploidy_case,
            runtime_attr_germline_cnv_caller = runtime_attr_germline_cnv_caller,
            runtime_attr_germline_cnv_caller_case = runtime_attr_germline_cnv_caller_case,
            runtime_attr_postprocess_germline_cnv_calls = runtime_attr_postprocess_germline_cnv_calls,
            runtime_attr_collect_sample_quality_metrics = runtime_attr_collect_sample_quality_metrics,
            runtime_attr_collect_model_quality_metrics = runtime_attr_collect_model_quality_metrics,
            runtime_attr_merge_contig_ploidy_calls = runtime_attr_merge_contig_ploidy_calls
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
