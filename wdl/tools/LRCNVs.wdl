version 1.0

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

import "../utils/Helpers.wdl"
import "../utils/Structs.wdl"

workflow LRCNVs {
    meta {
        description: "Workflow for creating a GATK GermlineCNVCaller denoising model and generating calls given a list of samples with long-read sequencing reads."
    }

    parameter_meta {
        intervals: "GATK-style intervals used to collect depth profiles."
        sample_ids: "Identifier for each sample."
        depth_profiles: "Median depth at each interval, rounded to an integer. One file per sample. See the TSV format here https://gatk.broadinstitute.org/hc/en-us/articles/35967568802843-CollectReadCounts."
        cohort_id: "Identifier for the cohort used for denoising model generation."
        contig_ploidy_priors: "File containing contig ploidy priors. See input comment."
        num_intervals_per_scatter: "Number of intervals to process in each scatter."
        ref_fa: "Reference sequences FASTA file."
        ref_fai: "Reference sequences FASTA file index."
        ref_dict: "Reference sequences dictionary."
        gatk_docker: "Docker image for the GATK tool."
    }

    input {
        File intervals
        Array[String]+ sample_ids
        Array[String]+ depth_profiles
        String cohort_id
        # A TSV with the prior probability of each ploidy state of each contig.
        # e.g.
        # CONTIG_NAME PLOIDY_PRIOR_0 PLOIDY_PRIOR_1 PLOIDY_PRIOR_2 PLOIDY_PRIOR_3
        #        chr1            0.0           0.01           0.98           0.01
        #        chr2            0.0           0.01           0.98           0.01
        # ...
        #        chrX           0.01           0.49           0.49           0.01
        #        chrY          0.495          0.495           0.01            0.0
        File contig_ploidy_priors
        Int num_intervals_per_scatter = 10000
        File ref_fa
        File ref_fai
        File ref_dict
        String gatk_docker = "broadinstitute/gatk:4.6.2.0"

        File? gatk4_jar_override 
        Int? preemptible_attempts

        # AnnotateIntervals
        File? mappability_track_bed
        File? mappability_track_bed_idx
        File? segmental_duplication_track_bed
        File? segmental_duplication_track_bed_idx
        Int? feature_query_lookahead
        Int? mem_gb_for_annotate_intervals

        # FilterIntervals
        File? blacklist_intervals
        Int? low_count_filter_count_threshold
        Float? low_count_filter_percentage_of_samples
        Float? extreme_count_filter_minimum_percentile
        Float? extreme_count_filter_maximum_percentile
        Float? extreme_count_filter_percentage_of_samples
        Int? mem_gb_for_filter_intervals

        # DeterminGermlineContigPloidyCohortMode
        Float? ploidy_mean_bias_standard_deviation
        Float? ploidy_mapping_error_rate
        Float? ploidy_global_psi_scale
        Float? ploidy_sample_psi_scale
        Int? mem_gb_for_determine_germline_contig_ploidy
        Int? cpu_for_determine_germline_contig_ploidy

        # GermlineCNVCallerCohortMode
        Float? gcnv_p_alt
        Float? gcnv_p_active
        Float? gcnv_cnv_coherence_length
        Float? gcnv_class_coherence_length
        Int? gcnv_max_copy_number
        Int? mem_gb_for_germline_cnv_caller
        Int? cpu_for_germline_cnv_caller

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
    }

    call AnnotateIntervals {
        input:
            intervals = intervals,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            mappability_track_bed = mappability_track_bed,
            mappability_track_bed_idx = mappability_track_bed_idx,
            segmental_duplication_track_bed = segmental_duplication_track_bed,
            segmental_duplication_track_bed_idx = segmental_duplication_track_bed_idx,
            feature_query_lookahead = feature_query_lookahead,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_annotate_intervals,
            preemptible_attempts = preemptible_attempts
    }

    call FilterIntervals {
        input:
            intervals = intervals,
            annotated_intervals = AnnotateIntervals.annotated_intervals,
            blacklist_intervals = blacklist_intervals,
            read_count_files = depth_profiles,
            low_count_filter_count_threshold = low_count_filter_count_threshold,
            low_count_filter_percentage_of_samples = low_count_filter_percentage_of_samples,
            extreme_count_filter_minimum_percentile = extreme_count_filter_minimum_percentile,
            extreme_count_filter_maximum_percentile = extreme_count_filter_maximum_percentile,
            extreme_count_filter_percentage_of_samples = extreme_count_filter_percentage_of_samples,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_filter_intervals,
            preemptible_attempts = preemptible_attempts
    }

    call DetermineGermlineContigPloidyCohortMode {
        input:
            cohort_id = cohort_id,
            intervals = FilterIntervals.filtered_intervals,
            read_count_files = depth_profiles,
            contig_ploidy_priors = contig_ploidy_priors,
            gatk4_jar_override = gatk4_jar_override,
            gatk_docker = gatk_docker,
            mem_gb = mem_gb_for_determine_germline_contig_ploidy,
            cpu = cpu_for_determine_germline_contig_ploidy,
            mean_bias_standard_deviation = ploidy_mean_bias_standard_deviation,
            mapping_error_rate = ploidy_mapping_error_rate,
            global_psi_scale = ploidy_global_psi_scale,
            sample_psi_scale = ploidy_sample_psi_scale,
            preemptible_attempts = preemptible_attempts
    }

    call ScatterIntervals {
        input:
            interval_list = FilterIntervals.filtered_intervals,
            num_intervals_per_scatter = num_intervals_per_scatter,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    scatter (scatter_index in range(length(ScatterIntervals.scattered_interval_lists))) {
        call GermlineCNVCallerCohortMode {
            input:
                scatter_index = scatter_index,
                cohort_id = cohort_id,
                read_count_files = depth_profiles,
                contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar,
                intervals = ScatterIntervals.scattered_interval_lists[scatter_index],
                annotated_intervals = AnnotateIntervals.annotated_intervals,
                gatk4_jar_override = gatk4_jar_override,
                gatk_docker = gatk_docker,
                mem_gb = mem_gb_for_germline_cnv_caller,
                cpu = cpu_for_germline_cnv_caller,
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
                preemptible_attempts = preemptible_attempts
        }
    }

    Array[Array[File]] call_tars_sample_by_shard = transpose(GermlineCNVCallerCohortMode.gcnv_call_tars)

    scatter (sample_index in range(length(sample_ids))) {
        call PostprocessGermlineCNVCalls {
            input:
                sample_id = sample_ids[sample_index],
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
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts
        }

        call CollectSampleQualityMetrics {
            input:
                genotyped_segments_vcf = PostprocessGermlineCNVCalls.genotyped_segments_vcf,
                sample_id = sample_ids[sample_index],
                maximum_number_events = maximum_number_events_per_sample,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts
        }
    }

    call CollectModelQualityMetrics {
        input:
            gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar,
            gatk_docker = gatk_docker,
            preemptible_attempts = preemptible_attempts
    }

    output {
        File annotated_intervals = AnnotateIntervals.annotated_intervals
        File filtered_intervals = FilterIntervals.filtered_intervals
        File contig_ploidy_model_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_model_tar
        File contig_ploidy_calls_tar = DetermineGermlineContigPloidyCohortMode.contig_ploidy_calls_tar
        Array[File] gcnv_model_tars = GermlineCNVCallerCohortMode.gcnv_model_tar
        Array[Array[File]] gcnv_calls_tars = GermlineCNVCallerCohortMode.gcnv_call_tars
        Array[File] gcnv_tracking_tars = GermlineCNVCallerCohortMode.gcnv_tracking_tar
        Array[File] genotyped_intervals_vcfs = PostprocessGermlineCNVCalls.genotyped_intervals_vcf
        Array[File] genotyped_segments_vcfs = PostprocessGermlineCNVCalls.genotyped_segments_vcf
        Array[File] sample_qc_status_files = CollectSampleQualityMetrics.qc_status_file
        Array[String] sample_qc_status_strings = CollectSampleQualityMetrics.qc_status_string
        File model_qc_status_file = CollectModelQualityMetrics.qc_status_file
        String model_qc_string = CollectModelQualityMetrics.qc_status_string
        Array[File] denoised_copy_ratios = PostprocessGermlineCNVCalls.denoised_copy_ratios
    }
}

task AnnotateIntervals {
    input {
        File intervals
        File ref_fa
        File ref_fai
        File ref_dict
        File? mappability_track_bed
        File? mappability_track_bed_idx
        File? segmental_duplication_track_bed
        File? segmental_duplication_track_bed_idx
        Int? feature_query_lookahead
        File? gatk4_jar_override

        # Runtime parameters
        String gatk_docker
        Int? mem_gb
        Int? disk_space_gb
        Boolean use_ssd = false
        Int? cpu
        Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = ceil(machine_mem_mb * 0.8)
    
    # Determine output filename
    String base_filename = basename(intervals, ".interval_list")

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
            --output ~{base_filename}.annotated.tsv
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, ceil(size(ref_fa, "GB")) + 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File annotated_intervals = "~{base_filename}.annotated.tsv"
    }
}

task FilterIntervals {
    input {
      File intervals
      File annotated_intervals
      File? blacklist_intervals
      Array[File]? read_count_files
      Int? low_count_filter_count_threshold
      Float? low_count_filter_percentage_of_samples
      Float? extreme_count_filter_minimum_percentile
      Float? extreme_count_filter_maximum_percentile
      Float? extreme_count_filter_percentage_of_samples
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = ceil(machine_mem_mb * 0.8)
    
    # Determine output filename
    String base_filename = basename(intervals, ".intervals")

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
            --output ~{base_filename}.filtered.interval_list
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 50]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File filtered_intervals = "~{base_filename}.filtered.interval_list"
    }
}

task ScatterIntervals {
    input {
      File interval_list
      Int num_intervals_per_scatter
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 2]) * 1000
    Int command_mem_mb = ceil(machine_mem_mb * 0.8)

    # If optional output_dir not specified, use "out";
    String output_dir_ = select_first([output_dir, "out"])

    String base_filename = basename(interval_list, ".interval_list")

    command <<<
        set -euo pipefail

        # IntervalListTools will fail if the output directory does not exist, so we create it
        mkdir ~{output_dir_}
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        # IntervalListTools behaves differently when scattering to a single or multiple shards, so we do some handling in bash

        # IntervalListTools tries to equally divide intervals across shards to give at least INTERVAL_COUNT in each and
        # puts remainder intervals in the last shard, so integer division gives the number of shards
        # (unless NUM_INTERVALS < num_intervals_per_scatter and NUM_SCATTERS = 0, in which case we still want a single shard)
        NUM_INTERVALS=$(grep -v '@' ~{interval_list} | wc -l)
        NUM_SCATTERS=$(echo $((NUM_INTERVALS / ~{num_intervals_per_scatter})))

        if [ $NUM_SCATTERS -le 1 ]; then
            # if only a single shard is required, then we can just rename the original interval list
            >&2 echo "Not running IntervalListTools because only a single shard is required. Copying original interval list..."
            cp ~{interval_list} ~{output_dir_}/~{base_filename}.scattered.0001.interval_list
        else
            gatk --java-options "-Xmx~{command_mem_mb}m" IntervalListTools \
                --INPUT ~{interval_list} \
                --SUBDIVISION_MODE INTERVAL_COUNT \
                --SCATTER_CONTENT ~{num_intervals_per_scatter} \
                --OUTPUT ~{output_dir_}

            # output files are named output_dir_/temp_0001_of_N/scattered.interval_list, etc. (N = number of scatters);
            # we rename them as output_dir_/base_filename.scattered.0001.interval_list, etc.
            ls -v ~{output_dir_}/*/scattered.interval_list | \
                cat -n | \
                while read n filename; do mv $filename ~{output_dir_}/~{base_filename}.scattered.$(printf "%04d" $n).interval_list; done
            rm -rf ~{output_dir_}/temp_*_of_*
        fi
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        Array[File] scattered_interval_lists = glob("~{output_dir_}/~{base_filename}.scattered.*.interval_list")
    }
}

task PostprocessGermlineCNVCalls {
    input {
      String sample_id
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

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = ceil(machine_mem_mb * 0.8)

    String genotyped_intervals_vcf_filename = "genotyped-intervals-~{sample_id}.vcf.gz"
    String genotyped_segments_vcf_filename = "genotyped-segments-~{sample_id}.vcf.gz"
    String denoised_copy_ratios_filename = "denoised_copy_ratios-~{sample_id}.tsv"

    Array[String] allosomal_contigs_args = if defined(allosomal_contigs) then prefix("--allosomal-contig ", select_first([allosomal_contigs])) else []

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}

        sharded_interval_lists_array=(~{sep=" " sharded_interval_lists})

        # untar calls to CALLS_0, CALLS_1, etc directories and build the command line
        # also copy over shard config and interval files
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

        # untar models to MODEL_0, MODEL_1, etc directories and build the command line
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

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File genotyped_intervals_vcf = genotyped_intervals_vcf_filename
        File genotyped_segments_vcf = genotyped_segments_vcf_filename
        File denoised_copy_ratios = denoised_copy_ratios_filename
    }
}

task CollectSampleQualityMetrics {
    input {
      File genotyped_segments_vcf
      String sample_id
      Int maximum_number_events

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 1]) * 1000

    command <<<
        set -euo pipefail

        NUM_SEGMENTS=$(gunzip -c ~{genotyped_segments_vcf} | grep -v '#' | wc -l)
        if [ $NUM_SEGMENTS -lt ~{maximum_number_events} ]; then
            echo "PASS" >> ~{sample_id}.qcStatus.txt
        else 
            echo "EXCESSIVE_NUMBER_OF_EVENTS" >> ~{sample_id}.qcStatus.txt
        fi
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 20]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File qc_status_file = "~{sample_id}.qcStatus.txt"
        String qc_status_string = read_string("~{sample_id}.qcStatus.txt")
    }
}

task CollectModelQualityMetrics {
    input {
      Array[File] gcnv_model_tars

     # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts
    }

    Int machine_mem_mb = select_first([mem_gb, 1]) * 1000

    command <<<
        sed -e 
        qc_status="PASS"

        gcnv_model_tar_array=(~{sep=" " gcnv_model_tars})
        for index in ${!gcnv_model_tar_array[@]}; do
            gcnv_model_tar=${gcnv_model_tar_array[$index]}
            mkdir MODEL_$index
            tar xzf $gcnv_model_tar -C MODEL_$index
            ard_file=MODEL_$index/mu_ard_u_log__.tsv

            #check whether all values for ARD components are negative
            NUM_POSITIVE_VALUES=$(awk '{ if (index($0, "@") == 0) {if ($1 > 0.0) {print $1} }}' MODEL_$index/mu_ard_u_log__.tsv | wc -l)
            if [ $NUM_POSITIVE_VALUES -eq 0 ]; then
                qc_status="ALL_PRINCIPAL_COMPONENTS_USED"
                break
            fi
        done
        echo $qc_status >> qcStatus.txt
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 40]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 1])
        preemptible: select_first([preemptible_attempts, 5])
    }

    output {
        File qc_status_file = "qcStatus.txt"
        String qc_status_string = read_string("qcStatus.txt")
    }
}

task DetermineGermlineContigPloidyCohortMode {
    input {
      String cohort_id
      File? intervals
      Array[File] read_count_files
      File contig_ploidy_priors
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Model parameters
      Float? mean_bias_standard_deviation
      Float? mapping_error_rate
      Float? global_psi_scale
      Float? sample_psi_scale
    }

    # We do not expose Hybrid ADVI parameters -- the default values are decent

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = ceil(machine_mem_mb * 0.8)

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])

    command <<<
        set -xeuo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

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

        tar czf ~{cohort_id}-contig-ploidy-model.tar.gz -C ~{output_dir_}/~{cohort_id}-model .
        tar czf ~{cohort_id}-contig-ploidy-calls.tar.gz -C ~{output_dir_}/~{cohort_id}-calls .
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File contig_ploidy_model_tar = "~{cohort_id}-contig-ploidy-model.tar.gz"
        File contig_ploidy_calls_tar = "~{cohort_id}-contig-ploidy-calls.tar.gz"
    }
}

task GermlineCNVCallerCohortMode {
    input {
      Int scatter_index
      String cohort_id
      Array[File] read_count_files
      File contig_ploidy_calls_tar
      File intervals
      File? annotated_intervals
      String? output_dir
      File? gatk4_jar_override

      # Runtime parameters
      String gatk_docker
      Int? mem_gb
      Int? disk_space_gb
      Boolean use_ssd = false
      Int? cpu
      Int? preemptible_attempts

      # Caller parameters
      Float? p_alt
      Float? p_active
      Float? cnv_coherence_length
      Float? class_coherence_length
      Int? max_copy_number

      # Denoising model parameters
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

      # Hybrid ADVI parameters
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
    }

    Int machine_mem_mb = select_first([mem_gb, 7]) * 1000
    Int command_mem_mb = ceil(machine_mem_mb * 0.8)

    # If optional output_dir not specified, use "out"
    String output_dir_ = select_first([output_dir, "out"])
    Int num_samples = length(read_count_files)

    String dollar = "$" #WDL workaround, see https://github.com/broadinstitute/cromwell/issues/1819

    command <<<
        set -euo pipefail

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk4_jar_override}
        export MKL_NUM_THREADS=~{default=8 cpu}
        export OMP_NUM_THREADS=~{default=8 cpu}

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

        tar czf ~{cohort_id}-gcnv-model-shard-~{scatter_index}.tar.gz -C ~{output_dir_}/~{cohort_id}-model .
        tar czf ~{cohort_id}-gcnv-tracking-shard-~{scatter_index}.tar.gz -C ~{output_dir_}/~{cohort_id}-tracking .

        CURRENT_SAMPLE=0
        NUM_SAMPLES=~{num_samples}
        NUM_DIGITS=${#NUM_SAMPLES}
        while [ $CURRENT_SAMPLE -lt $NUM_SAMPLES ]; do
            CURRENT_SAMPLE_WITH_LEADING_ZEROS=$(printf "%0${NUM_DIGITS}d" $CURRENT_SAMPLE)
            tar czf ~{cohort_id}-gcnv-calls-shard-~{scatter_index}-sample-$CURRENT_SAMPLE_WITH_LEADING_ZEROS.tar.gz -C ~{output_dir_}/~{cohort_id}-calls/SAMPLE_$CURRENT_SAMPLE .
            CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
        done

        rm -rf contig-ploidy-calls
    >>>

    runtime {
        docker: gatk_docker
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        File gcnv_model_tar = "~{cohort_id}-gcnv-model-shard-~{scatter_index}.tar.gz"
        Array[File] gcnv_call_tars = glob("~{cohort_id}-gcnv-calls-shard-~{scatter_index}-sample-*.tar.gz")
        File gcnv_tracking_tar = "~{cohort_id}-gcnv-tracking-shard-~{scatter_index}.tar.gz"
        File calling_config_json = "~{output_dir_}/~{cohort_id}-calls/calling_config.json"
        File denoising_config_json = "~{output_dir_}/~{cohort_id}-calls/denoising_config.json"
        File gcnvkernel_version_json = "~{output_dir_}/~{cohort_id}-calls/gcnvkernel_version.json"
        File sharded_interval_list = "~{output_dir_}/~{cohort_id}-calls/interval_list.tsv"
    }
}
