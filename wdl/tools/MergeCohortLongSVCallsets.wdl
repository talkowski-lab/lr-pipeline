# Derived from yuliamostovoy/lr_callset_integration, a fork of fabio-cunial/callset_integration_phase2:
# https://github.com/yuliamostovoy/lr_callset_integration/blob/main/wdl/SV_Integration_WorkflowD_Ultralong_Bnd.wdl

version 1.0

import "../utils/MergeCohortLongSVClass.wdl"
import "../utils/Structs.wdl"

workflow MergeCohortLongSVCallsets {
    meta {
        description: [
            "This tool merges the per-sample SVs longer than the MergeSampleSVCallsets length range, and the per-sample breakends, into one cohort callset for each class. Each class runs through the same bcftools merge, truvari-collapse sharding, Truvari (https://github.com/ACEnglish/truvari) collapse and concatenation stages as MergeCohortSVCallsets, without scoring or regenotyping.",
            "Nothing is returned as a workflow output. Each class writes under its own subdirectory of 'remote_outdir', named after the class, and that class's cohort callset is '15_concat/truvari_collapsed.bcf' there."
        ]
    }

    parameter_meta {
        remote_indir: "The 'remote_outdir' of MergeSampleSVCallsets."
        remote_outdir: "GCS directory each class's stages are written under, without a trailing slash."
        suffixes: "Classes to merge: 'ultralong' for the longer SVs and 'bnd' for the breakends."
        intrasample_subdir: "Subdirectory of 'remote_indir' holding the per-sample calls of each class. Leave empty to read them from 'remote_indir' itself."
        chromosomes: "Chromosomes to process, in output order."
        n_expected_samples: "Number of samples to merge. Derived from the per-sample files of each class when omitted."
        truvari_chunk_min_records: "Minimum number of records in each truvari-collapse shard."
        truvari_collapse_refdist: "Distance, in bp, that shard boundaries keep from any record, so that records Truvari could collapse together fall in one shard."
        consistency_checks: "Whether to verify that sharding kept every record: 1 for yes, 0 for no."
        truvari_matching_parameters: "Truvari collapse matching arguments."
        max_resolve: "Maximum length of a symbolic SV whose sequence Truvari resolves from the reference before collapsing."
        use_bed: "Whether Truvari collapse is restricted to each shard's intervals with a BED."
        chunk_ids_per_file: "Number of truvari-collapse shards processed on each VM."
        concat_all_naive: "Whether the genome-wide concatenation uses bcftools concat --naive: 1 for yes, 0 for no."
        ref_fa: "From references."
        ref_fai: "From references."
    }

    input {
        String remote_indir
        String remote_outdir
        Array[String] suffixes = ["ultralong", "bnd"]
        String intrasample_subdir = "01_intrasample"
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
        Int? n_expected_samples

        Int truvari_chunk_min_records = 2000
        Int truvari_collapse_refdist = 1000
        Int consistency_checks = 1
        String truvari_matching_parameters = "--refdist 500 --pctseq 0.95 --pctsize 0.95 --pctovl 0.0"
        Int max_resolve = 100000
        Boolean use_bed = false
        Int chunk_ids_per_file = 100
        Int concat_all_naive = 1

        File ref_fa
        File ref_fai

        String sv_integration_docker

        RuntimeAttr? runtime_attr_write_sample_list
        RuntimeAttr? runtime_attr_merge_calls
        RuntimeAttr? runtime_attr_shard_chromosome
        RuntimeAttr? runtime_attr_derive_chunk_ids
        RuntimeAttr? runtime_attr_collapse_shards
        RuntimeAttr? runtime_attr_concat_chromosome_shards
        RuntimeAttr? runtime_attr_concat_chromosomes
    }

    String indir_root = sub(remote_indir, "/+$", "")
    String indir = if intrasample_subdir == "" then indir_root else indir_root + "/" + intrasample_subdir
    String outdir = sub(remote_outdir, "/+$", "")

    scatter (suffix in suffixes) {
        call MergeCohortLongSVClass.MergeCohortLongSVClass {
            input:
                suffix = suffix,
                remote_indir = indir,
                remote_outdir_suffix = outdir + "/" + suffix,
                chromosomes = chromosomes,
                n_expected_samples = n_expected_samples,
                truvari_chunk_min_records = truvari_chunk_min_records,
                truvari_collapse_refdist = truvari_collapse_refdist,
                consistency_checks = consistency_checks,
                truvari_matching_parameters = truvari_matching_parameters,
                max_resolve = max_resolve,
                use_bed = use_bed,
                chunk_ids_per_file = chunk_ids_per_file,
                concat_all_naive = concat_all_naive,
                ref_fa = ref_fa,
                ref_fai = ref_fai,
                sv_integration_docker = sv_integration_docker,
                runtime_attr_write_sample_list = runtime_attr_write_sample_list,
                runtime_attr_merge_calls = runtime_attr_merge_calls,
                runtime_attr_shard_chromosome = runtime_attr_shard_chromosome,
                runtime_attr_derive_chunk_ids = runtime_attr_derive_chunk_ids,
                runtime_attr_collapse_shards = runtime_attr_collapse_shards,
                runtime_attr_concat_chromosome_shards = runtime_attr_concat_chromosome_shards,
                runtime_attr_concat_chromosomes = runtime_attr_concat_chromosomes
        }
    }

    output {
    }
}
