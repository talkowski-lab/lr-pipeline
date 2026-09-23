# Long-Read Annotation
This document describes each retired WDL workflow, including its purpose, inputs and outputs. Archived workflows are retained as historical reference only. They are not active pipeline entry points, and they are excluded from Dockstore registration and from the validation and style checks that run over `wdl/`.

This file is generated from the `meta` and `parameter_meta` blocks of each workflow by [`generate_workflows_doc.py`](../../.github/scripts/generate_workflows_doc.py). Edit those blocks rather than this document. Inputs described as `From references.` are the shared reference files listed in [references](../../docs/references.md).


## Annotations


### [AnnotateSingletonReads](../wdl/annotation/AnnotateSingletonReads.wdl)
This utility flags variants that look like single-read artifacts. Working one contig at a time it recomputes `AC`, then adds a `SINGLE_READ_SUPPORT` FILTER to any variant whose allele count is at or below two and whose alternate allele is supported by exactly one read in exactly one sample, and concatenates the per-contig results.

Inputs:
- `File vcf`: Cohort VCF to flag.
- `File vcf_idx`: Index for `vcf`.
- `Array[String] contigs`: Contigs to process.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File singleton_filtered_vcf`: VCF whose single-read-supported variants carry the `SINGLE_READ_SUPPORT` FILTER.
- `File singleton_filtered_vcf_idx`: Index for `singleton_filtered_vcf`.

### [AnnotateTruvariRemap](../wdl/annotation/AnnotateTruvariRemap.wdl)
This tool remaps insertion sequences with minimap2 (via Truvari) in order to flag insertions whose inserted sequence aligns elsewhere in the reference. Each insertion above a minimum length is realigned per contig and assessed against alignment-score and coverage thresholds, emitting a TSV of the remap results.

Inputs:
- `File vcf`: VCF whose insertions are remapped.
- `File vcf_idx`: Index for VCF.
- `File ref_fa`: From references.
- `Array[File] ref_bwa_idx`: BWA indices for `ref_fa`, from references.
- `Array[String] contigs`: Contigs to process.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String type_field`: INFO field giving each variant's allele type, used to select insertions to remap. (default `allele_type`)
- `String type_ins`: Value of `type_field` identifying an insertion. (default `ins`)
- `Int min_length`: Minimum insertion length to remap.
- `Int max_length`: Maximum insertion length to remap.
- `Int mm2_threshold`: Minimum minimap2 alignment score to flag an insertion.
- `Float cov_threshold`: Minimum alignment coverage to flag an insertion.
- `String prefix`: Prefix for output file names.
- `String remap_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotations_tsv_remap`: TSV of insertion remap results.


## Annotation Utilities


### [AnnotateCallsetOverlapWithPlotting](../wdl/annotation_utils/AnnotateCallsetOverlapWithPlotting.wdl)
This utility is an earlier form of `AnnotateCallsetOverlap` that also produces summary statistics and plots. It matches a callset VCF against SNV/indel and SV truth VCFs using exact, Truvari and `bedtools closest` rounds, each of which can be enabled on its own, and can normalize the callset, derive variant attributes and compare VEP annotations before matching.

Inputs:
- `File vcf`: Callset VCF being annotated.
- `File vcf_idx`: Index for `vcf`.
- `File truth_snv_indel_vcf`: Truth VCF containing SNVs & indels to match against.
- `File truth_snv_indel_vcf_idx`: Index for `truth_snv_indel_vcf`.
- `File truth_sv_vcf`: Truth VCF containing SVs to match against.
- `File truth_sv_vcf_idx`: Index for `truth_sv_vcf`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to evaluate.
- `Int? records_per_shard`: Number of variants to keep within a single shard during matching.
- `Boolean normalize_vcf`: Whether to normalize and split multiallelics around the VEP call. (default `false`)
- `Boolean create_variant_attributes`: Whether to derive variant attributes on the callset before matching. (default `false`)
- `Boolean compare_annotations`: Whether to compare VEP annotations between the callset and the truth callsets. (default `true`)
- `Boolean do_exact`: Whether to run the exact-match round. (default `true`)
- `Boolean do_truvari`: Whether to run the Truvari matching round. (default `true`)
- `Boolean do_bedtools_closest`: Whether to run the `bedtools closest` matching round. (default `true`)
- `Int min_sv_length_truvari`: Minimum length for a callset variant to enter the Truvari matching round.
- `Int min_sv_length_truth_truvari`: Minimum length for a truth variant to enter the Truvari matching round.
- `Int min_sv_length_bedtools_closest`: Minimum length for a callset variant to enter the `bedtools closest` matching round.
- `Int min_sv_length_truth_bedtools_closest`: Minimum length for a truth variant to enter the `bedtools closest` matching round.
- `String type_field`: INFO field in the callset VCF giving each variant's allele type. (default `allele_type`)
- `String length_field`: INFO field in the callset VCF giving each variant's allele length. (default `allele_length`)
- `String source_tag_truth_snv_indel_vcf`: Label used to tag matches against the SNV & indel truth VCF. (default `SNV_indel`)
- `String source_tag_truth_sv_vcf`: Label used to tag matches against the SV truth VCF. (default `SV`)
- `String normalize_check_ref`: `bcftools norm` `--check-ref` mode used when normalizing. (default `w`)
- `String skip_vep_categories`: VEP consequence categories excluded when comparing annotations. (default empty)
- `String af_field_sv_truth`: INFO field in the SV truth VCF holding the allele frequency. (default `AF`)
- `String ac_field_sv_truth`: INFO field in the SV truth VCF holding the allele count. (default `AC`)
- `String an_field_sv_truth`: INFO field in the SV truth VCF holding the allele number. (default `AN`)
- `String? args_string_vcf`: `bcftools view` arguments used to pre-subset the callset VCF.
- `String? args_string_truth_snv_indel_vcf`: `bcftools view` arguments used to pre-subset the SNV & indel truth VCF.
- `String? args_string_truth_sv_vcf`: `bcftools view` arguments used to pre-subset the SV truth VCF.
- `String? rename_id_string_vcf`: Expression used to rename variant IDs in the callset VCF prior to matching.
- `String? rename_id_string_truth_snv_indel_vcf`: Expression used to rename variant IDs in the SNV & indel truth VCF prior to matching.
- `String? rename_id_string_truth_sv_vcf`: Expression used to rename variant IDs in the SV truth VCF prior to matching.
- `Boolean? rename_id_strip_chr_vcf`: Whether to strip the `chr` prefix when renaming callset variant IDs.
- `Boolean? rename_id_strip_chr_truth_snv_indel_vcf`: Whether to strip the `chr` prefix when renaming SNV & indel truth variant IDs.
- `Boolean? rename_id_strip_chr_truth_sv_vcf`: Whether to strip the `chr` prefix when renaming SV truth variant IDs.
- `String prefix`: Prefix for output file names.
- `String benchmark_annotations_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (37).

Outputs:
- `File annotations_tsv_benchmark`: TSV mapping callset variants to their matched truth variants, match type, and the truth callset's AC/AF/AN and genotype-count fields.
- `File? benchmark_annotations_summary_tsv`: Summary counts of matched and unmatched variants.
- `File? benchmark_annotations_stats_tsv`: Match statistics underlying the plots.
- `File? benchmark_annotations_plots_tarball`: Tarball of the generated benchmarking plots.

### [AnnotateExternalAFs](../wdl/annotation_utils/AnnotateExternalAFs.wdl)
This utility annotates a cohort VCF with allele frequencies drawn from external reference BEDs. Each contig is subset, duplications are converted to insertions with their original type kept in a tag, the reference BEDs are matched with `bedtools closest`, and the selected matches are written back as INFO fields before the contigs are concatenated.

Inputs:
- `File vcf`: Cohort VCF to annotate.
- `File vcf_index`: Index for `vcf`.
- `Array[File] ref_beds`: External reference BEDs supplying the allele frequencies to transfer.
- `Array[String] ref_prefixes`: INFO field prefix for each entry in `ref_beds`, in the same order.
- `Array[String] contigs`: Contigs to annotate.
- `String prefix`: Prefix for output file names.
- `String pipeline_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File external_af_annotated_vcf`: VCF annotated with the external allele frequencies.
- `File external_af_annotated_vcf_index`: Index for `external_af_annotated_vcf`.

### [AnnotateExternalVariants](../wdl/annotation_utils/AnnotateExternalVariants.wdl)
This utility matches an evaluation VCF against a truth VCF by structural variant type. Both callsets are converted to BED and split into deletions, duplications and insertions, compared with `bedtools` both within type and across the duplication and insertion types, and the per-type results are combined into one TSV of matched variants.

Inputs:
- `File vcf_eval`: VCF being evaluated.
- `File vcf_eval_idx`: Index for `vcf_eval`.
- `File vcf_truth`: Truth VCF to evaluate against.
- `File vcf_truth_idx`: Index for `vcf_truth`.
- `Array[String] population`: Population labels whose allele frequencies are carried across from the truth callset.
- `String prefix`: Prefix for output file names.
- `String sv_pipeline_docker`, `String sv_base_mini_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File matched_variants_tsv`: TSV pairing each evaluation variant with its matched truth variant.

### [AnnotateSvCallerSupport](../wdl/annotation_utils/AnnotateSvCallerSupport.wdl)
This utility annotates each SV in a cohort VCF with the set of raw callers that independently support it. For every sample it matches the cohort calls against that sample's per-caller VCFs (Kanpig, cuteSV, Sniffles, Delly, pbsv, Sawfish, dipcall and hapdiff) using reciprocal-overlap, size- and sequence-similarity and a breakpoint window, then merges the support back into the cohort VCF. It outputs the annotated VCF and a TSV of per-caller match counts.

Inputs:
- `File sv_vcf`: Cohort SV VCF to annotate.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[File] kanpig_vcfs`: Per-sample Kanpig VCFs.
- `Array[File] kanpig_vcf_idxs`: Indexes for `kanpig_vcfs`.
- `Array[String] sample_ids`: Samples to process.
- `Array[File?]? sample_sv_stats`: Optional per-sample BED listing the callers supporting each variant.
- `Array[File?]? cutesv_vcfs`: Per-sample cuteSV VCFs.
- `Array[File?]? cutesv_vcf_idxs`: Indexes for `cutesv_vcfs`.
- `Array[File?]? sniffles_vcfs`: Per-sample Sniffles VCFs.
- `Array[File?]? sniffles_vcf_idxs`: Indexes for `sniffles_vcfs`.
- `Array[File?]? delly_vcfs`: Per-sample Delly VCFs.
- `Array[File?]? delly_vcf_idxs`: Indexes for `delly_vcfs`.
- `Array[File?]? pbsv_vcfs`: Per-sample pbsv VCFs.
- `Array[File?]? pbsv_vcf_idxs`: Indexes for `pbsv_vcfs`.
- `Array[File?]? sawfish_vcfs`: Per-sample Sawfish VCFs.
- `Array[File?]? sawfish_vcf_idxs`: Indexes for `sawfish_vcfs`.
- `Array[File?]? dipcall_vcfs`: Per-sample dipcall VCFs.
- `Array[File?]? dipcall_vcf_idxs`: Indexes for `dipcall_vcfs`.
- `Array[File?]? hapdiff_vcfs`: Per-sample hapdiff VCFs.
- `Array[File?]? hapdiff_vcf_idxs`: Indexes for `hapdiff_vcfs`.
- `Int truvari_breakpoint_window`: Breakpoint window, in bp, for matching a raw call. (default `500`)
- `Float truvari_reciprocal_overlap`: Minimum reciprocal overlap for matching a raw call. (default `0.0`)
- `Float truvari_sequence_similarity`: Minimum sequence similarity for matching a raw call. (default `0.7`)
- `Float truvari_size_similarity`: Minimum size similarity for matching a raw call. (default `0.7`)
- `Boolean fuzzy_match_vcf_to_stats`: Whether to match cohort records to `sample_sv_stats` by proximity rather than by exact variant ID. (default `true`)
- `Int fuzzy_match_breakpoint_window`: Breakpoint window, in bp, for fuzzy-matching a raw call to per-caller stats. (default `500`)
- `Boolean match_gt_kanpig`: Whether a Kanpig record must have a matching genotype to count as support. (default `true`)
- `Boolean match_gt_non_kanpig`: Whether a non-Kanpig caller record must have a matching genotype to count as support. (default `true`)
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `File? null_file`: Placeholder file used where an optional per-caller input is absent.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File sv_added_vcf`: Cohort VCF annotated with raw-caller support.
- `File sv_added_vcf_idx`: Index for the annotated VCF.
- `File sv_match_counts_tsv`: TSV of per-caller match counts.

### [AnnotateVcfCleared](../wdl/annotation_utils/AnnotateVcfCleared.wdl)
This utility is a variant of `AnnotateVcf` that, before applying the annotation TSVs, clears existing annotations and optionally swaps in records from an untrimmed VCF in order to restore full REF/ALT alleles. It then adds the specified INFO fields and outputs the annotated VCF.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File? subset_untrimmed_vcf`: Untrimmed VCF whose records are swapped in to restore full REF/ALT alleles.
- `File? subset_untrimmed_vcf_idx`: Index for `subset_untrimmed_vcf`.
- `Array[File] annotations_tsvs`: Annotation TSVs to apply, each as a set of INFO fields.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Array[Boolean]? sort_tsvs`: Per-TSV flag indicating whether to sort the TSV before annotation.
- `Array[String]? subset_vcf_strings`: Per-TSV `bcftools view` arguments used to pre-subset the VCF.
- `Array[String]? awk_tsv_conditions`: Per-TSV `awk` condition used to filter the TSV rows applied.
- `Array[Array[String]] info_names`: INFO field names added by each annotation TSV.
- `Array[Array[String]] info_descriptions`: INFO field header descriptions for each annotation TSV.
- `Array[Array[String]] info_types`: INFO field types for each annotation TSV.
- `Array[Array[String]] info_numbers`: INFO field `Number` values for each annotation TSV.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotated_vcf`: Annotated VCF.
- `File annotated_vcf_idx`: Index for the annotated VCF.

### [BenchmarkSTRs](../wdl/annotation_utils/BenchmarkSTRs.wdl)
This utility benchmarks per-sample TRGT tandem-repeat genotypes against a Vamos callset. Each sample's TRGT VCF is compared with the shared Vamos VCF, per-sample match statistics are collected, and the results are aggregated into genotype-concordance matrices and plots of sequence similarity, edit distance and length difference.

Inputs:
- `Array[String] sample_ids`: Sample IDs in the cohort, aligned to `trgt_vcfs`.
- `Array[File] trgt_vcfs`: Per-sample TRGT VCFs whose loci are matched against the callset.
- `Array[File] trgt_vcf_idx`: Index for `trgt_vcfs`.
- `File vamos_vcf`: Vamos callset the TRGT genotypes are benchmarked against.
- `File vamos_vcf_index`: Index for `vamos_vcf`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to benchmark.
- `String output_prefix`: Prefix for output file names.
- `Boolean include_all_regions`: Whether to additionally benchmark every locus rather than only non-reference genotypes. (default `false`)
- `String benchmark_strs_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `Array[File] benchmark_strs_per_sample_stats`: Per-sample match statistics.
- `Array[String] benchmark_strs_processed_sample_ids`: IDs of the samples that were successfully benchmarked.
- `File benchmark_strs_processed_samples_file`: File listing the samples that were successfully benchmarked.
- `File benchmark_strs_aggregated_match_data_non_ref`: Aggregated match data across samples, restricted to non-reference genotypes.
- `File benchmark_strs_genotype_concordance_matrix_non_ref`: Genotype concordance matrix for non-reference genotypes.
- `File benchmark_strs_similarity_plot_non_ref`: Sequence-similarity plot for non-reference genotypes.
- `File benchmark_strs_edit_distance_plot_non_ref`: Edit-distance plot for non-reference genotypes.
- `File benchmark_strs_length_difference_plot_non_ref`: Length-difference plot for non-reference genotypes.
- `File benchmark_strs_length_diff_vs_locus_size_non_ref`: Length difference against locus size for non-reference genotypes.
- `File? benchmark_strs_aggregated_match_data_all`: Aggregated match data across samples over all loci.
- `File? benchmark_strs_genotype_concordance_matrix_all`: Genotype concordance matrix over all loci.
- `File? benchmark_strs_similarity_plot_all`: Sequence-similarity plot over all loci.
- `File? benchmark_strs_edit_distance_plot_all`: Edit-distance plot over all loci.
- `File? benchmark_strs_length_difference_plot_all`: Length-difference plot over all loci.
- `File? benchmark_strs_length_diff_vs_locus_size_all`: Length difference against locus size over all loci.
- `File? benchmark_strs_edit_distance_to_reference_all`: Edit distance to the reference allele over all loci.
- `File? benchmark_strs_length_difference_to_reference_all`: Length difference from the reference allele over all loci.

### [CombineVcfs](../wdl/annotation_utils/CombineVcfs.wdl)
This utility concatenates two VCFs holding different variants for the same samples. The sample lists are checked for a match, each contig is subset from both inputs, and the per-contig results are concatenated into one VCF.

Inputs:
- `File a_vcf`: First VCF to combine.
- `File a_vcf_idx`: Index for `a_vcf`.
- `File b_vcf`: Second VCF to combine.
- `File b_vcf_idx`: Index for `b_vcf`.
- `Array[String] contigs`: Contigs to process.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File concat_vcf`: Combined VCF.
- `File concat_vcf_idx`: Index for the combined VCF.

### [CompareBams](../wdl/annotation_utils/CompareBams.wdl)
This utility compares two unaligned BAMs by read identity, sequence length, and sequence content. It reports total read counts, the number of reads whose IDs match across BAMs, the number of matched-ID pairs with identical sequence lengths, and the number with identical sequences (compared via MD5). It also emits a per-read TSV covering all reads from both files.

Inputs:
- `File bam1`: First unaligned BAM.
- `File bam2`: Second unaligned BAM.
- `String bam1_name`: Label for `bam1`, used as column/metric prefix in outputs.
- `String bam2_name`: Label for `bam2`, used as column/metric prefix in outputs.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File comparison_tsv`: TSV with columns `metric` and `value` reporting `{bam1_name}_total_reads`, `{bam2_name}_total_reads`, `matched_id_reads`, `matched_id_reads_same_sequence_length`, and `matched_id_reads_same_sequence`.
- `File per_read_tsv`: TSV with columns `read_id`, `{bam1_name}_len`, `{bam2_name}_len` for all reads across both BAMs. Length is empty for reads absent from that BAM.

### [CompareVcfSamples](../wdl/annotation_utils/CompareVcfSamples.wdl)
This utility compares the sample list of a VCF against a supplied list of sample IDs, reporting how many samples are shared, how many appear only in the VCF and how many appear only in the supplied list, along with the IDs in each category.

Inputs:
- `File vcf`: VCF whose samples are compared.
- `File vcf_idx`: Index for `vcf`.
- `Array[String] sample_ids`: Sample IDs to compare the VCF against.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File samples_summary_counts`: Counts of shared and unique samples.
- `File samples_common`: Sample IDs present in both the VCF and the supplied list.
- `File samples_vcf_only`: Sample IDs present only in the VCF.
- `File samples_sample_list_only`: Sample IDs present only in the supplied list.

### [CreateBiallelicVcf](../wdl/annotation_utils/CreateBiallelicVcf.wdl)
This utility normalizes a VCF into a streamlined biallelic callset. It splits multiallelic records and left-aligns variants against the reference, sorts the result, adds the `allele_length` and `allele_type` INFO fields, and rewrites each variant ID to `CHROM-POS-REF-ALT` for SNVs or `CHROM-POS-TYPE-LENGTH` otherwise, suffixing any colliding IDs to keep them unique. It outputs the biallelic VCF.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF to process.
- `File ref_fa`: Reference FASTA used for normalization.
- `File ref_fai`: Index for `ref_fa`.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File biallelic_vcf`: Normalized, sorted biallelic VCF with streamlined variant IDs and `allele_length`/`allele_type` annotations.
- `File biallelic_vcf_idx`: Index for the biallelic VCF.

### [CreateCohortDepthFiles](../wdl/annotation_utils/CreateCohortDepthFiles.wdl)
This utility ports GATK-SV's `MakeBincovMatrix` and `PloidyEstimation` workflows to build a cohort binned-coverage matrix and per-sample ploidy estimate from per-sample `MosDepth` per-base coverage BEDs. Since mosdepth's per-base output is run-length-encoded at irregular interval widths rather than GATK-SV's fixed-width `CollectReadCounts` bins, each sample's per-base BED is first binned at `bin_size` by taking the median depth per bin (dropping any trailing partial bin), matching the binning convention used by `CreateSampleReadCounts`; because every sample is binned identically, the format-detection/shift logic in upstream `MakeBincovMatrix` (which has to distinguish raw bincov BEDs from GATK `CollectReadCounts` output) is dropped as dead code. The binned files are then run through GATK-SV's `SetBins`/`MakeBincovMatrixColumns`/`ZPaste` logic to build the bincov matrix, and through `BuildPloidyMatrix` (re-binning the bincov matrix to `ploidy_bin_size`, summing depths) and GATK-SV's `estimatePloidy.R` to estimate ploidy. GATK-SV's `estimatePloidy.R` and `estimated_CN_denoising.py` are vendored under `scripts/helper/` and built into the `utils` image, so workflow has no dependency on GATK-SV docker images. Matrix outputs remain separate; `ploidy_plots` tarball contains only PNG figures from `estimatePloidy.R` and `cn_denoising_plots.pdf`. Unlike upstream `MakeBincovMatrix`, this does not support merging into a pre-existing batch's bincov matrix, since only a single one-shot cohort matrix was needed.

`estimatePloidy.R` hardcodes a 24-contig human karyotype (`chr1`..`chr22`, `chrX`, `chrY`, in that exact order) for sex assignment and per-contig ploidy expectations via positional indexing, and its 'X'/'Y' exclusion checks compare against bare `X`/`Y` rather than `chr`-prefixed names (a no-op against GRCh38-style contig names, with limited practical effect here since sample-batching/PCA (`-k`) is never invoked). `mosdepth_bed_files` must therefore be restricted to exactly those 24 contigs, in that order, or ploidy estimates will be silently wrong.

Inputs:
- `Array[String] sample_ids`: Cohort sample IDs, parallel to `mosdepth_bed_files`.
- `Array[File] mosdepth_bed_files`: Per-sample combined mosdepth per-base coverage BEDs, restricted to `chr1`-`chr22`, `chrX`, `chrY` in that order (see caveat above).
- `Int bin_size`: Size, in bp, of each coverage bin in the bincov matrix (GATK-SV convention default: 1000). (default `100`)
- `Int ploidy_bin_size`: Size, in bp, of each bin in the ploidy matrix (GATK-SV convention default: 1000000). (default `1000000`)
- `Int random_seed`: Seed for the draw, so the selection is reproducible. (default `42`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File binned_coverage`: Cohort binned-coverage matrix, bgzipped and tabix-indexed.
- `File binned_coverage_idx`: Index for `binned_coverage`.
- `File median_coverage`: Per-sample median coverage matrix.
- `File binned_estimated_ecn`: Per-sample, per-`ploidy_bin_size`-bin estimated copy number.
- `File estimated_cn`: Per-sample, per-chromosome estimated copy number.
- `File ploidy_plots`: Tarball containing only ploidy PNG and PDF figures.

### [CreateDepthProfile](../wdl/annotation_utils/CreateDepthProfile.wdl)
This utility builds a read-depth profile across one genomic window. Each sample's mosdepth BED is queried for the window and the extracted depths are combined into a single matrix with one column per sample.

Inputs:
- `Array[String] sample_ids`: Sample IDs in the cohort, aligned to `mosdepth_bed_files`.
- `Array[File] mosdepth_bed_files`: Per-sample mosdepth depth BEDs.
- `Array[File] mosdepth_bed_idx`: Indexes for `mosdepth_bed_files`.
- `String contig`: Contig containing the window.
- `Int window_start`: Start position of the window.
- `Int window_end`: End position of the window.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File region_depth_profile`: Depth matrix over the window, with one column per sample.

### [CreateFastqFromS3Reads](../wdl/annotation_utils/CreateFastqFromS3Reads.wdl)
This utility downloads BAM or FASTQ files from S3 in parallel, converts BAMs to FASTQ format preserving methylation tags, and merges all outputs into a single FASTQ.gz file.

Inputs:
- `Array[String] addresses`: S3 addresses of files to download. Supports `.bam`, `.fastq.gz`, and `.fastq` inputs.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File merged_fastq_gz`: Merged FASTQ.gz file containing reads from all input files.

### [CreatePedigreeAncestryFilesAoUPhase1](../wdl/annotation_utils/CreatePedigreeAncestryFilesAoUPhase1.wdl)
This utility generates a minimal pedigree file and an ancestry-assignment file from a list of sample IDs and their sexes, assigning every sample the `afr` ancestry of the All of Us Phase 1 cohort. It outputs both files.

Inputs:
- `Array[String] sample_ids`: Sample IDs to include.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File ped`: Generated pedigree file.
- `File ancestry`: Generated ancestry-assignment file.

### [DownloadAWSFile](../wdl/annotation_utils/DownloadAWSFile.wdl)
This utility downloads a single file from S3 and copies it to GCS, mirroring the S3 path structure relative to a configurable base prefix.

Inputs:
- `String aws_path`: S3 URI of the file to download.
- `String gcs_folder`: GCS destination folder.
- `String base_path`: S3 base prefix to strip when constructing the destination GCS path.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `String gcs_path`: GCS URI of the transferred file.

### [DropGenotypes](../wdl/annotation_utils/DropGenotypes.wdl)
This utility strips all genotype (sample) columns from a VCF, optionally sharding by record count for speed. It outputs the resulting sites-only VCF.

Inputs:
- `File vcf`: VCF whose genotypes are dropped.
- `File vcf_idx`: Index for VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File dropped_vcf`: Sites-only VCF.
- `File dropped_vcf_idx`: Index for the sites-only VCF.

### [ExtractBamRegion](../wdl/annotation_utils/ExtractBamRegion.wdl)
This utility extracts one genomic region from a BAM into a smaller indexed BAM, for inspecting or sharing a locus without moving the whole file.

Inputs:
- `File bam`: BAM to extract from.
- `File bai`: Index for `bam`.
- `Int start`: Start position of the region.
- `Int end`: End position of the region.
- `String chrom`: Contig containing the region.
- `String gatk_docker`, `String sv_pipeline_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File regional_bam`: BAM holding only the requested region.
- `File regional_bai`: Index for `regional_bam`.

### [FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl)
This utility fills missing FORMAT fields in one VCF using the values from a second, more complete VCF covering the same sites. It supports selectively copying named format fields plus toggles for filling alternate and reference genotypes, unphasing genotypes and adding PL. Sites are matched on CHROM/POS/REF/ALT, optionally also requiring a matching ID, and filling can be restricted to variants whose INFO field matches a given value. Either input can first be run through `bcftools norm`, sharded by record count so normalization never runs over a whole-contig VCF at once; normalized shards are re-concatenated with sorting (since normalization can shift a variant's position, e.g. when splitting a multiallelic) before being re-binned for matching. It outputs the refilled VCF.

Inputs:
- `File unfilled_vcf`: VCF whose FORMAT fields are filled.
- `File unfilled_vcf_idx`: Index for `unfilled_vcf`.
- `File filled_vcf`: VCF providing the FORMAT field values.
- `File filled_vcf_idx`: Index for `filled_vcf`.
- `String contig`: Contig to process.
- `File? ref_fa`: Reference FASTA used for normalization. Required if either normalize input is `true`.
- `File? ref_fai`: Index for `ref_fa`. Required if either normalize input is `true`.
- `Int? records_per_shard_normalize`: Number of variants per shard when normalizing. When set, normalization runs in parallel shards that are re-concatenated and sorted afterward.
- `Int? shard_bin_size_fill`: Region-bin size, in bp, used when sharding the contig for matching/filling.
- `Array[String] transfer_format_fields`: FORMAT fields to fill from `filled_vcf`.
- `Array[String] drop_format_fields`: FORMAT fields to drop entirely from the output (e.g. fields known to be unreliable). Cannot include `GT`.
- `Boolean fill_alt_gts`: Whether to overwrite a sample's GT in `unfilled_vcf` with `filled_vcf`'s GT when `filled_vcf`'s GT is alt-containing, regardless of the current GT in `unfilled_vcf`.
- `Boolean fill_ref_gts`: Whether to overwrite a sample's GT in `unfilled_vcf` with `filled_vcf`'s GT when `filled_vcf`'s GT is non-alt (hom-ref or no-call), regardless of the current GT in `unfilled_vcf`.
- `Boolean match_by_id`: Whether matching also requires equal variant IDs, in addition to CHROM/POS/REF/ALT.
- `Boolean unphase_gts`: Whether to unphase genotypes while filling.
- `Boolean add_missing_pl_via_ad`: Whether to add a `PL` FORMAT field derived from `AD` for genotypes that lack one.
- `Boolean expand_ad_across_alleles`: Whether to expand a fully-missing `AD` into one missing value per allele, which GLNexus writes as a bare '.' rather than '.,.'.
- `Boolean split_rnc_across_alleles`: Whether to split a merged `RNC` code into one character per allele copy, which GLNexus writes as 'MI' rather than 'M,I'.
- `Boolean normalize_unfilled_vcf`: Whether to normalize `unfilled_vcf` with `bcftools norm` before matching.
- `Boolean normalize_filled_vcf`: Whether to normalize `filled_vcf` with `bcftools norm` before matching.
- `String? subset_unfilled_vcf_field`: INFO field on `unfilled_vcf` used to limit which variants are filled. Requires `subset_unfilled_vcf_value`.
- `String? subset_unfilled_vcf_value`: Value that `subset_unfilled_vcf_field` must equal for a variant to be filled. Variants that don't match are left unfilled.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (11).

Outputs:
- `File refilled_vcf`: VCF with FORMAT fields filled.
- `File refilled_vcf_idx`: Index for the refilled VCF.

### [FillFormatFieldsBcfTools](../wdl/annotation_utils/FillFormatFieldsBcfTools.wdl)
This utility transfers FORMAT fields from a filled VCF back onto an unfilled VCF using `bcftools`, optionally unphasing genotypes, adding a `PL` field and adjusting the `EV` header Number. It is the `bcftools` counterpart to `FillFormatFields`.

Inputs:
- `File unfilled_vcf`: VCF whose FORMAT fields are filled.
- `File unfilled_vcf_idx`: Index for `unfilled_vcf`.
- `File filled_vcf`: VCF providing the FORMAT field values.
- `File filled_vcf_idx`: Index for `filled_vcf`.
- `Array[String] format_fields`: FORMAT fields to transfer from the filled VCF.
- `String? include_field`: INFO field used to limit which variants are filled. Requires `include_value`.
- `String? include_value`: Value that `include_field` must equal for a variant to be filled.
- `Boolean modify_ev_number`: Whether to rewrite the `EV` header Number so multi-caller values validate. (default `false`)
- `Boolean unphase_gts`: Whether to unphase genotypes while filling. (default `false`)
- `Boolean add_pl`: Whether to add a `PL` FORMAT field. (default `false`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File refilled_vcf`: VCF with FORMAT fields filled.
- `File refilled_vcf_idx`: Index for the refilled VCF.

### [FilterDuplicateZeroDepthReferenceBlocks](../wdl/annotation_utils/FilterDuplicateZeroDepthReferenceBlocks.wdl)
This utility cleans a single-sample gVCF by removing exact duplicate zero-depth, non-alt records, except that it retains one representative when removal would leave its start uncovered. It preserves gVCF coverage: a retained duplicate block is shortened by updating its `END` to one base before the next non-duplicate record when that record begins inside the block. This prevents cleanup from overlapping a distinct record or creating a coverage gap that GLNexus would genotype as `./.`. Singleton records, distinct records at the same coordinate, alternate genotypes, and records with non-zero or missing `MIN_DP` are retained unchanged.

Inputs:
- `File gvcf`: Single-sample gVCF to clean.
- `File gvcf_idx`: Index for `gvcf`.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File cleaned_vcf`: Cleaned gVCF.
- `File cleaned_vcf_idx`: Index for `cleaned_vcf`.

### [FilterTRGTCalls](../wdl/annotation_utils/FilterTRGTCalls.wdl)
This utility filters a TRGT tandem-repeat VCF, optionally dropping calls below a minimum repeat-unit length or length difference, or above a maximum catalog length. It outputs the filtered VCF.

Inputs:
- `File trgt_vcf`: TRGT VCF to filter.
- `File trgt_vcf_idx`: Index for the TRGT VCF.
- `Int? min_repeat_unit`: Minimum repeat-unit length to retain a call.
- `Int? min_length_diff`: Minimum length difference from the reference to retain a call.
- `Int? max_catalog_length`: Maximum catalog locus length to retain a call.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File trgt_filtered_vcf`: Filtered TRGT VCF.
- `File trgt_filtered_vcf_idx`: Index for the filtered VCF.

### [FindUntrimmedAlleles](../wdl/annotation_utils/FindUntrimmedAlleles.wdl)
This utility identifies variants in a VCF whose REF and ALT alleles retain untrimmed shared bases, producing a subset VCF of those records for use in restoring full allele representations downstream. It outputs the subset VCF.

Inputs:
- `File vcf`: VCF to scan.
- `File vcf_idx`: Index for `vcf`.
- `Array[String] contigs`: Contigs to scan within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during scanning.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File untrimmed_vcf`: VCF holding the variants whose alleles are not left-trimmed.
- `File untrimmed_vcf_idx`: Index for `untrimmed_vcf`.

### [GQCalculateCounts](../wdl/annotation_utils/GQCalculateCounts.wdl)
This utility computes GQ-stratified count tables used to derive GQ filtering cutoffs, from both a trio de novo analysis and a truth-set concordance analysis. Counts are bucketed by variant type, allele-length bin and supporting caller. For structural variants (`abs(allele_length) >= 50`), the `CALLER` column expands each call by its supporting callers using the `EV`/`BEV` FORMAT fields written by `AnnotateSvCallerSupport`: `kanpig`-backed calls are recorded under `CALLER=kanpig` with their own GQ, calls backed by other callers are split into one row per caller carrying an allelic depth in `EV` (with a per-caller GQ recomputed from that depth), and calls with no `BEV` are recorded with a blank `CALLER`. It outputs one TSV per analysis.

Inputs:
- `Array[File] vcfs`: Cohort VCFs to analyze.
- `Array[File] vcf_idxs`: Indexes for the cohort VCFs.
- `Array[File]? truth_vcfs`: Truth-set VCFs, one per input VCF, for the concordance analysis.
- `Array[File]? truth_vcf_idxs`: Indexes for the truth-set VCFs.
- `Array[Int] length_bins`: Allele-length bin boundaries defining the size buckets. (default `[0, 1, 10, 30, 50, 100, 500, 5000, 50000]`)
- `String? subset_vcf_string`: Optional `bcftools view` argument string to pre-subset each VCF.
- `File? ped`: Pedigree used to identify trios for the de novo analysis.
- `File? swap_samples_truth`: Optional sample-swap list applied to the truth VCFs.
- `Boolean run_trio_qc`: Whether to run the trio de novo analysis. (default `true`)
- `Boolean run_truth_qc`: Whether to run the truth-set concordance analysis. (default `true`)
- `Boolean skip_trv`: Whether to skip tandem-repeat variants. (default `true`)
- `Boolean drop_kanpig_supported_gq`: Whether a Kanpig-supported call also expands its `EV` callers, so co-supporting callers contribute their own genotype-quality rows. (default `false`)
- `Int min_fuzzy_match`: Minimum variant length to perform fuzzy matching for truth concordance. (default `20`)
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching deletions during truth concordance. (default `500`)
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for matching deletions during truth concordance. (default `0.7`)
- `Float del_size_similarity`: Minimum size similarity for matching deletions during truth concordance. (default `0.7`)
- `Int ins_breakpoint_window`: Breakpoint window, in bp, for matching insertions during truth concordance. (default `200`)
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for matching insertions during truth concordance. (default `0.0`)
- `Float ins_size_similarity`: Minimum size similarity for matching insertions during truth concordance. (default `0.5`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File? trio_denovo_tsv`: GQ-stratified trio de novo count table.
- `File? truth_concordance_tsv`: GQ-stratified truth-set concordance count table.

### [GQCutoffs](../wdl/annotation_utils/GQCutoffs.wdl)
This utility derives genotype-quality cutoffs by comparing a callset against trio and truth-set expectations. Trio children are identified from a PED, the callset and truth VCFs are subset to those samples, and precision and recall are tabulated across genotype-quality thresholds and variant length bins to produce a cutoff table.

Inputs:
- `Array[File] vcfs`: Per-contig callset VCFs to evaluate.
- `Array[File] vcf_idxs`: Index for `vcfs`.
- `Array[File] truth_vcfs`: Truth VCFs the callset is compared against.
- `Array[File] truth_vcf_idxs`: Index for `truth_vcfs`.
- `String? subset_vcf_string`: `bcftools view` arguments used to pre-subset the callset.
- `File ped`: Six-column PED used to identify trio children.
- `File? swap_samples_truth`: Optional sample-swap list applied to the truth VCFs.
- `Boolean skip_trv`: Whether to skip tandem-repeat variants. (default `true`)
- `Array[Int] length_bins`: Allele-length bin boundaries defining the size buckets. (default `[0, 1, 2, 6, 10, 30, 50, 100, 500, 5000, 50000]`)
- `Int min_length_heuristic_comparison`: Minimum variant length at which heuristic matching replaces exact matching. (default `30`)
- `Float del_size_similarity`: Minimum size similarity for matching deletions. (default `0.8`)
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for matching deletions. (default `0.8`)
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching deletions. (default `500`)
- `Float ins_size_similarity`: Minimum size similarity for matching insertions. (default `0.8`)
- `Int ins_breakpoint_window`: Breakpoint window, in bp, for matching insertions. (default `100`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File gq_cutoffs_tsv`: Table of genotype-quality cutoffs per variant class and length bin.

### [IntegrateHGSVCReference](../wdl/annotation_utils/IntegrateHGSVCReference.wdl)
This utility merges the HGSVC SNV, indel and SV reference callsets into one VCF. Each input is optionally sample-swapped, checked for a consistent sample list, subset per contig and tagged with its own source label before the three are merged.

Inputs:
- `File snv_vcf`: HGSVC SNV callset.
- `File snv_vcf_idx`: Index for `snv_vcf`.
- `File indel_vcf`: HGSVC indel callset.
- `File indel_vcf_idx`: Index for `indel_vcf`.
- `File sv_vcf`: HGSVC SV callset.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Sample IDs expected in every input callset.
- `File? sample_swap_list`: Two-column file mapping original sample IDs to their replacements.
- `String snv_source_tag`: Source label applied to variants from `snv_vcf`.
- `String snv_source_tag_description`: Header description for `snv_source_tag`.
- `String indel_source_tag`: Source label applied to variants from `indel_vcf`.
- `String indel_source_tag_description`: Header description for `indel_source_tag`.
- `String sv_source_tag`: Source label applied to variants from `sv_vcf`.
- `String sv_source_tag_description`: Header description for `sv_source_tag`.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File integrated_reference_vcf`: Merged reference callset.
- `File integrated_reference_vcf_idx`: Index for `integrated_reference_vcf`.

### [MakeDepthMetrics](../wdl/annotation_utils/MakeDepthMetrics.wdl)
This utility builds the cohort depth inputs used by depth-based CNV calling. Per-sample binned read counts are merged into one bgzipped matrix, and each sample's median coverage is computed from its mosdepth per-base BED with DuckDB and collected into a single table.

Inputs:
- `Array[String] sample_ids`: Sample IDs in the cohort, aligned to the per-sample inputs.
- `Array[File] binned_read_counts`: Binned read-counts file for the sample.
- `Array[File] mosdepth_per_base`: Per-contig per-base coverage (when `bin_size` is unset).
- `File duckdb`: DuckDB binary used to compute each sample's median coverage.
- `String output_prefix`: Prefix for output file names.
- `String unzip_docker`, `String sv_base_mini_docker`: Container images.

Outputs:
- `File merged_bincov`: Merged read-depth evidence and its tabix index for depth genotyping.
- `File merged_bincov_index`: Index for `merged_bincov`.
- `File median_cov`: Per-sample median coverage table.

### [MergeSites](../wdl/annotation_utils/MergeSites.wdl)
This utility merges redundant records at the site level within a VCF by collapsing near-identical deletions and insertions. Deletions are collapsed using size-, reciprocal-overlap, sequence- and sample-similarity thresholds plus a breakpoint distance, insertions using size-, sequence- and sample-similarity plus a breakpoint distance, while all other variants pass through untouched. It outputs the merged VCF.

Inputs:
- `File vcf`: VCF to merge.
- `File vcf_idx`: Index for VCF.
- `Float del_sample_similarity`: Minimum sample similarity for collapsing deletions. (default `0.5`)
- `Float ins_sample_similarity`: Minimum sample similarity for collapsing insertions. (default `0.5`)
- `Int del_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing deletions. (default `500`)
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for collapsing deletions. (default `0.0`)
- `Float del_sequence_similarity`: Minimum sequence similarity for collapsing deletions. (default `0.5`)
- `Float del_size_similarity`: Minimum size similarity for collapsing deletions. (default `0.5`)
- `Int del_size_max`: Maximum deletion size to collapse, or `-1` for no maximum. (default `50000`)
- `Int del_size_min`: Minimum deletion size to collapse. (default `0`)
- `Int ins_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing insertions. (default `200`)
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for collapsing insertions. (default `0.0`)
- `Float ins_sequence_similarity`: Minimum sequence similarity for collapsing insertions. (default `0.5`)
- `Float ins_size_similarity`: Minimum size similarity for collapsing insertions. (default `0.5`)
- `Int ins_size_max`: Maximum insertion size to collapse, or `-1` for no maximum. (default `50000`)
- `Int ins_size_min`: Minimum insertion size to collapse. (default `0`)
- `Int? shard_bin_size`: If set, shards each contig into regions of roughly this many base pairs, run in parallel.
- `File? ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File merged_vcf`: Site-merged VCF.
- `File merged_vcf_idx`: Index for the merged VCF.

### [MergeTRs](../wdl/annotation_utils/MergeTRs.wdl)
This utility merges a tandem-repeat callset into a base VCF one contig at a time and concatenates the contigs. It is an earlier form of `IntegrateTRs`.

Inputs:
- `File vcf`: Base VCF the tandem repeats are merged into.
- `File vcf_idx`: Index for `vcf`.
- `File tr_vcf`: Tandem-repeat VCF to integrate.
- `File tr_vcf_idx`: Index for the TR VCF.
- `Array[String] contigs`: Contigs to process.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File merged_vcf`: VCF combining the base and tandem-repeat callsets.
- `File merged_vcf_idx`: Index for `merged_vcf`.

### [MergeVEPAF](../wdl/annotation_utils/MergeVEPAF.wdl)
This utility combines a VEP-annotated VCF with an allele-frequency-annotated VCF, transferring the VEP consequence field onto the allele-frequency callset one contig at a time and merging the results.

Inputs:
- `File af_annotation_vcf`: VCF carrying the allele-frequency annotations.
- `File af_annotation_vcf_idx`: Index for `af_annotation_vcf`.
- `File vep_annotation_vcf`: VCF carrying the VEP annotations.
- `File vep_annotation_vcf_idx`: Index for `vep_annotation_vcf`.
- `Array[String] contigs`: Contigs to process.
- `String vep_info_field_name`: INFO field holding the VEP consequence string.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File merged_vcf`: VCF carrying both annotation sets.
- `File merged_vcf_idx`: Index for `merged_vcf`.

### [MergeVcfs](../wdl/annotation_utils/MergeVcfs.wdl)
This utility merges multiple per-contig VCFs covering the same contig into one, handling tandem-repeat and non-tandem-repeat variants separately. Non-TR variants are merged with Truvari using reciprocal-overlap, sequence-, size- and sample-similarity, a breakpoint distance and size bounds, while TR variants are merged on their identifiers, with optional region sharding. It outputs the merged VCF and a merge-summary TSV.

Inputs:
- `Array[File] contig_vcfs`: Per-callset VCFs for the contig being merged.
- `Array[File] contig_vcf_idxs`: Indexes for `contig_vcfs`.
- `String contig`: Contig being merged.
- `Int min_truvari_match`: Minimum variant length for Truvari matching. (default `20`)
- `Int truvari_breakpoint_window`: Maximum breakpoint distance, in bp, for merging non-TR variants. (default `500`)
- `Float truvari_reciprocal_overlap`: Minimum reciprocal overlap for merging non-TR variants. (default `0.0`)
- `Float truvari_sample_similarity`: Minimum sample similarity for merging non-TR variants. (default `0.0`)
- `Float truvari_sequence_similarity`: Minimum sequence similarity for merging non-TR variants. (default `0.7`)
- `Float truvari_size_similarity`: Minimum size similarity for merging non-TR variants. (default `0.7`)
- `Int truvari_size_max`: Maximum variant length Truvari will consider when collapsing. (default `50000`)
- `Int truvari_size_min`: Minimum variant length Truvari will consider when collapsing. (default `20`)
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the contig.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (10).

Outputs:
- `File merged_vcf`: Merged VCF.
- `File merged_vcf_idx`: Index for the merged VCF.
- `File merge_summary_tsv`: TSV summarizing the merge.

### [NormalizeAlleleTypes](../wdl/annotation_utils/NormalizeAlleleTypes.wdl)
This utility reclassifies `allele_type` values and records the original type in a new `allele_subtype` field. Variants with `allele_type=dup` are tested for tandemness against their duplication source (from `INFO/ORIGIN`) using two criteria: size similarity between the insertion length and the ORIGIN region length must meet the `dup_size_similarity` threshold, and the insertion POS must fall within the ORIGIN region or within `dup_breakpoint_window` bases of its breakpoints. All get `allele_subtype=tandem_dup`; those passing keep `allele_type=dup`, while those failing are set to `allele_type=ins`. Variants with `allele_type` of `complex_dup`, `dup_interspersed`, `inv_dup`, `alu_ins`, `line_ins`, `sva_ins` or `numt` are set to `allele_type=ins`, and those with `alu_del`, `line_del` or `sva_del` are set to `allele_type=del`, each recording the original value in `allele_subtype`. REF/ALT/POS are never modified. Records with other `allele_type` values are passed through unchanged. Supports optional record-count sharding.

Inputs:
- `File vcf`: VCF to transform.
- `File vcf_idx`: Index for `vcf`.
- `Int? records_per_shard`: Number of records per shard for parallel processing.
- `Int dup_breakpoint_window`: Maximum distance (bp) between insertion POS and ORIGIN breakpoints to pass the breakpoint check. (default `10`)
- `Float dup_size_similarity`: Minimum size similarity ratio (relative to the larger of the two lengths) between insertion and ORIGIN lengths. (default `0.9`)
- `Int min_dup_size`: Minimum insertion size (bp) to consider for the tandem check. (default `50`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File transformed_vcf`: VCF with revised `allele_type`/`allele_subtype`.
- `File transformed_vcf_idx`: Index for `transformed_vcf`.

### [NormalizeDuplicationOrigins](../wdl/annotation_utils/NormalizeDuplicationOrigins.wdl)
This utility resolves the relative `ORIGIN` coordinates of duplications and NUMTs into absolute genomic coordinates and annotates them back onto the VCF. `ORIGIN` values prefixed with `flank_` encode coordinates relative to a flanking window and are converted to genome-absolute positions; values already in absolute form are kept as-is. When multiple comma-separated `ORIGIN` values are present - whether flank-relative, absolute, or mixed - each is processed individually and the resulting absolute values are written back in their original order. It outputs the VCF with absolute-origin annotations.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for `vcf`.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.
- `Boolean modify_origin_header_number`: Whether to rewrite the `ORIGIN` header Number so multi-valued entries validate. (default `false`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File absolute_origin_vcf`: VCF with absolute `ORIGIN` coordinates.
- `File absolute_origin_vcf_idx`: Index for `absolute_origin_vcf`.

### [PALMERToVcf](../wdl/annotation_utils/PALMERToVcf.wdl)
This utility converts a sample's PALMER mobile-element calls into a VCF. Each mobile-element type is converted separately, and the per-type records are concatenated and sorted into one indexed VCF.

Inputs:
- `Array[File] PALMER_calls`: PALMER call files, one per entry in `mei_types`.
- `Array[String] mei_types`: Mobile-element type for each entry in `PALMER_calls`, in the same order.
- `String sample`: ID of the sample being processed.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String pipeline_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File PALMER_combined_vcf`: PALMER calls converted to VCF.
- `File PALMER_combined_vcf_idx`: Index for `PALMER_combined_vcf`.

### [ParseSVFormatFields](../wdl/annotation_utils/ParseSVFormatFields.wdl)
This utility reports the per-caller genotype qualities behind each structural variant call. For one contig it extracts each sample from the cohort VCF, looks that sample's record up in every per-caller VCF, and writes one row per call and supporting caller.

Inputs:
- `File cohort_vcf`: Cohort VCF whose calls are parsed.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[String] sample_ids`: Sample IDs to process.
- `Array[File] sample_sv_stats`: Per-sample BED listing the callers supporting each variant.
- `Array[File?] cutesv_vcfs`: Per-sample cuteSV VCFs.
- `Array[File?] cutesv_vcf_idxs`: Indexes for `cutesv_vcfs`.
- `Array[File?] sniffles_vcfs`: Per-sample Sniffles VCFs.
- `Array[File?] sniffles_vcf_idxs`: Indexes for `sniffles_vcfs`.
- `Array[File?] delly_vcfs`: Per-sample Delly VCFs.
- `Array[File?] delly_vcf_idxs`: Indexes for `delly_vcfs`.
- `Array[File?] pbsv_vcfs`: Per-sample pbsv VCFs.
- `Array[File?] pbsv_vcf_idxs`: Indexes for `pbsv_vcfs`.
- `Array[File?] sawfish_vcfs`: Per-sample Sawfish VCFs.
- `Array[File?] sawfish_vcf_idxs`: Indexes for `sawfish_vcfs`.
- `Array[File?] dipcall_vcfs`: Per-sample dipcall VCFs.
- `Array[File?] dipcall_vcf_idxs`: Indexes for `dipcall_vcfs`.
- `Array[File?] hapdiff_vcfs`: Per-sample hapdiff VCFs.
- `Array[File?] hapdiff_vcf_idxs`: Indexes for `hapdiff_vcfs`.
- `String contig`: Contig being processed.
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File gq_calls_tsv`: TSV with one row per call and supporting caller, carrying that caller's genotype quality.

### [PopulateSVFormatFields](../wdl/annotation_utils/PopulateSVFormatFields.wdl)
This utility fills per-caller support FORMAT fields on a cohort structural-variant VCF. Each sample's call is looked up in every per-caller VCF, the supporting callers and their genotype qualities are written back onto the record, and per-caller counts are tabulated by variant type and length bucket.

Inputs:
- `File cohort_vcf`: Cohort VCF to fill.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[String] sample_ids`: Sample IDs to process.
- `Array[File] sample_sv_stats`: Per-sample BED listing the callers supporting each variant.
- `Array[File?] cutesv_vcfs`: Per-sample cuteSV VCFs.
- `Array[File?] cutesv_vcf_idxs`: Indexes for `cutesv_vcfs`.
- `Array[File?] sniffles_vcfs`: Per-sample Sniffles VCFs.
- `Array[File?] sniffles_vcf_idxs`: Indexes for `sniffles_vcfs`.
- `Array[File?] delly_vcfs`: Per-sample Delly VCFs.
- `Array[File?] delly_vcf_idxs`: Indexes for `delly_vcfs`.
- `Array[File?] pbsv_vcfs`: Per-sample pbsv VCFs.
- `Array[File?] pbsv_vcf_idxs`: Indexes for `pbsv_vcfs`.
- `Array[File?] sawfish_vcfs`: Per-sample Sawfish VCFs.
- `Array[File?] sawfish_vcf_idxs`: Indexes for `sawfish_vcfs`.
- `Array[File?] dipcall_vcfs`: Per-sample dipcall VCFs.
- `Array[File?] dipcall_vcf_idxs`: Indexes for `dipcall_vcfs`.
- `Array[File?] hapdiff_vcfs`: Per-sample hapdiff VCFs.
- `Array[File?] hapdiff_vcf_idxs`: Indexes for `hapdiff_vcfs`.
- `String merge_args`: Arguments passed to the per-sample merge step. (default `--merge id`)
- `Boolean fuzzy_match_vcf_to_stats`: Whether to match cohort records to `sample_sv_stats` by proximity rather than by exact variant ID. (default `true`)
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (8).

Outputs:
- `File sv_filled_vcf`: VCF with the per-caller support fields populated.
- `File sv_filled_vcf_idx`: Index for `sv_filled_vcf`.
- `File sv_caller_counts_tsv`: Per-caller call counts by variant type and length bucket.
- `File sv_caller_source_tsv`: Per-call listing of the callers that supported it.

### [PostProcessTRLociAoU](../wdl/annotation_utils/PostProcessTRLociAoU.wdl)
The AoU counterpart of `PostProcessTRLociHPRCHGSVC` for cohorts that have a single joint-genotyped TRGT VCF and no haplotype-resolved base VCFs, so no sequence-agreement phasing is performed. For each disease-associated `TRExplorerV1` (JSON `Diseases` a non-empty array), it locates the matching entry in `trgt_catalog_bed_gz` (`TRExplorerV1` as a substring of the BED `ID=`), then uses that entry's coordinates to check the input VCF: a TRV whose `POS`/`POS+len(REF)-1` equal the BED start+1/end is treated as already present and left untouched. Otherwise it recovers the matching `trgt_vcf` record (subset and reordered to the main-VCF sample set), keeps it only when its recomputed `INFO/AC>0`, and either replaces the best-overlapping `INFO/allele_type=trv` record or, when nothing overlaps, inserts it as a new locus. Recovered records receive canonical `IntegrateTRs` IDs, `SOURCE=TRExplorer`, optional `HOMOPOLYMER_TRV`, VRS/region/in-silico/metric annotations, and flow through the shared `ApplyTRLocusUpdates` (extended to accept insert map rows) for envelope and `gnomAD_STR` assembly. Genotypes are emitted unphased and `POSTHOC_BACKBONE_PHASED` is never set.

This workflow emits no phasing audit: with no base VCFs there is nothing to phase against, so the shared `ApplyTRLocusUpdates` phasing summary (a header-only stub here) is deliberately not surfaced.

Inputs:
- `File vcf`: VCF to post-process.
- `File vcf_idx`: Index for `vcf`.
- `String contig`: Contig represented by `vcf`.
- `File trgt_vcf`: TRGT VCF whose loci are matched against the callset.
- `File trgt_vcf_idx`: Index for `trgt_vcf`.
- `File gnomad_tr_json`: TRExplorer catalog JSON; only entries with a non-empty `Diseases` array are eligible.
- `File trgt_catalog_bed_gz`: TRGT catalog BED (gzipped) whose column-4 `ID=` values bridge each `TRExplorerV1` to canonical coordinates.
- `Boolean run_flag_homopolymer_trvs`: Flag recovered TRVs whose shortest `MOTIFS` element has length one.
- `Boolean replace_gnomad_str`: Assemble `INFO/gnomAD_STR` from the catalog-match report.
- `File seqrepo_tar`: From references.
- `File simple_repeats_bed`: From references.
- `File seg_dup_bed`: From references.
- `File repeat_masked_bed`: From references.
- `String cadd_ht`: From references.
- `String pangolin_ht`: From references.
- `String phylop_ht`: From references.
- `String revel_ht`: From references.
- `String spliceai_ht`: From references.
- `String annotate_in_silico_predictors_script`: Path to the Hail script that performs the lookups (defaults to this repository's copy on `main`). (default `https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/annotation/annotate_insilico_predictors.py`)
- `String genome_build`: Reference genome build passed to the in-silico predictor annotation script. (default `GRCh38`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String vrs_docker`, `String hail_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (10).

Outputs:
- `File trv_postprocessed_vcf`: Post-processed tandem-repeat VCF.
- `File trv_postprocessed_vcf_idx`: Index for `trv_postprocessed_vcf`.
- `File trv_subsetted_vcf`: Tandem-repeat VCF subset to the recovered loci.
- `File trv_subsetted_vcf_idx`: Index for `trv_subsetted_vcf`.
- `File trv_catalog_match_tsv`: One row per contig-relevant catalog entry, sharing the HPRC/HGSVC columns plus a trailing `status` (`already_in_input_vcf`, `replaced_from_trgt`, `added_from_trgt`, `trgt_ac0_skipped`, `no_trgt_match`, `no_catalog_bed_match`, `not_eligible`).

### [PostProcessTRLociHPRCHGSVC](../wdl/annotation_utils/PostProcessTRLociHPRCHGSVC.wdl)
This utility reconciles disease-associated `TRExplorerV1` catalog loci with one integrated contig VCF. Only catalog records whose `Diseases` value is a non-empty array are eligible; records with a missing, non-array, or empty value are ignored. It uses only literal `TRExplorerV1` substring matches against `INFO/TRID`, searches unmatched catalog loci in per-sample TRGT VCFs, merges recovered loci with TRGT, drops merged calls with `AC=0`, and replaces overlapping integrated TRVs. It recomputes allele-specific `INFO/AC` after replacement ploidy normalization and again before output; a zero-AC replacement never removes its overlapping input TRV. Replacement calls receive VRS, region, and in-silico annotations; these annotations run directly on only recovered calls and are not sharded.

For each replaced, non-reference heterozygous TRGT genotype, it finds sample's matching truth/base VCF, reconstructs reference-relative sequence for both phased base haplotypes across replacement locus, and compares those sequences with both possible TRGT genotype orientations. `aligned` compares base haplotype 1 to replacement haplotype 1 and base haplotype 2 to replacement haplotype 2; `unaligned` uses crossed haplotypes. It phases only a unique lower-distance orientation when that orientation's summed edit distance is at most `max_phase_edit_distance` and its length-weighted combined edit-distance percentage is at most `max_phase_edit_distance_pct`; equality passes. It writes phased GT, sets `PS` to locus `POS`, and flags locus with `POSTHOC_BACKBONE_PHASED`. Reference, homozygous-alt, missing, and unresolved heterozygous calls remain unphased. It clears and reapplies `gnomAD_STR`, refreshes TR envelope tags, assigns replacement IDs exactly as `IntegrateTRs.SetTrVariantIds` (`contig-POS-TRV-(len(REF)-1)`, with `_1`, `_2`, ... on duplicates), and emits catalog-match and per-genotype TRV-phasing audit TSVs.

Inputs:
- `File vcf`: VCF to post-process.
- `File vcf_idx`: Index for `vcf`.
- `String contig`: Contig represented by `vcf`.
- `Array[File] trgt_vcfs`: Per-sample TRGT VCFs whose loci are matched against the callset.
- `Array[File] trgt_vcf_idxs`: Indexes for `trgt_vcfs`.
- `Array[String] sample_ids`: Cohort sample IDs in exact main-VCF and TRGT merge order; each parallel TRGT VCF must contain only its corresponding sample.
- `Array[File] base_vcfs`: Per-sample phased base VCFs used to evaluate tandem-repeat phasing.
- `Array[File] base_vcf_idxs`: Indexes for `base_vcfs`.
- `Boolean run_flag_homopolymer_trvs`: Flag recovered TRVs whose shortest `MOTIFS` element has length one.
- `Boolean run_normalize_ploidy`: Whether to normalize ploidy by sex - clearing chrY female calls, making chrX/chrY male calls hemizygous, enforcing diploidy and right-aligning unphased calls (requires `ped`).
- `Boolean replace_gnomad_str`: Assemble `INFO/gnomAD_STR` from the catalog-match report.
- `File? ped`: Cohort pedigree, used when normalizing ploidy.
- `File? swap_samples_base`: Optional whitespace-delimited raw-to-canonical sample-ID map applied when assigning cohort samples to `base_vcfs`.
- `Int max_phase_edit_distance`: Maximum allowed summed edit distance across both haplotype pairs in a unique winning orientation. (default `10`)
- `Float max_phase_edit_distance_pct`: Maximum allowed length-weighted combined edit-distance percentage across both haplotype pairs in a unique winning orientation: `100 * (distance_1 + distance_2) / (max(len(replacement_haplotype_1), len(base_haplotype_1), 1) + max(len(replacement_haplotype_2), len(base_haplotype_2), 1))`. (default `10.0`)
- `File gnomad_tr_json`: TRExplorer catalog JSON.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File seqrepo_tar`: From references.
- `File simple_repeats_bed`: From references.
- `File seg_dup_bed`: From references.
- `File repeat_masked_bed`: From references.
- `String cadd_ht`: From references.
- `String pangolin_ht`: From references.
- `String phylop_ht`: From references.
- `String revel_ht`: From references.
- `String spliceai_ht`: From references.
- `String annotate_in_silico_predictors_script`: Path to the Hail script that performs the lookups (defaults to this repository's copy on `main`). (default `https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/annotation/annotate_insilico_predictors.py`)
- `String genome_build`: Reference genome build passed to the in-silico predictor annotation script. (default `GRCh38`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String trgt_docker`, `String vrs_docker`, `String hail_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (17).

Outputs:
- `File trv_postprocessed_vcf`: Post-processed tandem-repeat VCF.
- `File trv_postprocessed_vcf_idx`: Index for `trv_postprocessed_vcf`.
- `File trv_subsetted_vcf`: Tandem-repeat VCF subset to the recovered loci.
- `File trv_subsetted_vcf_idx`: Index for `trv_subsetted_vcf`.
- `File trv_catalog_match_tsv`: Catalog-to-input/TRGT match audit, including numeric matched TRGT allele count (`0` indicates `AC=0`); rows with an input substring match leave all TRGT columns blank.
- `File trv_phasing_summary_tsv`: One row per replacement record and sample. Columns are `base_trid`, `replace_trid`, `sample_id`, input `base_gt`/`replace_gt`, base and replacement haplotype sequences, `edit_dist_aligned`, `edit_dist_unaligned`, winning-orientation `edit_dist_pct`, configured maxima, final VCF `final_gt`, and concise `status`. Distances are `sum (base_hap1 pair, base_hap2 pair)`.

### [PreprocessGregorVcf](../wdl/annotation_utils/PreprocessGregorVcf.wdl)
This utility prepares a GREGoR callset for annotation. Each contig is optionally sharded, normalized against the reference, given variant attributes and renamed variant IDs, and emitted both with genotypes and as a sites-only VCF.

Inputs:
- `File vcf`: GREGoR callset to preprocess.
- `File vcf_idx`: Index for `vcf`.
- `Array[String] contigs`: Contigs to process.
- `Int? records_per_shard`: Number of variants to keep within a single shard during preprocessing.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (8).

Outputs:
- `Array[File] full_vcf`: Per-contig preprocessed VCFs retaining genotypes.
- `Array[File] full_vcf_idx`: Index for `full_vcf`.
- `Array[File] stripped_vcf`: Per-contig preprocessed VCFs with genotypes removed.
- `Array[File] stripped_vcf_idx`: Index for `stripped_vcf`.

### [RenameVcfInfoFields](../wdl/annotation_utils/RenameVcfInfoFields.wdl)
This utility renames INFO fields in a VCF, replacing each given field string and its header description with a new one, optionally sharding by record count. It outputs the VCF with renamed INFO fields.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF.
- `Array[String] current_info_strings`: INFO field strings to replace.
- `Array[String] replace_info_strings`: Replacement INFO field strings, aligned to `current_info_strings`.
- `Array[String] replace_info_descriptions`: Replacement header descriptions, aligned to `replace_info_strings`.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File renamed_vcf`: VCF with renamed INFO fields.
- `File renamed_vcf_idx`: Index for the renamed VCF.

### [ReplaceKanpigGT](../wdl/annotation_utils/ReplaceKanpigGT.wdl)
This utility replaces the genotypes in a cohort VCF with the corresponding per-sample Kanpig calls for one contig, restricted to variants above a minimum length, and reports how many calls matched.

Inputs:
- `File vcf`: Cohort VCF for a single contig.
- `File vcf_idx`: Index for the cohort VCF.
- `Array[String] sample_ids`: Samples whose genotypes are replaced, aligned by index to 'sample_vcfs'.
- `Array[File] sample_vcfs`: Per-sample Kanpig VCFs supplying the replacement genotypes.
- `Array[File] sample_vcf_idxs`: Indices for the per-sample Kanpig VCFs.
- `String contig`: Contig the per-sample VCFs are subset to before replacement.
- `Int min_sv_length`: Minimum variant length for a genotype to be replaced.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File replaced_vcf`: Cohort VCF with replaced sample calls.
- `File replaced_vcf_idx`: Index for the updated VCF.
- `File match_counts_tsv`: Counts of replaced and unmatched calls per sample.

### [ReplaceSampleCalls](../wdl/annotation_utils/ReplaceSampleCalls.wdl)
This utility replaces the genotype calls of samples in a cohort VCF with the calls from a set of per-sample VCFs. It outputs the updated cohort VCF.

Inputs:
- `Array[File] sample_vcfs`: Per-sample VCFs providing the replacement calls.
- `Array[File] sample_vcf_idxs`: Indexes for `sample_vcfs`.
- `File cohort_vcf`: Cohort VCF whose calls are replaced.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `String prefix`: Prefix for output file names.
- `String docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File replaced_vcf`: Cohort VCF with replaced sample calls.
- `File replaced_vcf_idx`: Index for the updated VCF.

### [SubsetTRGTToCatalog](../wdl/annotation_utils/SubsetTRGTToCatalog.wdl)
This utility subsets a merged TRGT VCF down to the loci present in a given TRGT catalog BED, per contig. It outputs the catalog-restricted TRGT VCF.

Inputs:
- `File trgt_full_merged_vcf`: Merged TRGT VCF to subset.
- `File trgt_full_merged_vcf_idx`: Index for the TRGT VCF.
- `File trgt_catalog_bed_gz`: bgzipped TRGT catalog BED of loci to retain.
- `Array[String] contigs`: Contigs to process.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File trgt_merged_vcf`: Catalog-restricted TRGT VCF.
- `File trgt_merged_vcf_idx`: Index for the subset VCF.

### [SubsetVcfToContigs](../wdl/annotation_utils/SubsetVcfToContigs.wdl)
This utility subsets a VCF to a chosen set of contigs and concatenates the result. It outputs the subset VCF.

Inputs:
- `File vcf`: VCF to subset.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to retain.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File subset_contigs_vcf`: Contig-subset VCF.
- `File subset_contigs_vcf_idx`: Index for the subset VCF.

### [SubsetVcfToPerSample](../wdl/annotation_utils/SubsetVcfToPerSample.wdl)
This utility extracts a separate single-sample VCF for each requested sample from a set of cohort VCFs, optionally dropping specified fields first. It outputs the per-sample VCFs.

Inputs:
- `Array[File] cohort_vcfs`: Cohort VCFs to extract from.
- `Array[File] cohort_vcf_idxs`: Indexes for `cohort_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Samples to extract.
- `String? drop_fields`: Fields to drop from each VCF before extraction.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `Array[File] subset_vcfs`: Per-sample VCFs.
- `Array[File] subset_vcf_idxs`: Indexes for the per-sample VCFs.

### [SubsetVcfToSamples](../wdl/annotation_utils/SubsetVcfToSamples.wdl)
This utility subsets a cohort VCF to a list of samples, one contig at a time, and concatenates the results.

Inputs:
- `File vcf`: Cohort VCF to subset.
- `File vcf_idx`: Index for `vcf`.
- `Array[String] samples`: Sample IDs to retain.
- `Array[String] contigs`: Contigs to process.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File subset_samples_vcf`: VCF containing only the requested samples.
- `File subset_samples_vcf_idx`: Index for `subset_samples_vcf`.

### [UpdateGenotypes](../wdl/annotation_utils/UpdateGenotypes.wdl)
This utility rewrites the genotypes of a base VCF. It can transfer genotypes from a phased VCF, unphase or drop selected samples, normalize ploidy so male chrX and chrY calls are hemizygous and female chrY calls are cleared, and optionally drop genotypes altogether.

Inputs:
- `File base_vcf`: VCF whose genotypes are updated.
- `File base_vcf_idx`: Index for `base_vcf`.
- `File? phased_vcf`: VCF providing the phasing information.
- `File? phased_vcf_idx`: Index for `phased_vcf`.
- `Array[String] contigs`: Contigs to process.
- `Int? shard_bin_size`: If set, shards each contig into regions of roughly this many base pairs, run in parallel.
- `File ped`: Six-column PED giving each sample's sex, used to normalize ploidy.
- `Boolean transfer_genotypes`: Whether to transfer genotypes from `phased_vcf` onto the base VCF. (default `false`)
- `Boolean drop_genotypes`: Whether to strip genotypes before concatenation. (default `false`)
- `Boolean decrement_trv_ids`: Whether to decrement the numeric suffix of tandem-repeat variant IDs. (default `false`)
- `Array[String]? unphase_samples`: Samples to unphase when `run_unphase_samples` is set (defaults to empty).
- `Array[String]? drop_samples`: Sample IDs removed from the base VCF before updating.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File genotyped_vcf`: VCF with the updated genotypes.
- `File genotyped_vcf_idx`: Index for `genotyped_vcf`.

### [ValidateTRGTWithCatalog](../wdl/annotation_utils/ValidateTRGTWithCatalog.wdl)
This utility finds TRGT calls whose coordinates or motifs disagree with the catalog they were genotyped against. Each contig is sharded, checked against the catalog, and the incongruent records are concatenated into one VCF.

Inputs:
- `File trgt_vcf`: TRGT callset to validate.
- `File trgt_vcf_idx`: Index for `trgt_vcf`.
- `Array[String] contigs`: Contigs to process.
- `File trgt_catalog_bed_gz`: TRGT catalog BED (gzipped) the calls are validated against.
- `Int? records_per_shard`: Number of variants to keep within a single shard during validation.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File incongruent_vcf`: VCF holding the calls that disagree with the catalog.
- `File incongruent_vcf_idx`: Index for `incongruent_vcf`.


## Tools


### [CreateCramIndex](../wdl/tools/CreateCramIndex.wdl)
This tool indexes a CRAM with samtools and copies the resulting index next to it in Cloud Storage, for CRAMs delivered without one.

Inputs:
- `File cram`: CRAM to index.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String gcs_output_dir`: Cloud Storage directory the index is written to.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `String crai_gcs_path`: Cloud Storage path of the written index.

### [MergeWithAgglovar](../wdl/tools/MergeWithAgglovar.wdl)
This tool merges structural variant VCFs with agglovar, which clusters records by reciprocal overlap, size similarity and breakpoint offset, with optional allele matching.

Inputs:
- `Array[File] vcfs`: VCFs to merge.
- `Array[File] vcf_idxs`: Index for `vcfs`.
- `String run_agglovar_merge_script`: Path to the agglovar merge script run by the task. (default `https://raw.githubusercontent.com/talkowski-lab/gnomad-lr/main/scripts/agglovar/run_agglovar_merge.py`)
- `Float? ro_min`: Minimum reciprocal overlap for two records to cluster.
- `Float? size_ro_min`: Minimum size reciprocal overlap for two records to cluster.
- `Int? offset_max`: Maximum breakpoint offset, in bp, for two records to cluster.
- `Float? offset_prop_max`: Maximum breakpoint offset as a proportion of variant length.
- `Boolean match_ref`: Whether the reference alleles must match. (default `false`)
- `Boolean match_alt`: Whether the alternate alleles must match. (default `false`)
- `Float? match_prop_min`: Minimum proportion of matching allele sequence.
- `String prefix`: Prefix for output file names.
- `String agglovar_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File merged_vcf`: Merged callset.
- `File merged_vcf_index`: Index for `merged_vcf`.

### [MergeWithTruvari](../wdl/tools/MergeWithTruvari.wdl)
This tool merges VCFs by combining them with `bcftools merge` and then collapsing redundant records with Truvari (https://github.com/ACEnglish/truvari). An optional preprocessing script can reshape the merged VCF before collapsing.

Inputs:
- `Array[File] vcfs`: VCFs to merge.
- `Array[File] vcf_idxs`: Index for `vcfs`.
- `String? truvari_params`: Arguments passed to `truvari collapse`.
- `String? bcftools_merge_params`: Arguments passed to `bcftools merge`.
- `File? preprocess_script`: Script run on the merged VCF before collapsing.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String merge_docker`, `String truvari_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File truvari_collapsed_vcf`: Merged and collapsed callset.
- `File truvari_collapsed_vcf_idx`: Index for `truvari_collapsed_vcf`.

### [PhaseCallsetCommon](../wdl/tools/PhaseCallsetCommon.wdl)
This tool statistically phases one contig of a callset with SHAPEIT. The VCF is split and filtered, given unique IDs and normalized, optionally deduplicated by phased fraction and resolved for variant collisions, then phased either with SHAPEIT4 alone or with a SHAPEIT4 common-variant scaffold that SHAPEIT5 fills in with rare variants. Phase sets from the input are transferred back onto the phased output.

Inputs:
- `File vcf`: Callset VCF to phase.
- `File vcf_idx`: Index for `vcf`.
- `String contig`: Contig being phased.
- `Int operation`: Collision-resolution mode: `0` removes an entire VCF record, `1` removes single alleles from a genotype.
- `String weight_tag`: ID of the field holding each record's collision weight, so preferred records survive a collision.
- `Int is_weight_format_field`: Where `weight_tag` is read from: `0` for the INFO field, `1` for the sample column.
- `Float default_weight`: Weight assigned when `weight_tag` is absent from a record.
- `Boolean do_shapeit5`: Whether to phase rare variants with SHAPEIT5 against a SHAPEIT4 common-variant scaffold, rather than phasing everything with SHAPEIT4.
- `Boolean remove_duplicates_by_phased_fraction`: Whether to drop duplicate records, keeping the copy phased in the most samples.
- `Float min_af_common`: Minimum allele frequency for a variant to enter the common-variant scaffold.
- `String variant_filter_args`: Arguments used to filter variants before phasing. (default `-i 'MAC>=2'`)
- `String filter_common_args`: Arguments used to select the common variants for the scaffold. (default `-i 'MAF>=0.001'`)
- `String chunk_extra_args`: Extra arguments passed when creating the SHAPEIT chunks. (default `--thread $(nproc) --window-size 2000000 --buffer-size 200000`)
- `String shapeit4_extra_args`: Extra arguments passed to SHAPEIT4. (default `--thread $(nproc) --use-PS 0.0001`)
- `String shapeit5_extra_args`: Extra arguments passed to SHAPEIT5. (default `--thread $(nproc)`)
- `File genetic_maps_tsv`: TSV mapping each contig to its genetic map.
- `File fix_variant_collisions_java`: Compiled Java program that resolves variant collisions.
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String pysam_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (14).

Outputs:
- `File uqids_split_vcf`: Split VCF with unique variant IDs and normalized records.
- `File uqids_split_vcf_idx`: Index for `uqids_split_vcf`.
- `File? removed_duplicates_split_vcf`: Split VCF after duplicate records were dropped.
- `File? removed_duplicates_split_vcf_idx`: Index for `removed_duplicates_split_vcf`.
- `File collisionless_split_vcf`: Split VCF after variant collisions were resolved.
- `File collisionless_split_vcf_idx`: Index for `collisionless_split_vcf`.
- `File ps_anchors_vcf`: VCF of the phase-set anchor variants.
- `File ps_anchors_vcf_idx`: Index for `ps_anchors_vcf`.
- `File shapeit_phased_vcf`: Statistically phased callset.
- `File shapeit_phased_vcf_idx`: Index for `shapeit_phased_vcf`.
- `File shapeit_phased_ps_transferred_vcf`: Phased callset with the input phase sets transferred back on.
- `File shapeit_phased_ps_transferred_vcf_idx`: Index for `shapeit_phased_ps_transferred_vcf`.

### [PhaseCallsetWithBackbone](../wdl/tools/PhaseCallsetWithBackbone.wdl)
This tool transfers phasing from a backbone VCF onto a callset. Both callsets are reduced to their overlapping samples and to SNVs, the callset is split, deduplicated and resolved for variant collisions, and each phase set is oriented to whichever assignment agrees with the backbone's phased heterozygous calls.

Inputs:
- `File vcf`: Callset VCF to phase.
- `File vcf_idx`: Index for `vcf`.
- `File base_vcf`: Phased backbone VCF supplying the haplotype assignments.
- `File base_vcf_idx`: Index for `base_vcf`.
- `Int operation`: Collision-resolution mode: `0` removes an entire VCF record, `1` removes single alleles from a genotype.
- `String weight_tag`: ID of the field holding each record's collision weight, so preferred records survive a collision.
- `Int is_weight_format_field`: Where `weight_tag` is read from: `0` for the INFO field, `1` for the sample column.
- `Float default_weight`: Weight assigned when `weight_tag` is absent from a record.
- `Boolean remove_duplicates_by_phased_fraction`: Whether to drop duplicate records, keeping the copy phased in the most samples.
- `String variant_filter_args`: Arguments used to filter variants before phasing. (default `-i 'MAC>=2'`)
- `File fix_variant_collisions_java`: Compiled Java program that resolves variant collisions.
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String pysam_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (9).

Outputs:
- `File uqids_split_vcf`: Split VCF with unique variant IDs and normalized records.
- `File uqids_split_vcf_idx`: Index for `uqids_split_vcf`.
- `File? removed_duplicates_split_vcf`: Split VCF after duplicate records were dropped.
- `File? removed_duplicates_split_vcf_idx`: Index for `removed_duplicates_split_vcf`.
- `File collisionless_split_vcf`: Split VCF after variant collisions were resolved.
- `File collisionless_split_vcf_idx`: Index for `collisionless_split_vcf`.
- `File base_prepared_vcf`: Backbone VCF reduced to the overlapping samples and SNVs.
- `File base_prepared_vcf_idx`: Index for `base_prepared_vcf`.
- `File base_transferred_vcf`: Callset with the backbone haplotypes transferred on.
- `File base_transferred_vcf_idx`: Index for `base_transferred_vcf`.

### [PreprocessStatisticalPhasing](../wdl/tools/PreprocessStatisticalPhasing.wdl)
This tool prepares a callset for statistical phasing. Kanpig score annotations are added so collision resolution can prefer those records, selected INFO annotations are optionally removed, and calls overlapping TRGT loci are dropped.

Inputs:
- `File vcf`: Callset VCF to prepare.
- `File vcf_idx`: Index for `vcf`.
- `Boolean remove_annotations`: Whether to remove the INFO fields named in `annotations_to_remove`.
- `String? annotations_to_remove`: Comma-separated INFO fields removed when `remove_annotations` is set.
- `String prefix`: Prefix for output file names.
- `String docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File annotated_vcf`: Annotated VCF.
- `File annotated_vcf_idx`: Index for the annotated VCF.
- `File filtered_vcf`: VCF with low-coverage genotypes set to missing.
- `File filtered_vcf_idx`: Index for `filtered_vcf`.

### [TRGTMerge](../wdl/tools/TRGTMerge.wdl)
This tool merges per-sample TRGT VCFs into a cohort callset with `trgt merge`, one contig at a time, then concatenates the contigs.

Inputs:
- `Array[File] vcfs`: Per-sample TRGT VCFs to merge.
- `Array[File] vcf_idxs`: Index for `vcfs`.
- `Array[String] contigs`: Contigs to process.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String trgt_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File trgt_merged_vcf`: Catalog-restricted TRGT VCF.
- `File trgt_merged_vcf_idx`: Index for the subset VCF.

### [TransferMethylationTags](../wdl/tools/TransferMethylationTags.wdl)
This tool transfers methylation base-modification tags (MM/ML) from unaligned BAMs onto an aligned BAM. It extracts the tags per read, then per contig re-attaches them to the aligned reads and sorts, merging the result into a single tagged BAM. It outputs the methylation-tagged BAM and a TSV of the transferred tags.

Inputs:
- `File aligned_bam`: Aligned BAM to receive the tags.
- `File aligned_bai`: Index for `aligned_bam`.
- `Array[String] contigs`: Contigs to process.
- `Array[String] unaligned_bam_paths`: Paths to the unaligned BAMs carrying the methylation tags.
- `Boolean gcs_paths`: Whether `unaligned_bam_paths` are GCS paths. (default `false`)
- `Boolean recreate_bam`: Whether to rebuild the aligned BAM from the tagged reads rather than tagging it in place. (default `false`)
- `String mm_tag`: Base-modification tag name. (default `MM`)
- `String ml_tag`: Modification-likelihood tag name. (default `ML`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File methylation_tagged_bam`: Aligned BAM with methylation tags transferred.
- `File methylation_tagged_bai`: Index for the tagged BAM.
- `File? methylation_tags`: TSV of the transferred methylation tags.

### [VcfDist](../wdl/tools/VcfDist.wdl)
This tool runs vcfdist (https://github.com/TimD1/vcfdist) in order to benchmark an evaluation VCF against a truth VCF per contig, computing alignment-based precision/recall and phasing accuracy. It outputs vcfdist's precision-recall, phasing, switch-flip, phase-block and supercluster reports.

Inputs:
- `File vcf_eval`: VCF being evaluated.
- `File vcf_eval_idx`: Index for `vcf_eval`.
- `File vcf_truth`: Truth VCF to evaluate against.
- `File vcf_truth_idx`: Index for `vcf_truth`.
- `File ref_fa`: From references.
- `Array[String] contigs`: Contigs to evaluate.
- `File? bed_regions`: BED of regions to restrict the evaluation to.
- `String? mode`: vcfdist evaluation mode.
- `Float? threshold`: vcfdist matching threshold.
- `String? vcfdist_args`: Additional arguments passed to vcfdist.
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String vcfdist_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `Array[File] vcfdist_phasing_summary_tsv`: Per-contig phasing summaries.
- `Array[File] vcfdist_switchflips_tsv`: Per-contig switch and flip errors.
- `Array[File] vcfdist_precision_recall_tsv`: Per-contig precision-recall curves.
- `Array[File] vcfdist_precision_recall_summary_tsv`: Per-contig precision-recall summaries.
- `Array[File] vcfdist_phase_blocks_tsv`: Per-contig phase blocks.
- `Array[File] vcfdist_superclusters_tsv`: Per-contig variant superclusters.
- `Array[File] vcfdist_query_tsv`: Per-contig query-variant results.
- `Array[File] vcfdist_truth_tsv`: Per-contig truth-variant results.
- `Array[File] vcfdist_summary_vcf`: Per-contig annotated summary VCFs.

### [VcfDistCohort](../wdl/tools/VcfDistCohort.wdl)
This tool runs vcfdist (https://github.com/TimD1/vcfdist) across a cohort by pairing each evaluation VCF with its corresponding truth VCF and benchmarking every assigned sample, then aggregating the per-sample results. It outputs cohort-level precision/recall and phasing summaries.

Inputs:
- `Array[File] eval_vcfs`: Evaluation VCFs, one per group.
- `Array[File] eval_vcf_idxs`: Indexes for `eval_vcfs`.
- `Array[File] truth_vcfs`: Truth VCFs, aligned to `eval_vcfs`.
- `Array[File] truth_vcf_idxs`: Indexes for `truth_vcfs`.
- `File ref_fa`: From references.
- `Array[String] contigs`: Contigs to evaluate.
- `Array[String]? subset_samples`: Samples to restrict the evaluation to.
- `String? vcfdist_args`: Additional arguments passed to vcfdist.
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String vcfdist_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (11).

Outputs:
- `File vcfdist_phasing_summary_tsv`: Cohort phasing summary.
- `File vcfdist_precision_recall_summary_tsv`: Cohort precision-recall summary.
- `File vcfdist_precision_recall_tsv`: Cohort precision-recall curves.
- `File vcfdist_switchflips_tsv`: Cohort switch and flip errors.
- `File vcfdist_phase_blocks_tsv`: Cohort phase blocks.
- `File vcfdist_missing_samples`: Samples with no matching truth VCF.

### [Whatshap](../wdl/tools/Whatshap.wdl)
This tool haplotags a sample's BAM against a phased VCF using WhatsHap (https://github.com/whatshap/whatshap), per contig, then merges the tagged reads into a single BAM. It outputs the haplotagged BAM and per-contig haplotag read lists.

Inputs:
- `File bam`: Aligned reads to haplotag.
- `File bai`: Index for `bam`.
- `File phased_vcf`: Phased VCF used to assign haplotypes.
- `File phased_vcf_idx`: Index for `phased_vcf`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to process.
- `String? extra_args`: Additional arguments passed to WhatsHap.
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String whatshap_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File haplotagged_bam`: Haplotagged BAM.
- `File haplotagged_bai`: Index for the haplotagged BAM.
- `Array[File] haplotag_lists`: Per-contig haplotag read assignments.
