# Long-Read Annotation

## Annotation Utilities

### [CombineVcfsAcrossContigs](../wdl/annotation_utils/CombineVcfsAcrossContigs.wdl)
This utility concatenates a set of per-contig VCFs into a single VCF, optionally dropping genotypes in the process. It outputs the combined VCF.

Inputs:
- `Array[File] vcfs`: Per-contig VCFs to concatenate.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs corresponding to `vcfs`.
- `Boolean drop_genotypes`: Whether to strip genotypes from the combined VCF (default `false`).

Outputs:
- `concat_vcf`: Combined VCF.
- `concat_vcf_idx`: Index for the combined VCF.


### [CompareBams](../wdl/annotation_utils/CompareBams.wdl)
This utility compares two unaligned BAMs by read identity, sequence length, and sequence content. It reports total read counts, the number of reads whose IDs match across BAMs, the number of matched-ID pairs with identical sequence lengths, and the number with identical sequences (compared via MD5). It also emits a per-read TSV covering all reads from both files.

Inputs:
- `File bam1`: First unaligned BAM.
- `File bam2`: Second unaligned BAM.
- `String bam1_name`: Label for `bam1`, used as column/metric prefix in outputs.
- `String bam2_name`: Label for `bam2`, used as column/metric prefix in outputs.
- `String prefix`: Prefix for output file names.

Outputs:
- `comparison_tsv`: TSV with columns `metric` and `value` reporting `{bam1_name}_total_reads`, `{bam2_name}_total_reads`, `matched_id_reads`, `matched_id_reads_same_sequence_length`, and `matched_id_reads_same_sequence`.
- `per_read_tsv`: TSV with columns `read_id`, `{bam1_name}_len`, `{bam2_name}_len` for all reads across both BAMs. Length is empty for reads absent from that BAM.


### [CreateBiallelicVcf](../wdl/annotation_utils/CreateBiallelicVcf.wdl)
This utility normalizes a VCF into a streamlined biallelic callset. It splits multiallelic records and left-aligns variants against the reference, sorts the result, adds the `allele_length` and `allele_type` INFO fields, and rewrites each variant ID to `CHROM-POS-REF-ALT` for SNVs or `CHROM-POS-TYPE-LENGTH` otherwise, suffixing any colliding IDs to keep them unique. It outputs the biallelic VCF.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF to process.
- `File ref_fa`: Reference FASTA used for normalization.
- `File ref_fai`: Index for `ref_fa`.

Outputs:
- `biallelic_vcf`: Normalized, sorted biallelic VCF with streamlined variant IDs and `allele_length`/`allele_type` annotations.
- `biallelic_vcf_idx`: Index for the biallelic VCF.


### [CreatePedigreeAndAncestryFiles](../wdl/annotation_utils/CreatePedigreeAndAncestryFiles.wdl)
This utility generates a minimal pedigree file and an ancestry-assignment file from a list of sample IDs and their sexes. It outputs both files.

Inputs:
- `Array[String] sample_ids`: Sample IDs to include.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.

Outputs:
- `ped`: Generated pedigree file.
- `ancestry`: Generated ancestry-assignment file.


### [DownloadAWSFile](../wdl/annotation_utils/DownloadAWSFile.wdl)
This utility downloads a single file from S3 and copies it to GCS, mirroring the S3 path structure relative to a configurable base prefix.

Inputs:
- `String aws_path`: S3 URI of the file to download.
- `String gcs_folder`: GCS destination folder.
- `String base_path`: S3 base prefix to strip when constructing the destination GCS path.

Outputs:
- `String gcs_path`: GCS URI of the transferred file.


### [CreateFastqFromS3Reads](../wdl/annotation_utils/CreateFastqFromS3Reads.wdl)
This utility downloads BAM or FASTQ files from S3 in parallel, converts BAMs to FASTQ format preserving methylation tags, and merges all outputs into a single FASTQ.gz file.

Inputs:
- `Array[String] addresses`: S3 addresses of files to download. Supports `.bam`, `.fastq.gz`, and `.fastq` inputs.
- `String prefix`: Prefix for output file names.

Outputs:
- `merged_fastq_gz`: Merged FASTQ.gz file containing reads from all input files.


### [FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl)
This utility fills missing FORMAT fields in one VCF using the values from a second, more complete VCF covering the same sites. It supports selectively copying named format fields plus toggles for filling alternate and reference genotypes, unphasing genotypes and adding PL, with optional region sharding. It outputs the refilled VCF.

Inputs:
- `File unfilled_vcf`: VCF whose FORMAT fields are filled.
- `File unfilled_vcf_idx`: Index for `unfilled_vcf`.
- `File filled_vcf`: VCF providing the FORMAT field values.
- `File filled_vcf_idx`: Index for `filled_vcf`.
- `String contig`: Contig to process.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the contig.
- `Array[String] format_fields`: FORMAT fields to fill.
- `String? include_field`: FORMAT field whose value gates whether a record is filled.
- `String? include_value`: Value of `include_field` required for a record to be filled.
- `Boolean fill_alt_gts`: Whether to fill alternate-allele genotypes (default `false`).
- `Boolean fill_ref_gts`: Whether to fill reference genotypes (default `false`).
- `Boolean unphase_gts`: Whether to unphase genotypes while filling (default `false`).
- `Boolean add_pl`: Whether to add a PL field (default `false`).

Outputs:
- `refilled_vcf`: VCF with FORMAT fields filled.
- `refilled_vcf_idx`: Index for the refilled VCF.


### [FilterTRGTCalls](../wdl/annotation_utils/FilterTRGTCalls.wdl)
This utility filters a TRGT tandem-repeat VCF, optionally dropping calls below a minimum repeat-unit length or length difference, or above a maximum catalog length. It outputs the filtered VCF.

Inputs:
- `File trgt_vcf`: TRGT VCF to filter.
- `File trgt_vcf_idx`: Index for the TRGT VCF.
- `Int? min_repeat_unit`: Minimum repeat-unit length to retain a call.
- `Int? min_length_diff`: Minimum length difference from the reference to retain a call.
- `Int? max_catalog_length`: Maximum catalog locus length to retain a call.

Outputs:
- `trgt_filtered_vcf`: Filtered TRGT VCF.
- `trgt_filtered_vcf_idx`: Index for the filtered VCF.


### [FindUntrimmedAlleles](../wdl/annotation_utils/FindUntrimmedAlleles.wdl)
This utility identifies variants in a VCF whose REF and ALT alleles retain untrimmed shared bases, producing a subset VCF of those records for use in restoring full allele representations downstream. It outputs the subset VCF.

Inputs:
- `File vcf`: VCF to scan.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to scan within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during scanning.

Outputs:
- `subset_untrimmed_vcf`: VCF of records with untrimmed alleles.
- `subset_untrimmed_vcf_idx`: Index for the subset VCF.


### [GQCalculateCounts](../wdl/annotation_utils/GQCalculateCounts.wdl)
This utility computes GQ-stratified count tables used to derive GQ filtering cutoffs, from both a trio de novo analysis and a truth-set concordance analysis. Counts are bucketed by variant type, allele-length bin and supporting caller. For structural variants (`abs(allele_length) >= 50`), the `CALLER` column expands each call by its supporting callers using the `EV`/`BEV` FORMAT fields written by [AnnotateSvCallerSupport](#annotatesvcallersupport): `kanpig`-backed calls are recorded under `CALLER=kanpig` with their own GQ, calls backed by other callers are split into one row per caller carrying an allelic depth in `EV` (with a per-caller GQ recomputed from that depth), and calls with no `BEV` are recorded with a blank `CALLER`. It outputs one TSV per analysis.

Inputs:
- `Array[File] vcfs`: Cohort VCFs to analyze.
- `Array[File] vcf_idxs`: Indexes for the cohort VCFs.
- `Array[File]? truth_vcfs`: Truth-set VCFs, one per input VCF, for the concordance analysis.
- `Array[File]? truth_vcf_idxs`: Indexes for the truth-set VCFs.
- `Array[Int] length_bins`: Allele-length bin boundaries defining the size buckets.
- `String? subset_vcf_string`: Optional `bcftools view` argument string to pre-subset each VCF.
- `File? ped`: Pedigree used to identify trios for the de novo analysis.
- `File? swap_samples_truth`: Optional sample-swap list applied to the truth VCFs.
- `Boolean run_trio_qc`: Whether to run the trio de novo analysis.
- `Boolean run_truth_qc`: Whether to run the truth-set concordance analysis.
- `Boolean skip_trv`: Whether to skip tandem-repeat variants.
- `Int min_fuzzy_match`: Minimum variant length to perform fuzzy matching for truth concordance (default `20`).
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching deletions during truth concordance (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for matching deletions during truth concordance (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for matching deletions during truth concordance (default `0.7`).
- `Int ins_breakpoint_window`: Breakpoint window, in bp, for matching insertions during truth concordance (default `200`).
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for matching insertions during truth concordance (default `0.0`).
- `Float ins_size_similarity`: Minimum size similarity for matching insertions during truth concordance (default `0.5`).

Outputs:
- `trio_denovo_tsv`: GQ-stratified trio de novo count table.
- `truth_concordance_tsv`: GQ-stratified truth-set concordance count table.


### [MergeSites](../wdl/annotation_utils/MergeSites.wdl)
This utility merges redundant records at the site level within a VCF by collapsing near-identical deletions and insertions. Deletions are collapsed using size-, reciprocal-overlap, sequence- and sample-similarity thresholds plus a breakpoint distance, insertions using size-, sequence- and sample-similarity plus a breakpoint distance, while all other variants pass through untouched. It outputs the merged VCF.

Inputs:
- `File vcf`: VCF to merge.
- `File vcf_idx`: Index for VCF.
- `Int del_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing deletions (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for collapsing deletions (default `0.0`).
- `Float del_sample_similarity`: Minimum sample similarity for collapsing deletions (default `0.5`).
- `Float del_sequence_similarity`: Minimum sequence similarity for collapsing deletions (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for collapsing deletions (default `0.7`).
- `Int del_size_max`: Maximum deletion size to collapse, or `-1` for no maximum (default `-1`).
- `Int del_size_min`: Minimum deletion size to collapse (default `0`).
- `Int ins_breakpoint_window`: Maximum breakpoint distance, in bp, for collapsing insertions (default `200`).
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for collapsing insertions (default `0.0`).
- `Float ins_sample_similarity`: Minimum sample similarity for collapsing insertions (default `0.5`).
- `Float ins_sequence_similarity`: Minimum sequence similarity for collapsing insertions (default `0.7`).
- `Float ins_size_similarity`: Minimum size similarity for collapsing insertions (default `0.7`).
- `Int ins_size_max`: Maximum insertion size to collapse, or `-1` for no maximum (default `-1`).
- `Int ins_size_min`: Minimum insertion size to collapse (default `0`).

Outputs:
- `merged_vcf`: Site-merged VCF.
- `merged_vcf_idx`: Index for the merged VCF.


### [MergeVcfs](../wdl/annotation_utils/MergeVcfs.wdl)
This utility merges multiple per-contig VCFs covering the same contig into one, handling tandem-repeat and non-tandem-repeat variants separately. Non-TR variants are merged with Truvari using reciprocal-overlap, sequence-, size- and sample-similarity, a breakpoint distance and size bounds, while TR variants are merged on their identifiers, with optional region sharding. It outputs the merged VCF and a merge-summary TSV.

Inputs:
- `Array[File] contig_vcfs`: Per-callset VCFs for the contig being merged.
- `Array[File] contig_vcf_idxs`: Indexes for `contig_vcfs`.
- `String contig`: Contig being merged.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the contig.
- `Int min_truvari_match`: Minimum variant length for Truvari matching (default `20`).
- `Int truvari_breakpoint_window`: Maximum breakpoint distance, in bp, for merging non-TR variants (default `500`).
- `Float truvari_reciprocal_overlap`: Minimum reciprocal overlap for merging non-TR variants (default `0.0`).
- `Float truvari_sample_similarity`: Minimum sample similarity for merging non-TR variants (default `0.0`).
- `Float truvari_sequence_similarity`: Minimum sequence similarity for merging non-TR variants (default `0.7`).
- `Float truvari_size_similarity`: Minimum size similarity for merging non-TR variants (default `0.7`).
- `Int size_min`: Minimum variant size to merge (default `20`).
- `Int size_max`: Maximum variant size to merge (default `50000`).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `merged_vcf`: Merged VCF.
- `merged_vcf_idx`: Index for the merged VCF.
- `merge_summary_tsv`: TSV summarizing the merge.


### [NormalizeDuplicationOrigins](../wdl/annotation_utils/NormalizeDuplicationOrigins.wdl)
This utility resolves the relative `ORIGIN` coordinates of duplications and NUMTs into absolute genomic coordinates and annotates them back onto the VCF. `ORIGIN` values prefixed with `flank_` encode coordinates relative to a flanking window and are converted to genome-absolute positions; values already in absolute form are kept as-is. When multiple comma-separated `ORIGIN` values are present - whether flank-relative, absolute, or mixed - each is processed individually and the resulting absolute values are written back in their original order. It outputs the VCF with absolute-origin annotations.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.

Outputs:
- `absolute_origin_vcf`: VCF with absolute `ORIGIN` coordinates.
- `absolute_origin_vcf_idx`: Index for the annotated VCF.


### [RenameVcfInfoFields](../wdl/annotation_utils/RenameVcfInfoFields.wdl)
This utility renames INFO fields in a VCF, replacing each given field string and its header description with a new one, optionally sharding by record count. It outputs the VCF with renamed INFO fields.

Inputs:
- `File vcf`: VCF to process.
- `File vcf_idx`: Index for VCF.
- `Array[String] current_info_strings`: INFO field strings to replace.
- `Array[String] replace_info_strings`: Replacement INFO field strings, aligned to `current_info_strings`.
- `Array[String] replace_info_descriptions`: Replacement header descriptions, aligned to `replace_info_strings`.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.

Outputs:
- `renamed_vcf`: VCF with renamed INFO fields.
- `renamed_vcf_idx`: Index for the renamed VCF.


### [ReplaceSampleCalls](../wdl/annotation_utils/ReplaceSampleCalls.wdl)
This utility replaces the genotype calls of samples in a cohort VCF with the calls from a set of per-sample VCFs. It outputs the updated cohort VCF.

Inputs:
- `Array[File] sample_vcfs`: Per-sample VCFs providing the replacement calls.
- `Array[File] sample_vcf_idxs`: Indexes for `sample_vcfs`.
- `File cohort_vcf`: Cohort VCF whose calls are replaced.
- `File cohort_vcf_idx`: Index for the cohort VCF.

Outputs:
- `replaced_vcf`: Cohort VCF with replaced sample calls.
- `replaced_vcf_idx`: Index for the updated VCF.


### [SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)
This utility splits a VCF into per-contig VCFs, optionally also producing genotype-free copies and applying fixups such as adding missing INFO header lines, modifying SNV IDs and renaming dbSNP/dbVar contigs. It outputs the per-contig VCFs and their no-genotype counterparts.

Inputs:
- `File vcf`: VCF to split.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to split the VCF into.
- `Boolean create_no_geno`: Whether to also produce genotype-free copies (default `false`).
- `Boolean modify_snv_ids`: Whether to rewrite SNV variant IDs (default `false`).
- `Boolean rename_dbsnp_contigs`: Whether to rename contigs to dbSNP naming (default `false`).
- `Boolean rename_dbvar_contigs`: Whether to rename contigs to dbVar naming (default `false`).
- `Array[String]? missing_info_header_fields`: INFO header lines to add if missing.

Outputs:
- `contig_vcfs`: Per-contig VCFs.
- `contig_vcf_idxs`: Indexes for the per-contig VCFs.
- `contig_no_geno_vcfs`: Per-contig genotype-free VCFs.
- `contig_no_geno_vcf_idxs`: Indexes for the genotype-free VCFs.


### [SubsetTRGTToCatalog](../wdl/annotation_utils/SubsetTRGTToCatalog.wdl)
This utility subsets a merged TRGT VCF down to the loci present in a given TRGT catalog BED, per contig. It outputs the catalog-restricted TRGT VCF.

Inputs:
- `File trgt_full_merged_vcf`: Merged TRGT VCF to subset.
- `File trgt_full_merged_vcf_idx`: Index for the TRGT VCF.
- `File trgt_catalog_bed_gz`: bgzipped TRGT catalog BED of loci to retain.
- `Array[String] contigs`: Contigs to process.

Outputs:
- `trgt_merged_vcf`: Catalog-restricted TRGT VCF.
- `trgt_merged_vcf_idx`: Index for the subset VCF.


### [SubsetVcfToContigs](../wdl/annotation_utils/SubsetVcfToContigs.wdl)
This utility subsets a VCF to a chosen set of contigs and concatenates the result. It outputs the subset VCF.

Inputs:
- `File vcf`: VCF to subset.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to retain.

Outputs:
- `subset_contigs_vcf`: Contig-subset VCF.
- `subset_contigs_vcf_idx`: Index for the subset VCF.


### [SubsetVcfToPerSample](../wdl/annotation_utils/SubsetVcfToPerSample.wdl)
This utility extracts a separate single-sample VCF for each requested sample from a set of cohort VCFs, optionally dropping specified fields first. It outputs the per-sample VCFs.

Inputs:
- `Array[File] cohort_vcfs`: Cohort VCFs to extract from.
- `Array[File] cohort_vcf_idxs`: Indexes for `cohort_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Samples to extract.
- `String? drop_fields`: Fields to drop from each VCF before extraction.

Outputs:
- `subset_vcfs`: Per-sample VCFs.
- `subset_vcf_idxs`: Indexes for the per-sample VCFs.


### [AnnotateSvCallerSupport](../wdl/annotation_utils/AnnotateSvCallerSupport.wdl)
This utility annotates each SV in a cohort VCF with the set of raw callers that independently support it. For every sample it matches the cohort calls against that sample's per-caller VCFs (Kanpig, cuteSV, Sniffles, Delly, pbsv, Sawfish, dipcall and hapdiff) using reciprocal-overlap, size- and sequence-similarity and a breakpoint window, then merges the support back into the cohort VCF. It outputs the annotated VCF and a TSV of per-caller match counts.

Inputs:
- `File sv_vcf`: Cohort SV VCF to annotate.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[String] sample_ids`: Samples to process.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.
- `Array[File] kanpig_vcfs`: Per-sample Kanpig VCFs.
- `Array[File] kanpig_vcf_idxs`: Indexes for `kanpig_vcfs`.
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
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `Int truvari_breakpoint_window`: Breakpoint window, in bp, for matching a raw call (default `500`).
- `Float truvari_reciprocal_overlap`: Minimum reciprocal overlap for matching a raw call (default `0.0`).
- `Float truvari_sequence_similarity`: Minimum sequence similarity for matching a raw call (default `0.7`).
- `Float truvari_size_similarity`: Minimum size similarity for matching a raw call (default `0.7`).
- `Int fuzzy_match_breakpoint_window`: Breakpoint window, in bp, for fuzzy-matching a raw call to per-caller stats (default `500`).

Outputs:
- `sv_added_vcf`: Cohort VCF annotated with raw-caller support.
- `sv_added_vcf_idx`: Index for the annotated VCF.
- `sv_match_counts_tsv`: TSV of per-caller match counts.


## Tools

### [TransferMethylationTags](../wdl/tools/TransferMethylationTags.wdl)
This tool transfers methylation base-modification tags (MM/ML) from unaligned BAMs onto an aligned BAM. It extracts the tags per read, then per contig re-attaches them to the aligned reads and sorts, merging the result into a single tagged BAM. It outputs the methylation-tagged BAM and a TSV of the transferred tags.

Inputs:
- `Array[String] unaligned_bam_paths`: Paths to the unaligned BAMs carrying the methylation tags.
- `File aligned_bam`: Aligned BAM to receive the tags.
- `File aligned_bai`: Index for `aligned_bam`.
- `Array[String] contigs`: Contigs to process.
- `Boolean gcs_paths`: Whether `unaligned_bam_paths` are GCS paths (default `false`).
- `String mm_tag`: Base-modification tag name (default `MM`).
- `String ml_tag`: Modification-likelihood tag name (default `ML`).

Outputs:
- `methylation_tagged_bam`: Aligned BAM with methylation tags transferred.
- `methylation_tagged_bai`: Index for the tagged BAM.
- `methylation_tags`: TSV of the transferred methylation tags.
