# Workflows

## Annotations

### [AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)
This workflow leverages [AnnotateVcf](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/AnnotateVcf.wdl) from the GATK-SV pipeline in order to annotate internal allele frequencies based on sample sexes and ancestries. It runs on all variants in the input VCF, including SVs.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `File sample_pop_assignments`: Two column file containing sample IDs in the first column and ancestry labels in the second column.
- `File ped_file`: Six column file containing the cohort pedigree, with specifications described in [this article](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).
- `File lps_tsv`: LPS TSV file output by the [`trgt-lps` tool](https://github.com/PacificBiosciences/trgt-lps), which contains information regarding the longest polymer sequence within each sample for each loci genotyped by TRGT.
- `Array[String]? strip_info_fields`: Comma-separated list of `INFO` fields to remove from the input VCF prior to annotation.
- `Int records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File par_bed`: From [references](references.md).

Outputs:
- `af_annotated_vcf`: Annotated VCF.
- `af_annotated_vcf_idx`: Index for annotated VCF.


### [AnnotateAgeMetrics](../wdl/annotation/AnnotateAgeMetrics.wdl)
This workflow computes the age distribution of carriers for every variant in the input VCF. For each sample it derives an age from a date-of-birth table relative to a fixed reference date, then tabulates the number of heterozygous and homozygous carriers of each allele that fall into a set of user-defined age bins, along with overflow `smaller` and `larger` bins for ages outside the configured range. It emits a TSV of these per-allele age-bin counts.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `File age_data`: CSV file with `person_id` and `date_of_birth` columns, used to derive each sample's age.
- `Array[Int] age_bins`: Age-bin edges, in years, into which carrier ages are binned.
- `String reference_date`: Reference date (`YYYY-MM-DD`) against which each sample's age is computed.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.

Outputs:
- `annotations_tsv_age`: TSV of per-allele carrier counts across the age bins.


### [AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl)
This workflow ingests a callset VCF and two truth VCFs - one of SNVs & indels and one of SVs - and finds matching variants across them, annotating each matched callset variant with the truth callset's AC/AF/AN and genotype-count fields. This enables benchmarking annotations against an existing cohort (e.g. gnomAD) and surfacing variants that are outliers relative to it.

The workflow undergoes multiple rounds of variant matching in order to determine matched pairs:
1. Exact match across CHROM, POS, REF and ALT.
2. Truvari match with overlap percentages of 90%, 70% and 50%.
3. Matching based on `bedtools closest`, finetuned for SVs. Here the callset and truth variants are split by type and converted to a symbolic representation, after which separate `bedtools closest` passes are run - one tuned for deletions and duplications via reciprocal positional overlap, and one tuned for insertions via breakpoint proximity - so that each callset variant is paired with the nearest same-type truth variant above the per-callset minimum SV-length thresholds.

> **Note:** When converting to symbolic representation, only canonical DUPs (allele_type = `DUP` exactly) are treated as DUP; other DUP subtypes (e.g., `dup_interspersed`, `inv_dup`) are treated as insertions.

Inputs:
- `File vcf`: Callset VCF being annotated.
- `File vcf_idx`: Index for `vcf`.
- `File truth_snv_indel_vcf`: Truth VCF containing SNVs & indels to match against.
- `File truth_snv_indel_vcf_idx`: Index for `truth_snv_indel_vcf`.
- `File truth_sv_vcf`: Truth VCF containing SVs to match against.
- `File truth_sv_vcf_idx`: Index for `truth_sv_vcf`.
- `Array[String] contigs`: Contigs to evaluate.
- `Int min_sv_length_truvari_vcf`: Minimum length for a callset variant to enter the Truvari matching round.
- `Int min_sv_length_truvari_truth_vcf`: Minimum length for a truth variant to enter the Truvari matching round.
- `Int min_sv_length_bedtools_closest_vcf`: Minimum length for a callset variant to enter the `bedtools closest` matching round.
- `Int min_sv_length_bedtools_closest_truth_vcf`: Minimum length for a truth variant to enter the `bedtools closest` matching round.
- `Int? shard_bin_size_exact_match`: If set, shards the exact-match round into contig regions each containing roughly this many combined callset + truth records, run in parallel.
- `Boolean do_exact`: Whether to run the exact matching round (default `true`).
- `Boolean do_truvari`: Whether to run the Truvari matching round (default `true`).
- `Boolean do_bedtools_closest`: Whether to run the `bedtools closest` matching round (default `true`).
- `String type_field_vcf`: INFO field in the callset VCF giving each variant's allele type (default `allele_type`).
- `String length_field_vcf`: INFO field in the callset VCF giving each variant's allele length (default `allele_length`).
- `String source_tag_truth_snv_indel_vcf`: Label used to tag matches against the SNV & indel truth VCF (default `SNV_indel`).
- `String source_tag_truth_sv_vcf`: Label used to tag matches against the SV truth VCF (default `SV`).
- `String? args_string_vcf`: `bcftools view` arguments used to pre-subset the callset VCF.
- `String? args_string_truth_snv_indel_vcf`: `bcftools view` arguments used to pre-subset the SNV & indel truth VCF.
- `String? args_string_truth_sv_vcf`: `bcftools view` arguments used to pre-subset the SV truth VCF.
- `String? rename_id_string_vcf`: Expression used to rename variant IDs in the callset VCF prior to matching.
- `String? rename_id_string_truth_snv_indel_vcf`: Expression used to rename variant IDs in the SNV & indel truth VCF prior to matching.
- `String? rename_id_string_truth_sv_vcf`: Expression used to rename variant IDs in the SV truth VCF prior to matching.
- `Boolean? rename_id_strip_chr_vcf`: Whether to strip the `chr` prefix when renaming callset variant IDs.
- `Boolean? rename_id_strip_chr_truth_snv_indel_vcf`: Whether to strip the `chr` prefix when renaming SNV & indel truth variant IDs.
- `Boolean? rename_id_strip_chr_truth_sv_vcf`: Whether to strip the `chr` prefix when renaming SV truth variant IDs.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `annotations_tsv_benchmark`: TSV mapping callset variants to their matched truth variants, match type, and the truth callset's AC/AF/AN and genotype-count fields.
- `annotations_header_benchmark`: Header listing the extra annotation columns present in `annotations_tsv_benchmark`.


### [AnnotateDbSNP](../wdl/annotation/AnnotateDbSNP.wdl)
This workflow annotates each variant in the input VCF with its dbSNP reference SNP identifier (rsID). It matches variants against a per-contig dbSNP VCF on CHROM, POS, REF and ALT, emitting a TSV mapping each matched variant to its `dbSNP_ID`.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? shard_bin_size`: If set, shards each contig into regions each containing roughly this many combined VCF + dbSNP records, run in parallel.
- `File dbsnp_vcf`: From [references](references.md).
- `File dbsnp_vcf_idx`: From [references](references.md).

Outputs:
- `annotations_tsv_dbsnp`: TSV mapping variants to their dbSNP identifiers.


### [AnnotateDbVaR](../wdl/annotation/AnnotateDbVaR.wdl)
This workflow annotates structural variants in the input VCF with matching records from dbVar. It restricts to variants at or above a minimum length, converts them to a symbolic representation, and matches deletions, duplications and insertions separately against a per-contig dbVar VCF using type-specific size-similarity, reciprocal-overlap and breakpoint-window thresholds. It emits a TSV linking matched variants to their dbVar records.

> **Note:** When converting to symbolic representation, only canonical DUPs (allele_type = `DUP` exactly) are treated as DUP; other DUP subtypes (e.g., `dup_interspersed`, `inv_dup`) are treated as insertions.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int min_length`: Minimum variant length to consider for matching.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching deletions (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for matching deletions (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for matching deletions (default `0.7`).
- `Int dup_breakpoint_window`: Breakpoint window, in bp, for matching duplications (default `500`).
- `Float dup_reciprocal_overlap`: Minimum reciprocal overlap for matching duplications (default `0.7`).
- `Float dup_size_similarity`: Minimum size similarity for matching duplications (default `0.7`).
- `Int ins_breakpoint_window`: Breakpoint window, in bp, for matching insertions (default `200`).
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for matching insertions (default `0.0`).
- `Float ins_size_similarity`: Minimum size similarity for matching insertions (default `0.5`).
- `File dbvar_vcf`: From [references](references.md).
- `File dbvar_vcf_idx`: From [references](references.md).

Outputs:
- `annotations_tsv_dbvar`: TSV mapping variants to their matched dbVar records.


### [AnnotateGnomADSTR](../wdl/annotation/AnnotateGnomADSTR.wdl)
This workflow annotates tandem-repeat variants in the input VCF with overlapping loci from the gnomAD V4 tandem-repeat catalog. It subsets to tandem-repeat calls and matches each against the catalog using a minimum reciprocal-overlap threshold, emitting a TSV linking calls to their gnomAD TR locus.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File gnomad_tr_json`: From [references](references.md).
- `Float trv_reciprocal_overlap`: Minimum reciprocal overlap between a call and a catalog locus to be matched (default `0.7`).

Outputs:
- `annotations_tsv_gnomad_str`: TSV mapping tandem-repeat calls to their gnomAD TR loci.


### [AnnotateGQMetrics](../wdl/annotation/AnnotateGQMetrics.wdl)
This workflow computes binned distributions of genotype-quality metrics across the carriers of each variant. For every configured FORMAT field it counts the genotypes whose value falls into each bin - optionally restricted to a variant filter and respecting whether larger or smaller values of that field are better - and can additionally bin allele-balance values. It emits a per-variant TSV of these distribution counts.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Array[String] gq_fields`: FORMAT fields whose values are binned, one per field to annotate.
- `Array[Array[Int]] gq_bins`: Bin edges for each field in `gq_fields`.
- `Array[String] gq_variant_filters`: Per-field expression restricting which variants the field is binned over.
- `Array[Boolean] gq_larger_field`: Per-field flag indicating whether larger values of the field are better.
- `Boolean ab_annotation`: Whether to additionally compute allele-balance distributions.
- `Array[Float] ab_bins`: Bin edges for allele-balance values.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.

Outputs:
- `annotations_tsv_gq`: TSV of per-variant genotype-quality (and optional allele-balance) distribution counts.


### [AnnotateIndelTRs](../wdl/annotation/AnnotateIndelTRs.wdl)
This workflow flags short insertions and deletions that represent tandem repeats. Using the str-analysis `filter_vcf_to_tandem_repeats` tool, it inspects each indel's sequence and marks it as a tandem repeat when it meets a minimum total repeat length, minimum number of repeats and minimum repeat-unit length, emitting a TSV of the flagged variants.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCF before tandem-repeat filtering (default `-i 'INFO/allele_type!="trv" && INFO/TR_ENVELOPED!=1'`).
- `Int min_tandem_repeat_length`: Minimum total tandem-repeat length for an indel to be flagged (default `9`).
- `Int min_repeats`: Minimum number of repeats for an indel to be flagged (default `3`).
- `Int min_repeat_unit_length`: Minimum repeat-unit length for an indel to be flagged (default `1`).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `annotations_tsv_trs`: TSV of indels flagged as tandem repeats.


### [AnnotateInSilicoPredictors](../wdl/annotation/AnnotateInSilicoPredictors.wdl)
This workflow annotates SNVs and indels with precomputed in-silico predictor scores - CADD, Pangolin, PhyloP, REVEL and SpliceAI - drawn from the gnomAD V4 Hail Tables. It shards the VCF and uses a Hail-based script to look up each variant's scores, emitting a TSV of per-variant predictions.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String annotate_in_silico_predictors_script`: Path to the Hail script that performs the lookups (defaults to the `lr-annotation` repository copy).
- `String genome_build`: Genome build to annotate against (default `GRCh38`).
- `String cadd_ht`: From [references](references.md).
- `String pangolin_ht`: From [references](references.md).
- `String phylop_ht`: From [references](references.md).
- `String revel_ht`: From [references](references.md).
- `String spliceai_ht`: From [references](references.md).

Outputs:
- `annotations_tsv_insilico`: TSV of per-variant in-silico predictor scores.


### [AnnotateL1MEAID](../wdl/annotation/AnnotateL1MEAID.wdl)
This workflow first runs _RepeatMasker_ on the insertions in an input VCF. It then uses its output to run [L1ME-AID](https://github.com/Markloftus/L1ME-AID) and [INTACT_MEI](https://github.com/xzhuo/INTACT_MEI) in order to identify, annotate and filter mobile element insertion (MEI) calls. It restricts to insertions at or above a minimum length and emits a TSV of the resulting MEI annotations.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int min_length`: Minimum insertion length to consider for MEI annotation.

Outputs:
- `annotations_tsv_l1meaid`: TSV of L1ME-AID and INTACT_MEI MEI annotations.


### [AnnotateMEDs](../wdl/annotation/AnnotateMEDs.wdl)
This workflow annotates mobile element deletions (MEDs) by intersecting the deletions in the input VCF against a catalog of known mobile-element loci. Deletions are extracted to BED form and matched to the catalog using size-similarity, reciprocal-overlap, breakpoint-window and sequence-similarity thresholds, producing a TSV of the deletions identified as MEDs.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching (default `500`).
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for a deletion to match a catalog locus (default `0.0`).
- `Float del_sequence_similarity`: Minimum sequence similarity for a deletion to match a catalog locus (default `0.7`).
- `Float del_size_similarity`: Minimum size similarity for a deletion to match a catalog locus (default `0.7`).
- `File mei_catalog`: From [references](references.md).

Outputs:
- `annotations_tsv_meds`: TSV of deletions identified as mobile element deletions.


### [AnnotateMEIs](../wdl/annotation/AnnotateMEIs.wdl)
This workflow consolidates the mobile element insertion calls produced by the [AnnotateL1MEAID](#annotatel1meaid), [AnnotatePALMER](#annotatepalmer) and [AnnotateSVAN](#annotatesvan) workflows into a single harmonized set. It reconciles the three per-tool annotation TSVs, using the SVAN annotation header for typing, to produce a final TSV of MEI calls.

Inputs:
- `File annotations_tsv_l1meaid`: MEI annotation TSV output by `AnnotateL1MEAID`.
- `File annotations_tsv_palmer`: MEI annotation TSV output by `AnnotatePALMER`.
- `File annotations_tsv_svan`: MEI annotation TSV output by `AnnotateSVAN`.
- `File annotations_header_svan`: Annotation header output by `AnnotateSVAN`, used to type the consolidated fields.
- `Array[String] contigs`: Contigs to annotate.

Outputs:
- `annotations_tsv_meis`: Consolidated TSV of mobile element insertion calls.


### [AnnotatePALMER](../wdl/annotation/AnnotatePALMER.wdl)
This workflow leverages [PALMER](https://github.com/WeichenZhou/PALMER) in order to annotate MEI calls for a cohort in a given cohort VCF. It retains the genotypes present in the VCF, simply adding an INFO field `ME_TYPE` to insertions whose characteristics match those of the PALMER calls. Matching is performed per MEI type using type-specific reciprocal-overlap, size-similarity, sequence-similarity, breakpoint-window and minimum-shared-sample thresholds.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File PALMER_vcf`: VCF of PALMER MEI calls to match against.
- `File PALMER_vcf_idx`: Index for `PALMER_vcf`.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Array[String] mei_types`: MEI types to run on - must be a subset of [`ALU`, `SVA`, `LINE` or `HERVK`].
- `Int min_length`: Minimum insertion length to consider for annotation.
- `File rm_out`: _RepeatMasker_ output for the input VCF's insertions.
- `Int rm_buffer`: Padding, in bp, applied around RepeatMasker annotations when matching.
- `Int ins_breakpoint_window_{ALU,SVA,LINE,HERVK}`: Per-type breakpoint window, in bp, for matching (default `500`).
- `Float ins_reciprocal_overlap_{ALU,SVA,LINE,HERVK}`: Per-type minimum reciprocal overlap for matching (default `0.9`).
- `Float ins_sequence_similarity_{ALU,SVA,LINE,HERVK}`: Per-type minimum sequence similarity for matching (default `0.9`).
- `Float ins_size_similarity_{ALU,SVA,LINE,HERVK}`: Per-type minimum size similarity for matching (default `0.9`).
- `Int ins_min_shared_samples_{ALU,SVA,LINE,HERVK}`: Per-type minimum number of shared samples for matching (default `0`).
- `File ref_fai`: From [references](references.md).

Outputs:
- `annotations_tsv_palmer`: TSV of insertions annotated with their PALMER `ME_TYPE`.


### [AnnotateRegion](../wdl/annotation/AnnotateRegion.wdl)
This workflow annotates each variant with the genomic region class it falls within - simple repeat (`SR`), segmental duplication (`SD`), RepeatMasker region (`RM`) or unique sequence (`US`) - by intersecting it against the corresponding BED panels. It emits a TSV of per-variant `REGION` assignments.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File simple_repeats_bed`: From [references](references.md).
- `File seg_dup_bed`: From [references](references.md).
- `File repeat_masked_bed`: From [references](references.md).

Outputs:
- `annotations_tsv_region`: TSV of per-variant genomic-region assignments.


### [AnnotateSQMetrics](../wdl/annotation/AnnotateSQMetrics.wdl)
This workflow recomputes site-level quality metrics for each variant directly from its genotype-level data. It clears any stale allele-specific INFO fields and recalculates Hardy-Weinberg equilibrium, the inbreeding coefficient, the maximum p(allele balance), and the allele-specific quality approximation, quality-by-depth and variant depth from the per-sample DP, PL and AD fields, emitting a per-variant TSV.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.

Outputs:
- `annotations_tsv_sq`: TSV of recomputed site-level quality metrics.


### [AnnotateSVAN](../wdl/annotation/AnnotateSVAN.wdl)
This workflow leverages [SVAN](https://github.com/REPBIO-LAB/SVAN) in order to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It processes insertions and deletions separately, first running Tandem Repeat Finder (TRF) on the inserted or deleted sequence of each SV in the input VCF and then running SVAN over the result, before extracting and aligning the annotations into a single TSV. Before extraction, the `DUP_COORD` field produced by SVAN is reformatted: any `flank_`-prefixed relative coordinates are resolved to absolute genomic positions, preserving the original order of comma-separated values.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String type_field`: INFO field giving each variant's allele type, used to select insertions/deletions for annotation (default `allele_type`).
- `String type_ins`: Value of `type_field` identifying an insertion (default `ins`).
- `String type_del`: Value of `type_field` identifying a deletion (default `del`).
- `String length_field`: INFO field giving each variant's allele length (default `allele_length`).
- `Int min_length`: Minimum insertion/deletion length to consider for annotation (default `0`).
- `Boolean annotate_ins`: Whether to annotate insertions (default `true`).
- `Boolean annotate_del`: Whether to annotate deletions (default `true`).
- `File vntr_bed`: From [references](references.md).
- `File exons_bed`: From [references](references.md).
- `File repeats_bed`: From [references](references.md).
- `File ref_fa`: From [references](references.md).
- `Array[File] ref_fa_indices`: From [references](references.md).
- `File mei_fa`: From [references](references.md).
- `Array[File] mei_fa_indices`: From [references](references.md).

Outputs:
- `annotations_tsv_svan`: TSV of SVAN annotations.
- `annotations_header_svan`: Header lines describing the SVAN annotation fields.


### [AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl)
This workflow leverages [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) in order to annotate predicted functional effects for SVs. It conditionally only runs SVs through this workflow, ignoring all SNVs and indels, converting each SV to a symbolic representation before annotating it against coding and noncoding panels and extracting the resulting `PREDICTED_` annotations into a TSV.

> **Note:** When converting to symbolic representation, all DUP allele types (including `dup_interspersed`, `inv_dup`, `complex_dup`, etc.) are treated as DUP.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int min_length`: Minimum length for a variant to be treated as an SV and annotated.
- `File coding_gtf`: From [references](references.md).
- `File noncoding_bed`: From [references](references.md).

Outputs:
- `annotations_tsv_svannotate`: TSV of SVAnnotate functional-effect annotations.
- `annotations_header_svannotate`: Header lines describing the SVAnnotate annotation fields.


### [AnnotateTruvariRemap](../wdl/annotation/AnnotateTruvariRemap.wdl)
This tool remaps insertion sequences with minimap2 (via Truvari) in order to flag insertions whose inserted sequence aligns elsewhere in the reference. Each insertion above a minimum length is realigned per contig and assessed against alignment-score and coverage thresholds, emitting a TSV of the remap results.

Inputs:
- `File vcf`: VCF whose insertions are remapped.
- `File vcf_idx`: Index for VCF.
- `Array[String] contigs`: Contigs to process.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String type_field`: INFO field giving each variant's allele type, used to select insertions to remap (default `allele_type`).
- `String type_ins`: Value of `type_field` identifying an insertion (default `ins`).
- `Int min_length`: Minimum insertion length to remap.
- `Int max_length`: Maximum insertion length to remap.
- `Int mm2_threshold`: Minimum minimap2 alignment score to flag an insertion.
- `Float cov_threshold`: Minimum alignment coverage to flag an insertion.
- `File ref_fa`: From [references](references.md).
- `Array[File] ref_bwa_indices`: BWA indices for `ref_fa`, from [references](references.md).

Outputs:
- `annotations_tsv_remap`: TSV of insertion remap results.


### [AnnotateVEPHail](../wdl/annotation/AnnotateVEPHail.wdl)
This workflow leverages [the Ensembl Variant Effect Predictor (VEP)](https://useast.ensembl.org/info/docs/tools/vep/index.html) in order to annotate predicted functional effects based on site-level information. It strips genotypes, scatters the VCF into shards, optionally normalizes and splits multiallelics around the VEP call, and uses Hail in order to run this annotation process in a more efficient and scalable manner before concatenating the per-shard annotations into a single TSV.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `String? subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCF before annotation.
- `String split_vcf_hail_script`: Path to the Hail script used to scatter the VCF (defaults to the `lr-annotation` repository copy).
- `String vep_annotate_hail_python_script`: Path to the Hail script used to run VEP (defaults to the `lr-annotation` repository copy).
- `String genome_build`: Genome build to annotate against (default `GRCh38`).
- `String vep_json_schema`: Hail type schema describing the structure of VEP's JSON output.
- `String normalize_check_ref`: `bcftools norm` `--check-ref` mode used when normalizing (default `w`).
- `Boolean normalize_vcf`: Whether to normalize and split multiallelics around the VEP call (default `false`).
- `Boolean localize_vcf`: Whether to localize the VCF before annotation (default `true`).
- `Boolean has_index`: Whether the input VCF is indexed (default `true`).
- `Boolean get_chromosome_sizes`: Whether to compute chromosome sizes (default `false`).
- `Boolean split_by_chromosome`: Whether to scatter the VCF by chromosome (default `false`).
- `Boolean split_into_shards`: Whether to scatter the VCF into fixed-size shards (default `false`).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).
- `File ref_fa_gz`: bgzipped `ref_fa`, from [references](references.md).
- `File ref_fai_gz`: Index for `ref_fa_gz`, from [references](references.md).
- `File ref_vep_cache`: From [references](references.md).

Outputs:
- `annotations_tsv_vep`: TSV of VEP functional-effect annotations.


### [AnnotateVRS](../wdl/annotation/AnnotateVRS.wdl)
This workflow annotates each variant with its GA4GH Variant Representation Specification (VRS) attributes using a seqrepo sequence repository. It runs `vrs-annotate` per contig to add the VRS INFO fields, then extracts them into an annotation TSV of five locating columns (CHROM, POS, REF, ALT, ID) followed by a column for each VRS field.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File seqrepo_tar`: From [references](references.md).

Outputs:
- `annotations_tsv_vrs`: TSV of per-variant VRS attributes (`VRS_Allele_IDs`, `VRS_Error`, `VRS_Starts`, `VRS_Ends`, `VRS_States`, `VRS_Lengths`, `VRS_RepeatSubunitLengths`).


## Annotation Utilities

### [AddEndTRs](../wdl/annotation_utils/AddEndTRs.wdl)
This utility adds an `END` INFO tag to the tandem-repeat records of a VCF, computed per contig, so that downstream tools correctly interpret the span of each TR call. It outputs the updated VCF.

Inputs:
- `File vcf`: VCF to update.
- `File vcf_idx`: Index for VCF to update.
- `Array[String] contigs`: Contigs to process within the input VCF.

Outputs:
- `vcf_with_end`: VCF with `END` tags added to tandem-repeat records.
- `vcf_with_end_idx`: Index for the updated VCF.


### [AnnotateAlleleType](../wdl/annotation_utils/AnnotateAlleleType.wdl)
This utility sets the `allele_type` INFO field on variants in a VCF using three annotation TSVs - one for mobile element deletions, one for mobile element insertions and one for duplications - applying each in turn. Each annotation source can have its values transformed via an optional prefix, suffix and lowercasing. It outputs the annotated VCF.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File med_tsv`: TSV of mobile element deletion allele types.
- `File mei_tsv`: TSV of mobile element insertion allele types.
- `File dup_tsv`: TSV of duplication allele types.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `String? med_prefix`: Prefix prepended to mobile element deletion allele-type values.
- `String? med_suffix`: Suffix appended to mobile element deletion allele-type values.
- `Boolean? med_lowercase`: Whether to lowercase mobile element deletion allele-type values.
- `String? mei_prefix`: Prefix prepended to mobile element insertion allele-type values.
- `String? mei_suffix`: Suffix appended to mobile element insertion allele-type values.
- `Boolean? mei_lowercase`: Whether to lowercase mobile element insertion allele-type values.
- `String? dup_prefix`: Prefix prepended to duplication allele-type values.
- `String? dup_suffix`: Suffix appended to duplication allele-type values.
- `Boolean? dup_lowercase`: Whether to lowercase duplication allele-type values.

Outputs:
- `allele_type_annotated_vcf`: VCF annotated with `allele_type`.
- `allele_type_annotated_vcf_idx`: Index for the annotated VCF.


### [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)
This utility applies a list of annotation TSVs to a VCF as new INFO fields, adding the specified field names, descriptions, types and numbers for each TSV in turn. Each annotation source can optionally have its TSV sorted, the VCF pre-subset and its TSV rows filtered beforehand. It outputs the annotated VCF.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[File] annotations_tsvs`: Annotation TSVs to apply, each as a set of INFO fields.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation. When set, each contig VCF is sharded by record count, annotated in parallel and concatenated.
- `Array[String]? strip_info_fields`: INFO fields to remove from the input VCF prior to annotation.
- `Array[Boolean] sort_tsvs`: Per-TSV flag indicating whether to sort the TSV before annotation (default empty).
- `Array[String] subset_vcf_strings`: Per-TSV `bcftools view` arguments used to pre-subset the VCF (default empty).
- `Array[String] awk_tsv_conditions`: Per-TSV `awk` condition used to filter the TSV rows applied (default empty).
- `Array[Array[String]] info_names`: INFO field names added by each annotation TSV.
- `Array[Array[String]] info_descriptions`: INFO field header descriptions for each annotation TSV.
- `Array[Array[String]] info_types`: INFO field types for each annotation TSV.
- `Array[Array[String]] info_numbers`: INFO field `Number` values for each annotation TSV.

Outputs:
- `annotated_vcf`: Annotated VCF.
- `annotated_vcf_idx`: Index for the annotated VCF.


### [AnnotateVcfCleared](../wdl/annotation_utils/AnnotateVcfCleared.wdl)
This utility is a variant of [AnnotateVcf](#annotatevcf) that, before applying the annotation TSVs, clears existing annotations and optionally swaps in records from an untrimmed VCF in order to restore full REF/ALT alleles. It then adds the specified INFO fields and outputs the annotated VCF.

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

Outputs:
- `annotated_vcf`: Annotated VCF.
- `annotated_vcf_idx`: Index for the annotated VCF.


### [CombineTRs](../wdl/annotation_utils/CombineTRs.wdl)
This utility combines tandem-repeat VCFs from multiple callers for a single sample into one VCF. It checks sample consistency, sets missing filters to pass, tags each caller's calls, assigns TR identifiers, deduplicates overlapping variants and priority-merges the callers per contig. It outputs the combined TR VCF.

Inputs:
- `Array[File] tr_vcfs`: Tandem-repeat VCFs to combine, one per caller.
- `Array[File] tr_vcf_idxs`: Indexes for `tr_vcfs`.
- `Array[String] tr_callers`: Caller name for each VCF in `tr_vcfs`, used for tagging and merge priority.
- `Array[String] contigs`: Contigs to process.
- `String sample_id`: Sample ID expected across all input VCFs.

Outputs:
- `trgt_combined_vcf`: Combined single-sample tandem-repeat VCF.
- `trgt_combined_vcf_idx`: Index for the combined VCF.


### [ConcatenateMosDepth](../wdl/annotation_utils/ConcatenateMosDepth.wdl)
This utility concatenates a sample's per-contig [MosDepth](#mosdepth) per-base coverage BED files into a single indexed BED. It outputs the combined per-base coverage BED and its index.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-contig mosdepth per-base coverage BED files to concatenate.

Outputs:
- `mosdepth_per_base_combined`: Combined per-base coverage BED.
- `mosdepth_per_base_combined_idx`: Index for the combined BED.


### [CountAnnotations](../wdl/annotation_utils/CountAnnotations.wdl)
This utility tallies annotation values across one or more VCFs to produce summary count tables, size-binned by allele class (SNV/DEL/INS/DUP/TRV). It always counts at the site level and can optionally count per sample, per allele, per functional gene consequence (from VEP/SVAnnotate `PREDICTED_*` fields), as raw per-variant value lists, and — when `create_plotting` is enabled — produce a separate set of AF-binned, region-aware Parquet tables for plotting (including a de novo transmission breakdown when a PED file is supplied and trios are found).

Inputs:
- `Array[File] vcfs`: VCFs whose annotations are counted.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[Int] length_bins_summary`: Size-bin edges used for the summary count tables (default `[0, 1, 50, 500]`).
- `Array[Int] length_bins_plotting`: Size-bin edges used for the plotting tables (default `[0, 1, 50, 100, 500, 5000, 50000]`).
- `Array[Float] af_bins_plotting`: Allele-frequency bin edges used for the plotting tables (default `[0.0, 0.01, 0.05, 0.1, 0.5]`).
- `Boolean create_per_sample`: Whether to additionally produce per-sample counts (default `false`).
- `Boolean create_per_allele`: Whether to additionally produce per-allele counts (default `false`).
- `Boolean create_list`: Whether to additionally produce raw per-variant value-list tables (default `false`).
- `Boolean create_functional`: Whether to additionally produce per-gene functional counts (default `false`).
- `Boolean create_plotting`: Whether to additionally produce AF-binned Parquet tables for plotting (default `false`).
- `Boolean use_ssd`: Whether to use SSD-backed local disks (default `false`).
- `Boolean split_by_region`: Whether to split each VCF by genomic region before counting (default `false`).
- `String subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCFs (default empty).
- `Int max_length`: Maximum variant length to count, or `-1` for no maximum (default `-1`).
- `Int min_length`: Minimum variant length to count, or `-1` for no minimum (default `-1`).
- `Int? records_per_shard`: Number of variants to keep within a single shard.
- `File? ped`: PED file used to identify trios for the de novo transmission breakdown (only used when `create_plotting` is enabled).

Outputs:
- `summary_sites_tsv`: Site-level annotation counts.
- `summary_samples_tsv`: Per-sample counts (when `create_per_sample`).
- `summary_alleles_tsv`: Per-allele counts (when `create_per_allele`).
- `summary_list_tsv`: Raw per-variant value lists (when `create_list`).
- `summary_functional_tsv`: Per-gene functional counts (when `create_functional`).
- `summary_functional_samples_tsv`: Per-gene per-sample functional counts (when `create_functional` and `create_per_sample`).
- `summary_functional_alleles_tsv`: Per-gene per-allele functional counts (when `create_functional` and `create_per_allele`).
- `plotting_sites_parquet`: Site-level AF/size-binned counts, as Parquet (when `create_plotting`).
- `plotting_samples_parquet`: Per-sample AF/size-binned counts, as Parquet (when `create_plotting`).
- `plotting_alleles_parquet`: Per-allele AF/size-binned counts, as Parquet (when `create_plotting`).
- `plotting_denovo_parquet`: Per-proband de novo transmission counts, as Parquet (when `create_plotting` and trios are found via `ped`).
- `plotting_variant_list_parquet`: Raw per-variant genotype-count list, as Parquet (when `create_plotting`).


### [CreateCohortMethylationFile](../wdl/annotation_utils/CreateCohortMethylationFile.wdl)
This utility builds cohort-level CpG methylation matrices from per-sample [MethylationProfiling](../wdl/tools/MethylationProfiling.wdl) BED outputs. For each contig, it merges every sample's combined and per-haplotype modification-score BEDs into a wide site-by-sample(/haplotype) matrix, filling `.` for sites missing in a given sample or haplotype. Samples can optionally be processed in shards (merged independently, then joined column-wise) to bound how many sample files are localized onto a single task at once.

Inputs:
- `Array[File] combined_beds`: Per-sample combined CpG pileup BEDs (`cpg_combined_bed` from MethylationProfiling).
- `Array[File] combined_bed_idxs`: Indexes for `combined_beds`.
- `Array[File] hap1_beds`: Per-sample haplotype 1 CpG pileup BEDs (`cpg_hap1_bed` from MethylationProfiling).
- `Array[File] hap1_bed_idxs`: Indexes for `hap1_beds`.
- `Array[File] hap2_beds`: Per-sample haplotype 2 CpG pileup BEDs (`cpg_hap2_bed` from MethylationProfiling).
- `Array[File] hap2_bed_idxs`: Indexes for `hap2_beds`.
- `Array[String] sample_ids`: Sample IDs, parallel to `combined_beds`/`hap1_beds`/`hap2_beds`.
- `Array[String] contigs`: Contigs to process.
- `Int? samples_per_shard`: Maximum number of samples merged per shard before shards are joined column-wise into the final matrix.

Outputs:
- `combined_methylation_beds`: Per-contig site-by-sample modification-score matrix BEDs.
- `haplotype_methylation_beds`: Per-contig site-by-haplotype modification-score matrix BEDs.


### [CreateCoverageFile](../wdl/annotation_utils/CreateCoverageFile.wdl)
This utility builds a binned coverage matrix across a cohort from per-sample mosdepth BED outputs. It tiles the genome into windows, computes the mean coverage and threshold-crossing counts within each bin for every sample, and concatenates the results into a single coverage TSV.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-sample mosdepth coverage BED files.
- `Array[File] mosdepth_bed_indices`: Indexes for `mosdepth_bed_files`.
- `Array[String] contigs`: Contigs over which to compute coverage.
- `Int window_size`: Size, in bp, of each genomic window.
- `Int bin_size`: Size, in bp, of each coverage bin within a window.
- `Array[Int] thresholds`: Coverage thresholds at which to count bins as covered.
- `File ref_fai`: From [references](references.md).

Outputs:
- `binned_coverage_tsv`: Binned coverage matrix across the cohort.


### [CreateIntervalsFile](../wdl/annotation_utils/CreateIntervalsFile.wdl)
This utility creates an interval file matching the fixed-width bins emitted by [CreateReadCountsFile](#createreadcountsfile). It writes 1-based, inclusive `contig:start-end` intervals in the requested contig order and omits each contig's trailing partial bin.

Inputs:
- `File ref_fai`: From [references](references.md).
- `Array[String] contigs`: Contigs for which to create intervals, in output order.
- `Int bin_size`: Size, in bp, of each interval. Use the same value supplied to `CreateReadCountsFile`.

Outputs:
- `intervals`: Fixed-width interval file.


### [CreateDepthFiles](../wdl/annotation_utils/CreateDepthFiles.wdl)
This utility ports GATK-SV's [MakeBincovMatrix](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/MakeBincovMatrix.wdl) and [PloidyEstimation](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/PloidyEstimation.wdl) workflows to build a cohort binned-coverage matrix and per-sample ploidy estimate from per-sample [MosDepth](#mosdepth) per-base coverage BEDs. Since mosdepth's per-base output is run-length-encoded at irregular interval widths rather than GATK-SV's fixed-width `CollectReadCounts` bins, each sample's per-base BED is first binned at `bin_size` by taking the median depth per bin (dropping any trailing partial bin), matching the binning convention used by [CreateReadCountsFile](#createreadcountsfile); because every sample is binned identically, the format-detection/shift logic in upstream `MakeBincovMatrix` (which has to distinguish raw bincov BEDs from GATK `CollectReadCounts` output) is dropped as dead code. The binned files are then run through GATK-SV's `SetBins`/`MakeBincovMatrixColumns`/`ZPaste` logic to build the bincov matrix, and through `BuildPloidyMatrix` (re-binning the bincov matrix to `ploidy_bin_size`, summing depths) and GATK-SV's `estimatePloidy.R` to estimate ploidy. GATK-SV's `estimatePloidy.R` and `estimated_CN_denoising.py` are vendored under [`scripts/helper/`](../scripts/helper/) and built into the `utils` image, so workflow has no dependency on GATK-SV docker images. Matrix outputs remain separate; `ploidy_plots` tarball contains only PNG figures from `estimatePloidy.R` and `cn_denoising_plots.pdf`. Unlike upstream `MakeBincovMatrix`, this does not support merging into a pre-existing batch's bincov matrix, since only a single one-shot cohort matrix was needed.

`estimatePloidy.R` hardcodes a 24-contig human karyotype (`chr1`..`chr22`, `chrX`, `chrY`, in that exact order) for sex assignment and per-contig ploidy expectations via positional indexing, and its "X"/"Y" exclusion checks compare against bare `X`/`Y` rather than `chr`-prefixed names (a no-op against GRCh38-style contig names, with limited practical effect here since sample-batching/PCA (`-k`) is never invoked). `mosdepth_bed_files` must therefore be restricted to exactly those 24 contigs, in that order, or ploidy estimates will be silently wrong.

Inputs:
- `Array[String] sample_ids`: Cohort sample IDs, parallel to `mosdepth_bed_files`.
- `Array[File] mosdepth_bed_files`: Per-sample combined mosdepth per-base coverage BEDs, restricted to `chr1`-`chr22`, `chrX`, `chrY` in that order (see caveat above).
- `Int bin_size`: Size, in bp, of each coverage bin in the bincov matrix (GATK-SV convention default: 1000).
- `Int ploidy_bin_size`: Size, in bp, of each bin in the ploidy matrix (GATK-SV convention default: 1000000).

Outputs:
- `binned_coverage`: Cohort binned-coverage matrix, bgzipped and tabix-indexed.
- `binned_coverage_idx`: Index for `binned_coverage`.
- `median_coverage`: Per-sample median coverage matrix.
- `estimated_cn`: Per-sample, per-chromosome estimated copy number.
- `binned_estimated_ecn`: Per-sample, per-`ploidy_bin_size`-bin estimated copy number.
- `ploidy_plots`: Tarball containing only ploidy PNG and PDF figures.


### [CreateMetadataFile](../wdl/annotation_utils/CreateMetadataFile.wdl)
This utility builds a cohort metadata file by combining a pedigree file with an ancestry-assignment file. It outputs the merged metadata file.

Inputs:
- `File ped_file`: Cohort pedigree file.
- `File ancestry_file`: Two-column file of sample IDs and ancestry labels.

Outputs:
- `metadata`: Merged cohort metadata file.


### [CreateReadCountsFile](../wdl/annotation_utils/CreateReadCountsFile.wdl)
This utility produces a binned read-counts file for a single sample from its per-contig mosdepth BED outputs, binning counts at a fixed resolution and merging across contigs. It outputs the binned read-counts file.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-contig mosdepth coverage BED files for the sample.
- `Array[File] mosdepth_bed_indices`: Indexes for `mosdepth_bed_files`.
- `Array[String] contigs`: Contigs over which to bin read counts.
- `Int bin_size`: Size, in bp, of each read-count bin.
- `String sample_id`: ID of the sample being processed.
- `File ref_dict`: From [references](references.md).

Outputs:
- `binned_read_counts`: Binned read-counts file for the sample.


### [DiagnoseSingletons](../wdl/annotation_utils/DiagnoseSingletons.wdl)
This utility counts each sample's called genotypes across variant type, allele-length range, genomic region, evidence source, and sample-level alternate-allele count. It classifies calls supported only by `hapdiff` and/or `dipcall` as assemblies, calls with any other `EV` value as alignments, and calls without either kind of evidence as other. Output columns use the format `variant_type - size_range - region - count_type - singleton_type`.

Inputs:
- `Array[File] vcfs`: VCFs whose sample calls are counted.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Int min_length`: Minimum absolute `INFO/allele_length` to count (default `50`).
- `Int? records_per_shard`: Number of variants to keep within a single shard.

Outputs:
- `singleton_counts_tsv`: Wide per-sample count table.


### [DropGenotypes](../wdl/annotation_utils/DropGenotypes.wdl)
This utility strips all genotype (sample) columns from a VCF, optionally sharding by record count for speed. It outputs the resulting sites-only VCF.

Inputs:
- `File vcf`: VCF whose genotypes are dropped.
- `File vcf_idx`: Index for VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.

Outputs:
- `dropped_vcf`: Sites-only VCF.
- `dropped_vcf_idx`: Index for the sites-only VCF.


### [ExtractSampleVcfs](../wdl/annotation_utils/ExtractSampleVcfs.wdl)
This utility extracts per-sample VCFs from a cohort VCF, splitting each sample's variants into a SNV/indel VCF and an SV VCF based on a minimum SV length. It outputs the per-sample SNV/indel and SV VCFs.

Inputs:
- `File cohort_vcf`: Cohort VCF to extract from.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[String] contigs`: Contigs to process within the cohort VCF.
- `Array[String] sample_ids`: Samples to extract.
- `Int min_sv_length`: Minimum length at which a variant is routed to the SV VCF rather than the SNV/indel VCF.

Outputs:
- `snv_indel_vcfs`: Per-sample SNV/indel VCFs.
- `snv_indel_vcf_idxs`: Indexes for the SNV/indel VCFs.
- `sv_vcfs`: Per-sample SV VCFs.
- `sv_vcf_idxs`: Indexes for the SV VCFs.


### [FlagLowCoverageRegions](../wdl/annotation_utils/FlagLowCoverageRegions.wdl)
This utility finds recurrent low-coverage regions from cohort mosdepth per-base BED files. It streams each sample independently, divides each chromosome into fixed bins anchored at position 0, and calculates each bin's base-weighted median coverage. Each sample's regular low-coverage cutoff is its median binned coverage multiplied by `median_coverage_cutoff`; bins at or below the cutoff are flagged. A cohort bin fails when its low-coverage sample proportion is at or above `sample_proportion_cutoff`.

Sample sex is read from a six-column PED (`1` for male and `2` for female). Male chrX and chrY use half the regular cutoff. Female chrY is excluded from binning, the sample median, sample histograms, and low-coverage calls. Samples missing from the PED or carrying unsupported sex codes fail explicitly. The workflow scatters once across samples, then aggregates flagged bins in one cohort task to avoid a nested scatter.

Inputs:
- `Array[File] mosdepth_bed_files`: One sorted, contiguous mosdepth per-base BED or BED.GZ per sample. Files must use four columns (`chrom`, `start`, `end`, `coverage`), cover each included chromosome from position 0, and have the same chromosome coordinate system.
- `Array[String] sample_ids`: Sample IDs corresponding by array index to `mosdepth_bed_files`.
- `File ped`: Six-column PED containing every input sample and its sex.
- `Int bin_size`: Fixed genomic bin size in bp.
- `Float median_coverage_cutoff`: Fraction of sample median coverage defining the inclusive regular low-coverage cutoff (default `0.2`; must be greater than `0` and at most `1`). For example, `0.2` gives a `12x` cutoff for a sample with median coverage `60x`.
- `Float sample_proportion_cutoff`: Inclusive minimum proportion of eligible samples with low coverage required to fail a cohort bin (default `0.5`; must be greater than `0` and at most `1`).

Outputs:
- `sample_histograms_tar`: Tarball of per-sample coverage histograms with weighted, data-driven bins. The displayed range extends through at least the 95th percentile and the `Q3 + 3 * IQR` upper fence; omitted extreme bins are counted, actual low-coverage bins are orange, and the sample median plus regular/sex-chromosome cutoffs are annotated.
- `chromosome_low_coverage_tar`: Tarball of per-chromosome plots showing number of samples flagged in each genomic bin.
- `cohort_low_coverage_tsv`: Gzipped TSV underlying the chromosome plots, with low-coverage and eligible-sample counts for each bin flagged by at least one sample. chrY eligibility includes male samples only.
- `failed_bins_bed`: Gzipped, naturally chromosome-sorted BED of cohort bins whose low-coverage sample proportion is greater than or equal to `sample_proportion_cutoff`.
- `sample_cutoffs_tsv`: TSV with `sample_id`, regular `cutoff`, and `median_coverage` for every input sample.


### [FillPhasedGenotypes](../wdl/annotation_utils/FillPhasedGenotypes.wdl)
This utility transfers phasing information from a phased VCF onto the genotypes of an unphased VCF over matching sites, optionally sharding each contig by region. It outputs the phased VCF.

Inputs:
- `File phased_vcf`: VCF providing the phasing information.
- `File phased_vcf_idx`: Index for `phased_vcf`.
- `File unphased_vcf`: VCF whose genotypes are phased.
- `File unphased_vcf_idx`: Index for `unphased_vcf`.
- `Array[String] contigs`: Contigs to process.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding each contig.

Outputs:
- `hiphase_phased_vcf`: Phased VCF.
- `hiphase_phased_vcf_idx`: Index for the phased VCF.


### [GenerateTRGTJson](../wdl/annotation_utils/GenerateTRGTJson.wdl)
This utility generates per-locus tandem repeat allele-frequency histograms, stratified by population and sex, from a multisample LPS (longest polymer sequence) table for use in the TR browser. It outputs a single combined histograms TSV.

Inputs:
- `File lps_tsv`: Multisample LPS table.
- `File metadata_tsv`: Sample metadata (population, sex) used to stratify the histograms.
- `Array[String] contigs`: Contigs to process within the LPS table.

Outputs:
- `trgt_histograms_tsv`: Combined per-locus allele-frequency histograms TSV.


### [IntegrateTRs](../wdl/annotation_utils/IntegrateTRs.wdl)
This utility integrates tandem-repeat calls into a base VCF for a cohort. It aligns samples between the base and TR VCFs, sets missing filters to pass, tags TR records with their source catalog, assigns TR identifiers and annotates the base VCF with the integrated TR calls. It outputs the TR-annotated VCF.

Inputs:
- `File vcf`: Base VCF to integrate TRs into.
- `File vcf_idx`: Index for the base VCF.
- `File tr_vcf`: Tandem-repeat VCF to integrate.
- `File tr_vcf_idx`: Index for the TR VCF.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Samples shared between the base and TR VCFs.
- `Array[File] tr_catalogs`: Catalogs from which the TR calls were derived.
- `Array[String] tr_catalog_ids`: Identifier for each catalog in `tr_catalogs`.

Outputs:
- `tr_annotated_vcf`: Base VCF annotated with integrated TR calls.
- `tr_annotated_vcf_idx`: Index for the annotated VCF.


### [IntegrateVcfs](../wdl/annotation_utils/IntegrateVcfs.wdl)
This utility integrates a SNV/indel VCF and an SV VCF into a single cohort VCF. Each input is normalized, harmonized to a common sample set and tagged with a source label and a size-based flag, after which the two are merged and the combined variants are renamed and filtered - for example to flag large SNVs/indels and small SVs. Sample IDs can optionally be swapped first. It outputs the integrated VCF.

Inputs:
- `File snv_indel_vcf`: SNV/indel VCF to integrate.
- `File snv_indel_vcf_idx`: Index for `snv_indel_vcf`.
- `File sv_vcf`: SV VCF to integrate.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Samples shared between the two VCFs.
- `String snv_indel_vcf_source_tag`: `SOURCE` value applied to the SNV/indel calls.
- `String snv_indel_vcf_size_flag`: Filter flag applied to out-of-range SNV/indel calls.
- `String snv_indel_vcf_size_flag_description`: Header description for `snv_indel_vcf_size_flag`.
- `String sv_vcf_source_tag`: `SOURCE` value applied to the SV calls.
- `String sv_vcf_size_flag`: Filter flag applied to out-of-range SV calls.
- `String sv_vcf_size_flag_description`: Header description for `sv_vcf_size_flag`.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.
- `Int min_sv_length`: Length boundary separating SNVs/indels from SVs (default `50`).
- `File? swap_samples_snv_indel`: Sample-ID swap map applied to the SNV/indel VCF.
- `File? swap_samples_sv`: Sample-ID swap map applied to the SV VCF.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `integrated_vcf`: Integrated cohort VCF.
- `integrated_vcf_idx`: Index for the integrated VCF.


### [MergeBackbonePhased](../wdl/annotation_utils/MergeBackbonePhased.wdl)
This utility merges a backbone-phased VCF with its no-TRGT counterpart (the same backbone-phasing run without TRGT calls included): for each still-unphased heterozygous genotype in `backbone_phased_vcf`, if a matching variant exists in `backbone_phased_notrgt_vcf` with a phased genotype, that phased `GT` (and `PS`) is pulled into the output. Region sharding is optional. It outputs the merged VCF and a per-sample TSV of heterozygous/unphased/pulled genotype counts.

Inputs:
- `File backbone_phased_vcf`: Backbone-phased VCF whose remaining unphased het genotypes are filled.
- `File backbone_phased_vcf_idx`: Index for `backbone_phased_vcf`.
- `File backbone_phased_notrgt_vcf`: Backbone-phased VCF (without TRGT calls) providing phased genotypes to pull from.
- `File backbone_phased_notrgt_vcf_idx`: Index for `backbone_phased_notrgt_vcf`.
- `String contig`: Contig to process.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the contig.

Outputs:
- `backbone_merged_vcf`: Merged VCF with phased genotypes pulled in where available.
- `backbone_merged_vcf_idx`: Index for the merged VCF.
- `backbone_merged_tsv`: Per-sample TSV of heterozygous, unphased, and post-pull unphased genotype counts.


### [PlotPhasingResults](../wdl/annotation_utils/PlotPhasingResults.wdl)
This utility evaluates backbone-phasing accuracy by comparing backbone-phased VCFs against base (truth) VCFs. It assigns samples to their base VCFs, compares phased genotypes per contig and aggregates the results into tables broken down by variants outside tandem repeats, TR-enveloped variants and TR variants. It outputs these summary tables plus per-VCF status tables.

Inputs:
- `Array[File] backbone_phased_vcfs`: Backbone-phased VCFs to evaluate.
- `Array[File] backbone_phased_vcf_idxs`: Indexes for `backbone_phased_vcfs`.
- `Array[File] base_vcfs`: Base (truth) VCFs to compare against.
- `Array[File] base_vcf_idxs`: Indexes for `base_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Int max_variants`: Maximum number of variants to evaluate, or `-1` for no limit (default `-1`).
- `Array[String]? subset_samples`: Samples to restrict the evaluation to.

Outputs:
- `outside_tr_table`: Phasing-accuracy table for variants outside tandem repeats.
- `tr_enveloped_table`: Phasing-accuracy table for TR-enveloped variants.
- `trv_table`: Phasing-accuracy table for tandem-repeat variants.
- `missing_samples`: Samples with no matching base VCF.
- `vcf_tables`: Per-VCF variant-status tables.


### [PostProcess](../wdl/annotation_utils/PostProcess.wdl)
This utility bundles every genotype-update and post-processing step applied to a near-final callset into one workflow, with a required `run_` Boolean guarding each step so that the input VCF is left untouched when all are set to `false`. The per-record steps are applied in a single pass over the VCF: each variant is first matched against `transfer_vcf` and has its genotypes transferred (when `run_transfer_genotypes` is set) using its unmodified properties, after which the remaining steps - unphasing, ploidy normalization, TR-ID decrementing, MEI pruning, homopolymer flagging, singleton filtering and same-coordinate sorting - run in order. Some steps require an accompanying field - `run_transfer_genotypes` needs `transfer_vcf`, `run_unphase_samples` needs `unphase_samples`, and `run_normalize_ploidy` needs `ped`. The per-record pass can optionally be region-sharded via `shard_bin_size`.

Inputs:
- `File vcf`: VCF to post-process.
- `File vcf_idx`: Index for VCF to post-process.
- `Array[String] contigs`: Contigs to process within the input VCF.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the per-record pass.
- `Boolean run_transfer_genotypes`: Whether to transfer genotypes from `transfer_vcf` onto heterozygous calls (run first; requires `transfer_vcf`).
- `Boolean run_unphase_samples`: Whether to unphase the samples in `unphase_samples` (requires `unphase_samples`).
- `Boolean run_normalize_ploidy`: Whether to normalize ploidy by sex - clearing chrY female calls, making chrX/chrY male calls hemizygous, enforcing diploidy and right-aligning unphased calls (requires `ped`).
- `Boolean run_decrement_trv_ids`: Whether to decrement tandem-repeat variant IDs.
- `Boolean run_prune_meis`: Whether to reclassify mobile elements whose length falls outside the expected bounds back to plain insertions/deletions.
- `Boolean run_flag_homopolymer_trvs`: Whether to flag tandem repeats with a length-1 shortest motif as `HOMOPOLYMER_TRV`.
- `Boolean run_sorting`: Whether to sort records sharing a coordinate by absolute allele length and variant ID.
- `Boolean run_filter_single_read_singletons`: Whether to apply the `SINGLE_READ_SUPPORT` filter to singleton calls.
- `Boolean run_filter_assembly_only_singletons`: Whether to apply the `ASSEMBLY_ONLY_SINGLETON` filter and emit a matching TSV.
- `File? transfer_vcf`: VCF whose genotypes are transferred when `run_transfer_genotypes` is set.
- `File? transfer_vcf_idx`: Index for `transfer_vcf`.
- `Array[String] unphase_samples`: Samples to unphase when `run_unphase_samples` is set (defaults to empty).
- `File? ped`: Cohort pedigree file, used for ploidy normalization when `run_normalize_ploidy` is set.

Outputs:
- `post_processed_vcf`: Post-processed VCF.
- `post_processed_vcf_idx`: Index for the post-processed VCF.
- `assembly_only_singletons_tsv`: Optional TSV containing one row per assembly-only singleton ALT allele; present only when `run_filter_assembly_only_singletons` is true.


### [QcAnnotations](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/QcAnnotations.wdl)
This workflow adapts the GATK-SV annotation QC pipeline in order to produce a quality-control report for an annotated callset. It collects VCF-wide site-level statistics per contig, converts the VCF to BED, plots the aggregated site metrics and - when comparison datasets are supplied - benchmarks the callset against them at the site level. It can additionally run a per-sample pass that collects per-sample variant lists, plots per-sample and per-family QC, and benchmarks samples against sample-level comparison datasets, before sanitizing all outputs into a single QC tarball.

Inputs:
- `Array[File] vcfs`: Annotated VCFs to QC.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs to QC.
- `File ped_file`: Cohort pedigree file, used for per-family QC.
- `Int variants_per_shard`: Number of variants per shard during VCF-wide collection.
- `Int samples_per_shard`: Number of samples per shard during per-sample collection.
- `Boolean create_variant_attributes`: Whether to compute variant-attribute breakdowns (default `false`).
- `Boolean run_per_sample`: Whether to run the per-sample QC pass (default `true`).
- `Int? subset_sample_count`: Number of samples to subset to before the per-sample pass.
- `String? subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCFs.
- `Int? random_seed`: Random seed for sample subsetting and downsampling.
- `Int? max_gq`: Maximum genotype-quality value used when binning QC metrics.
- `Int? downsample_qc_per_sample`: Number of variants to downsample to per sample during QC.
- `File? sample_renaming_tsv`: TSV mapping sample IDs to renamed IDs prior to QC.
- `Array[Array[String]]? site_level_comparison_datasets`: Site-level datasets to benchmark the callset against.
- `Array[Array[String]]? sample_level_comparison_datasets`: Sample-level datasets to benchmark the callset against.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `sv_vcf_qc_output`: Tarball of the QC plots and metrics.
- `vcf2bed_output`: Merged BED representation of the QC'd VCFs.


### [ResolveHaplotypeOverlaps](../wdl/annotation_utils/ResolveHaplotypeOverlaps.wdl)
This utility detects and resolves haplotype-level overlaps among non-TR, non-TR-enveloped variants in a phased cohort VCF. For each sample, it extracts the sample's non-ref calls (excluding `allele_type="trv"` and `INFO/TR_ENVELOPED` variants), then sweeps each haplotype's variant intervals to find all overlapping pairs. Overlapping pairs are resolved by keeping the variant that spans more reference sequence (larger `len(REF)`) - which always favors DELs over INS or SNVs. When two variants span the same reference length, the higher-GQ call wins; remaining ties are broken by `INFO/allele_length`, then type rank (DEL > INS > SNV), then QUAL, then input-file order. The loser's FORMAT fields (`GT`, `GQ`, `DP`, `EV`, `BEV`, `AD`, `PL`) are cleared in the output VCF. The workflow scatters per-sample detection across all samples, then applies the collected clears to the given contig (with optional record-count sharding) to produce the resolved VCF.

Inputs:
- `File vcf`: Phased cohort VCF to resolve.
- `File vcf_idx`: Index for `vcf`.
- `String contig`: Contig to process.
- `Int? records_per_shard`: When set, shards the contig into chunks of this many records for the clearing step.

Outputs:
- `overlap_resolved_vcf`: VCF with overlapping loser genotypes cleared.
- `overlap_resolved_vcf_idx`: Index for `overlap_resolved_vcf`.
- `overlap_tsv`: TSV of all detected overlap pairs, with columns `sample`, `haplotype`, `variant_id_retained`, `var_type_retained`, `size_bin_retained`, `variant_id_cleared`, `var_type_cleared`, `size_bin_cleared`.


### [SubsetTsvToColumns](../wdl/annotation_utils/SubsetTsvToColumns.wdl)
This utility subsets an annotation TSV to a chosen set of columns, optionally filtering rows to those whose columns match specified values. It outputs the subset TSV.

Inputs:
- `File annotations_tsv`: Annotation TSV to subset.
- `File annotations_header`: Header describing the TSV columns.
- `Array[String] subset_columns`: Columns to retain.
- `Array[Array[String]]? subset_column_values`: Per-column values to filter rows by.

Outputs:
- `subset_tsv`: Column-subset TSV.


### [TransformAlleleType](../wdl/annotation_utils/TransformAlleleType.wdl)
This utility reclassifies `allele_type` values and records the original type in a new `allele_subtype` field. Variants with `allele_type=dup` are tested for tandemness against their duplication source (from `INFO/ORIGIN`) using two criteria: size similarity between the insertion length and the ORIGIN region length must meet the `dup_size_similarity` threshold, and the insertion POS must fall within the ORIGIN region or within `dup_breakpoint_window` bases of its breakpoints. All get `allele_subtype=tandem_dup`; those passing keep `allele_type=dup`, while those failing are set to `allele_type=ins`. Variants with `allele_type` of `complex_dup`, `dup_interspersed`, `inv_dup`, `alu_ins`, `line_ins`, `sva_ins` or `numt` are set to `allele_type=ins`, and those with `alu_del`, `line_del` or `sva_del` are set to `allele_type=del`, each recording the original value in `allele_subtype`. REF/ALT/POS are never modified. Records with other `allele_type` values are passed through unchanged. Supports optional record-count sharding.

Inputs:
- `File vcf`: VCF to transform.
- `File vcf_idx`: Index for `vcf`.
- `Int? records_per_shard`: Number of records per shard for parallel processing.
- `Int dup_breakpoint_window`: Maximum distance (bp) between insertion POS and ORIGIN breakpoints to pass the breakpoint check (default `10`).
- `Float dup_size_similarity`: Minimum size similarity ratio (relative to the larger of the two lengths) between insertion and ORIGIN lengths (default `0.9`).
- `Int min_dup_size`: Minimum insertion size (bp) to consider for the tandem check (default `50`).

Outputs:
- `transformed_vcf`: VCF with revised `allele_type`/`allele_subtype`.
- `transformed_vcf_idx`: Index for `transformed_vcf`.


## Tools

### [Automop](../wdl/tools/Automop.wdl)
This tool runs `mop` (via FISS) to clean up unreferenced intermediate files in a Terra workspace, freeing storage. A dry-run mode reports what would be deleted without removing anything.

Inputs:
- `String workspace_namespace`: Terra workspace namespace to clean.
- `String workspace_name`: Terra workspace name to clean.
- `String user`: User running the cleanup.
- `Boolean dry_run`: Whether to report rather than perform deletions.

Outputs:
- `fissfc_log`: Log of the cleanup run.


### [BackbonePhase](../wdl/tools/BackbonePhase.wdl)
This tool transfers ("backbone") phasing from a set of base VCFs onto a target VCF for a single contig. It assigns each sample to its base VCF, computes the phase-flip orientation needed to make the target consistent with the backbone, and applies those flips. It outputs the phase-transferred VCF and a list of samples with no matching base VCF.

Inputs:
- `File vcf`: Target VCF to phase.
- `File vcf_idx`: Index for `vcf`.
- `Array[File] base_vcfs`: Base VCFs providing the backbone phasing.
- `Array[File] base_vcf_idxs`: Indexes for `base_vcfs`.
- `String contig`: Contig to process.
- `File? swap_samples_base`: Sample-ID swap map applied to the base VCFs.
- `Boolean allow_unphased_match_phase`: Whether to allow unphased genotypes to set the phase orientation (default `false`).

Outputs:
- `transferred_vcf`: Phase-transferred VCF.
- `transferred_vcf_idx`: Index for the transferred VCF.
- `missing_samples`: Samples with no matching base VCF.


### [HiFiCNV](../wdl/tools/HiFiCNV.wdl)
This tool runs PacBio [HiFiCNV](https://github.com/PacificBiosciences/HiFiCNV) on a sample's aligned HiFi BAM to call copy number variants from read depth. It outputs the CNV VCF, a copy-number bedgraph, a depth BigWig track and the tool's log.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `String sex`: Sex of sample (one of `M` or `F`), used to select the matching expected-CN file.
- `File ref_fa`: Reference sequences FASTA file.
- `File ref_fai`: Index for `ref_fa`.
- `File exclude_bed`: Regions to exclude from CNV calling (e.g. centromeres).
- `File expected_cn_male`: PAR regions and expected copy numbers for sex chromosomes, male.
- `File expected_cn_female`: PAR regions and expected copy numbers for sex chromosomes, female.

Outputs:
- `hificnv_vcf`: CNV calls VCF.
- `hificnv_vcf_idx`: Index for the CNV calls VCF.
- `hificnv_bedgraph`: Per-window copy number bedgraph.
- `hificnv_depth_bw`: Depth BigWig track.
- `hificnv_log`: HiFiCNV log file.


### [HiPhase](../wdl/tools/HiPhase.wdl)
This tool runs PacBio [HiPhase](https://github.com/PacificBiosciences/HiPhase) to jointly phase a sample's small-variant, SV and (optionally) TRGT VCFs against its aligned reads. It preprocesses and synchronizes the input VCFs per contig, phases them together and optionally haplotags the BAM. It outputs the phased VCF, per-contig phasing statistics and an optional haplotagged BAM.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `File small_vcf`: Small-variant (SNV/indel) VCF to phase.
- `File small_vcf_idx`: Index for `small_vcf`.
- `File sv_vcf`: SV VCF to phase.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `Array[String] contigs`: Contigs to phase.
- `File? trgt_vcf`: TRGT tandem-repeat VCF to additionally phase.
- `File? trgt_vcf_idx`: Index for `trgt_vcf`.
- `Int? trgt_min_repeat_unit`: Minimum repeat-unit length retained when filtering the TRGT VCF.
- `Boolean? trgt_normalize`: Whether to normalize the TRGT VCF before phasing.
- `Int? trgt_min_length_diff`: Minimum length difference retained when filtering the TRGT VCF.
- `Int? trgt_max_catalog_length`: Maximum catalog length retained when filtering the TRGT VCF.
- `String? hiphase_extra_args`: Additional arguments passed to HiPhase.
- `Boolean run_haplotagging`: Whether to also haplotag the BAM (default `false`).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `hiphase_vcf`: Phased VCF.
- `hiphase_vcf_idx`: Index for the phased VCF.
- `hiphase_haplotag_files`: Per-contig haplotag read assignments.
- `hiphase_stats`: Per-contig phasing statistics.
- `hiphase_blocks`: Per-contig phase blocks.
- `hiphase_summary`: Per-contig phasing summaries.
- `hiphase_haplotagged_bam`: Haplotagged BAM (only when `run_haplotagging`).
- `hiphase_haplotagged_bam_idx`: Index for the haplotagged BAM (only when `run_haplotagging`).


### [HiPhaseMerge](../wdl/tools/HiPhaseMerge.wdl)
This tool merges per-sample HiPhase-phased VCFs into a cohort VCF on a per-contig basis, optionally also merging the TRGT tandem-repeat calls separately - fixing TRGT `END`/`AL` headers and propagating phase-set tags. It outputs the merged integrated VCF and an optional merged TRGT VCF.

Inputs:
- `Array[File] phased_vcfs`: Per-sample HiPhase-phased VCFs to merge.
- `Array[File] phased_vcf_idxs`: Indexes for `phased_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Boolean merge_trgt`: Whether to additionally merge the TRGT tandem-repeat calls separately.
- `String merge_args`: Arguments controlling the VCF merge (default `--merge id`).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `hiphase_merged_integrated_vcf`: Merged integrated cohort VCF.
- `hiphase_merged_integrated_vcf_idx`: Index for the merged integrated VCF.
- `hiphase_merged_trgt_vcf`: Merged TRGT VCF (only when `merge_trgt`).
- `hiphase_merged_trgt_vcf_idx`: Index for the merged TRGT VCF (only when `merge_trgt`).


### [Kanpig](../wdl/tools/Kanpig.wdl)
This tool regenotypes a cohort SV VCF against each sample's aligned reads using [Kanpig](https://github.com/ACEnglish/kanpig). It subsets the cohort to the target samples, runs Kanpig per sample with sex-aware ploidy beds, and merges the per-sample genotypes back into both a raw and a processed cohort VCF. It outputs the regenotyped (processed) and raw Kanpig VCFs.

Inputs:
- `File cohort_vcf`: Cohort SV VCF to regenotype.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[File] bams`: Aligned reads, one per sample.
- `Array[File] bais`: Indexes for `bams`.
- `Array[String] sample_ids`: Samples to regenotype.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `String merge_args`: Arguments controlling the per-sample genotype merge (default `--merge id`).
- `String kanpig_params`: Parameters passed to Kanpig (default `--neighdist 500 --gpenalty 0.04 --hapsim 0.97`).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).
- `File ploidy_bed_male`: From [references](references.md).
- `File ploidy_bed_female`: From [references](references.md).

Outputs:
- `sv_kanpig_vcf`: Regenotyped (processed) cohort VCF.
- `sv_kanpig_vcf_idx`: Index for the processed VCF.
- `sv_kanpig_raw_vcf`: Raw Kanpig cohort VCF.
- `sv_kanpig_raw_vcf_idx`: Index for the raw VCF.


### [LRCNVs](../wdl/utils/LRCNVs.wdl)
This tool calls copy-number variants across a cohort using the GATK germline CNV (gCNV) pipeline in cohort mode. From per-sample depth profiles over a shared interval list it annotates and filters intervals, determines contig ploidy, fits the gCNV model across scattered interval shards, post-processes the per-sample calls into genotyped interval and segment VCFs, and collects sample- and model-level QC. It outputs the gCNV model, per-sample CNV VCFs, denoised copy ratios and QC status.

Inputs:
- `File intervals`: Interval list over which CNVs are called.
- `Array[String]+ entity_ids`: Sample IDs in the cohort.
- `Array[String]+ depth_profiles`: Per-sample read-depth profiles, aligned to `entity_ids`.
- `String cohort_entity_id`: Identifier for the cohort.
- `File contig_ploidy_priors`: Contig ploidy priors used to determine per-sample contig ploidy.
- `Int num_intervals_per_scatter`: Number of intervals processed per scatter shard (default `10000`).
- `File? gatk4_jar_override`: Override GATK4 jar.
- `File? mappability_track_bed`: Mappability track used to annotate intervals.
- `File? mappability_track_bed_idx`: Index for `mappability_track_bed`.
- `File? segmental_duplication_track_bed`: Segmental-duplication track used to annotate intervals.
- `File? segmental_duplication_track_bed_idx`: Index for `segmental_duplication_track_bed`.
- `Int? feature_query_lookahead`: Base pairs to look ahead when querying interval-annotation feature tracks.
- `File? blacklist_intervals`: Intervals to exclude from calling.
- `Int? low_count_filter_count_threshold`: Minimum read count for an interval to be considered well-covered in a sample.
- `Float? low_count_filter_percentage_of_samples`: Minimum percentage of samples that must meet `low_count_filter_count_threshold` for an interval to pass.
- `Float? extreme_count_filter_minimum_percentile`: Lower count percentile below which an interval is considered an outlier.
- `Float? extreme_count_filter_maximum_percentile`: Upper count percentile above which an interval is considered an outlier.
- `Float? extreme_count_filter_percentage_of_samples`: Minimum percentage of samples that must pass the extreme-count percentile bounds for an interval to pass.
- `Int ref_copy_number_autosomal_contigs`: Reference copy number for autosomes (default `2`).
- `Array[String]? allosomal_contigs`: Contigs treated as allosomal.
- `Int maximum_number_events_per_sample`: Maximum number of events permitted per sample (default `1000`).
- The workflow additionally exposes numerous optional gCNV model and contig-ploidy hyperparameters, prefixed `gcnv_` and `ploidy_`, that tune the underlying GATK tasks.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).
- `File ref_dict`: From [references](references.md).

Outputs:
- `annotated_intervals`: Intervals annotated with GC content and tracks.
- `filtered_intervals`: Intervals retained after filtering.
- `contig_ploidy_model_tar`: Fitted contig-ploidy model.
- `contig_ploidy_calls_tar`: Per-sample contig-ploidy calls.
- `gcnv_model_tars`: Fitted gCNV models, one per scatter shard.
- `gcnv_calls_tars`: Per-shard per-sample gCNV calls.
- `gcnv_tracking_tars`: Per-shard model-fitting tracking files.
- `genotyped_intervals_vcfs`: Per-sample genotyped interval VCFs.
- `genotyped_segments_vcfs`: Per-sample genotyped segment VCFs.
- `sample_qc_status_files`: Per-sample QC status files.
- `sample_qc_status_strings`: Per-sample QC status strings.
- `model_qc_status_file`: Model-level QC status file.
- `model_qc_string`: Model-level QC status string.
- `denoised_copy_ratios`: Per-sample denoised copy ratios.


### [MethylationProfiling](../wdl/tools/MethylationProfiling.wdl)
This tool generates CpG methylation pileups from a haplotagged BAM using [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools), producing combined and per-haplotype methylation BED tracks.

Inputs:
- `File bam`: Haplotagged aligned reads.
- `File bai`: Index for `bam`.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `cpg_combined_bed`: Combined methylation pileup BED.
- `cpg_combined_bed_idx`: Index for the combined BED.
- `cpg_hap1_bed`: Haplotype 1 methylation pileup BED.
- `cpg_hap1_bed_idx`: Index for the haplotype 1 BED.
- `cpg_hap2_bed`: Haplotype 2 methylation pileup BED.
- `cpg_hap2_bed_idx`: Index for the haplotype 2 BED.


### [MinimapAlignment](../wdl/tools/MinimapAlignment.wdl)
This workflow leverages [Minimap2](https://github.com/lh3/minimap2) in order to align a sample's maternal and paternal assemblies to a reference.

Inputs:
- `File assembly_mat`: Maternal assembly.
- `File assembly_pat`: Paternal assembly.
- `String sample_id`: ID of the sample being aligned.
- `String minimap_flags`: Parameters to use when running Minimap2 (default `-a -x asm20 --cs --eqx`).
- `Int minimap_threads`: Number of alignment threads (default `32`).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `minimap_assembled_bam_mat`: Aligned maternal-assembly BAM.
- `minimap_assembled_bai_mat`: Index for the maternal BAM.
- `minimap_assembled_paf_mat`: Maternal-assembly PAF alignment.
- `minimap_assembled_bam_pat`: Aligned paternal-assembly BAM.
- `minimap_assembled_bai_pat`: Index for the paternal BAM.
- `minimap_assembled_paf_pat`: Paternal-assembly PAF alignment.


### [MosDepth](../wdl/tools/MosDepth.wdl)
This tool runs [mosdepth](https://github.com/brentp/mosdepth) to compute sequencing depth over a sample's BAM per contig. By default it emits per-base coverage; when `bin_size` is set, it instead windows depth into fixed-size bins (`--by`, `--no-per-base`) and emits per-region coverage.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `Array[String] contigs`: Contigs over which to compute depth.
- `Int? bin_size`: If set, windows depth into bins of this size (bp) and disables per-base output.
- `File? ref_fa`: Reference FASTA, required for CRAM input.
- `File? ref_fai`: Index for `ref_fa`.

Outputs:
- `mosdepth_dist`: Per-contig cumulative coverage distributions.
- `mosdepth_summary`: Per-contig coverage summaries.
- `mosdepth_per_base`: Per-contig per-base coverage (when `bin_size` is unset).
- `mosdepth_per_base_csi`: Indexes for the per-base coverage.
- `mosdepth_regions_bed`: Per-contig windowed coverage BEDs (when `bin_size` is set).
- `mosdepth_regions_bed_csi`: Indexes for the windowed coverage BEDs.


### [PALMER](../wdl/tools/PALMER.wdl)
This workflow runs PALMER on a pair of aligned assembly haplotypes in order to generate MEI calls. It then convets the raw PALMER calls generated into a VCF, merges calls across the haplotypes to create a diploid VCF per haplotype and then finally integrates these into a final VCF containing multiple MEI types.

Inputs:
- `File? bam_pat`: Aligned assembly for paternal haplotype.
- `File? bai_pat`: Index for `bam_pat`.
- `File? bam_mat`: Aligned assembly for maternal haplotype.
- `File? bai_mat`: Index for `bam_mat`.
- `Array[File]? override_palmer_calls_pat`: Optional precomputed PALMER calls for the paternal haplotype, causing the workflow to bypass execution.
- `Array[File]? override_palmer_tsd_files_pat`: Optional precomputed PALMER TSD files for the paternal haplotype, causing the workflow to bypass execution.
- `Array[File]? override_palmer_calls_mat`: Optional precomputed PALMER calls for the maternal haplotype, causing the workflow to bypass execution.
- `Array[File]? override_palmer_tsd_files_mat`: Optional precomputed PALMER TSD files for the maternal haplotype, causing the workflow to bypass execution.
- `Array[String] contigs`: Contigs to run PALMER on.
- `String sample`: ID of the sample being processed.
- `String mode`: PALMER run mode.
- `Array[String] mei_types`: MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `Array[String]? truvari_collapse_params`: Per-MEI-type Truvari parameters used when merging calls across haplotypes (default `--pctsize 0.9 --pctovl 0.9 --pctseq 0.9 --refdist 500` for each type).
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `palmer_pat_calls`: Raw PALMER calls for the paternal haplotype, per MEI type.
- `palmer_pat_tsd_reads`: PALMER TSD reads for the paternal haplotype, per MEI type.
- `palmer_pat_vcfs`: Paternal-haplotype PALMER VCFs, per MEI type.
- `palmer_pat_vcf_idxs`: Indexes for the paternal-haplotype VCFs.
- `palmer_mat_calls`: Raw PALMER calls for the maternal haplotype, per MEI type.
- `palmer_mat_tsd_reads`: PALMER TSD reads for the maternal haplotype, per MEI type.
- `palmer_mat_vcfs`: Maternal-haplotype PALMER VCFs, per MEI type.
- `palmer_mat_vcf_idxs`: Indexes for the maternal-haplotype VCFs.
- `palmer_diploid_vcfs`: Diploid PALMER VCFs merged across haplotypes, per MEI type.
- `palmer_diploid_vcf_idxs`: Indexes for the diploid VCFs.
- `palmer_combined_vcf`: Final VCF combining all MEI types.
- `palmer_combined_vcf_idx`: Index for the combined VCF.


### [PALMERDiploid](../wdl/tools/PALMERDiploid.wdl)
This tool runs [PALMER](https://github.com/WeichenZhou/PALMER) on a single sample to generate mobile element insertion calls and convert them to a VCF. It shards the input BAM, runs PALMER per MEI type, merges the shard outputs and converts the raw calls into a per-type VCF, optionally bypassing execution when PALMER calls are supplied directly. It outputs the raw PALMER call and TSD files, per-type VCFs and a combined VCF.

Inputs:
- `File? bam`: Aligned reads to run PALMER on.
- `File? bai`: Index for `bam`.
- `Array[File]? override_palmer_calls`: Optional precomputed PALMER calls, causing the workflow to bypass execution.
- `Array[File]? override_palmer_tsd_files`: Optional precomputed PALMER TSD files, causing the workflow to bypass execution.
- `Array[String] contigs`: Contigs to run PALMER on.
- `String sample`: ID of the sample being processed.
- `String mode`: PALMER run mode.
- `Array[String] mei_types`: MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `palmer_calls`: Raw PALMER calls, per MEI type.
- `palmer_tsd_reads`: PALMER TSD reads, per MEI type.
- `palmer_diploid_vcfs`: Per-MEI-type PALMER VCFs.
- `palmer_diploid_vcf_idxs`: Indexes for the per-type VCFs.
- `palmer_combined_vcf`: Final VCF combining all MEI types.
- `palmer_combined_vcf_idx`: Index for the combined VCF.


### [PALMERMerge](../wdl/tools/PALMERMerge.wdl)
This tool merges multiple PALMER MEI VCFs into a single VCF per contig and concatenates the result across contigs. It outputs the merged PALMER VCF.

Inputs:
- `Array[File] vcfs`: PALMER VCFs to merge.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs to process.

Outputs:
- `palmer_merged_vcf`: Merged PALMER VCF.
- `palmer_merged_vcf_idx`: Index for the merged VCF.


### [PAV](../wdl/tools/PAV.wdl)
This tool runs [PAV](https://github.com/EichlerLab/pav) in batch mode across multiple samples' phased haplotype assemblies to call variants against the reference. It outputs per-sample VCFs, along with tarballs of the full PAV results and log directories.

Inputs:
- `Array[File] mat_haplotypes`: Maternal haplotype assemblies, one per sample.
- `Array[File] pat_haplotypes`: Paternal haplotype assemblies, one per sample.
- `Array[String] sample_ids`: Sample IDs, aligned by index to `mat_haplotypes`/`pat_haplotypes`.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `pav_results_tarball`: Tarball of the full PAV results directory.
- `pav_log_tarball`: Tarball of the full PAV log directory.
- `pav_vcfs`: Per-sample called VCFs.
- `pav_vcf_indices`: Indexes for the per-sample VCFs.
- `debug_sam`: Optional debug alignment file.
- `debug_temp`: Optional debug intermediate files.


### [RepeatMasker](../wdl/tools/RepeatMasker.wdl)
This workflow leverages [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) in order to annotate repeated and mobile-element content in the insertions of an input VCF. It extracts each insertion's inserted sequence to a FASTA, optionally restricted to a minimum length, and runs RepeatMasker over it.

Inputs:
- `File vcf`: VCF whose insertions are masked.
- `File vcf_idx`: Index for VCF.
- `Int? min_length`: Minimum insertion length to extract and mask.

Outputs:
- `rm_out`: RepeatMasker output table.
- `rm_fa`: FASTA of the masked insertion sequences.


### [TRGT](../wdl/tools/TRGT.wdl)
This workflow leverages [TRGT](https://github.com/PacificBiosciences/trgt) in order to genotype short-tandem repeats.

Inputs:
- `File bam`: Aligned reads.
- `File bai`: Index for aligned reads.
- `String sample_id`: ID of the sample being genotyped.
- `String sex`: Sex of sample (one of `M` or `F`).
- `String catalog_name`: Name of the repeat catalog used, included in the output VCF filename.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).
- `File repeat_catalog_trgt`: From [references](references.md).

Outputs:
- `trgt_vcf`: TRGT tandem-repeat genotype VCF.
- `trgt_vcf_idx`: Index for the TRGT VCF.


### [TRGTLPS](../wdl/tools/TRGTLPS.wdl)
This tool runs the [trgt-lps](https://github.com/PacificBiosciences/trgt-lps) tool per contig to compute the longest polymer sequence (LPS) within each TRGT-genotyped tandem-repeat locus for every sample, concatenating the results into a single TSV. It outputs the LPS TSV.

Inputs:
- `File vcf`: TRGT VCF to process.
- `File vcf_idx`: Index for the TRGT VCF.
- `Array[String] contigs`: Contigs to process.

Outputs:
- `trgt_lps_tsv`: TSV of per-locus longest polymer sequences.


### [Vamos](../wdl/tools/Vamos.wdl)
This tool runs [Vamos](https://github.com/ChaissonLab/vamos) in order to genotype tandem repeats against a Vamos repeat catalog, in read mode (from an aligned read BAM) and/or assembly mode (from per-haplotype assembly BAMs). It outputs the resulting Vamos VCFs.

Inputs:
- `File? read_bam`: Aligned reads to genotype in read mode.
- `File? read_bai`: Index for `read_bam`.
- `Array[File]? assembly_bams`: Per-haplotype assembly BAMs to genotype in assembly mode.
- `Array[File]? assembly_bais`: Indexes for `assembly_bams`.
- `File repeat_catalog_vamos`: Vamos repeat catalog to genotype against.
- `String sample_id`: ID of the sample being genotyped.

Outputs:
- `vamos_assembly_vcfs`: Per-haplotype assembly-mode Vamos VCFs.
- `vamos_assembly_vcf_idxs`: Indexes for the assembly-mode VCFs.
- `vamos_reads_vcf`: Read-mode Vamos VCF.
- `vamos_reads_vcf_idx`: Index for the read-mode VCF.


### [VcfDist](../wdl/tools/VcfDist.wdl)
This tool runs [vcfdist](https://github.com/TimD1/vcfdist) in order to benchmark an evaluation VCF against a truth VCF per contig, computing alignment-based precision/recall and phasing accuracy. It outputs vcfdist's precision-recall, phasing, switch-flip, phase-block and supercluster reports.

Inputs:
- `File vcf_eval`: VCF being evaluated.
- `File vcf_eval_idx`: Index for `vcf_eval`.
- `File vcf_truth`: Truth VCF to evaluate against.
- `File vcf_truth_idx`: Index for `vcf_truth`.
- `Array[String] contigs`: Contigs to evaluate.
- `File? bed_regions`: BED of regions to restrict the evaluation to.
- `String? mode`: vcfdist evaluation mode.
- `Float? threshold`: vcfdist matching threshold.
- `String? vcfdist_args`: Additional arguments passed to vcfdist.
- `File ref_fa`: From [references](references.md).

Outputs:
- `vcfdist_phasing_summary_tsv`: Per-contig phasing summaries.
- `vcfdist_switchflips_tsv`: Per-contig switch and flip errors.
- `vcfdist_precision_recall_tsv`: Per-contig precision-recall curves.
- `vcfdist_precision_recall_summary_tsv`: Per-contig precision-recall summaries.
- `vcfdist_phase_blocks_tsv`: Per-contig phase blocks.
- `vcfdist_superclusters_tsv`: Per-contig variant superclusters.
- `vcfdist_query_tsv`: Per-contig query-variant results.
- `vcfdist_truth_tsv`: Per-contig truth-variant results.
- `vcfdist_summary_vcf`: Per-contig annotated summary VCFs.


### [VcfDistCohort](../wdl/tools/VcfDistCohort.wdl)
This tool runs [vcfdist](https://github.com/TimD1/vcfdist) across a cohort by pairing each evaluation VCF with its corresponding truth VCF and benchmarking every assigned sample, then aggregating the per-sample results. It outputs cohort-level precision/recall and phasing summaries.

Inputs:
- `Array[File] eval_vcfs`: Evaluation VCFs, one per group.
- `Array[File] eval_vcf_idxs`: Indexes for `eval_vcfs`.
- `Array[File] truth_vcfs`: Truth VCFs, aligned to `eval_vcfs`.
- `Array[File] truth_vcf_idxs`: Indexes for `truth_vcfs`.
- `Array[String] contigs`: Contigs to evaluate.
- `Array[String]? subset_samples`: Samples to restrict the evaluation to.
- `String? vcfdist_args`: Additional arguments passed to vcfdist.
- `File ref_fa`: From [references](references.md).

Outputs:
- `vcfdist_phasing_summary_tsv`: Cohort phasing summary.
- `vcfdist_precision_recall_summary_tsv`: Cohort precision-recall summary.
- `vcfdist_precision_recall_tsv`: Cohort precision-recall curves.
- `vcfdist_switchflips_tsv`: Cohort switch and flip errors.
- `vcfdist_phase_blocks_tsv`: Cohort phase blocks.
- `vcfdist_missing_samples`: Samples with no matching truth VCF.


### [Whatshap](../wdl/tools/Whatshap.wdl)
This tool haplotags a sample's BAM against a phased VCF using [WhatsHap](https://github.com/whatshap/whatshap), per contig, then merges the tagged reads into a single BAM. It outputs the haplotagged BAM and per-contig haplotag read lists.

Inputs:
- `File bam`: Aligned reads to haplotag.
- `File bai`: Index for `bam`.
- `File phased_vcf`: Phased VCF used to assign haplotypes.
- `File phased_vcf_idx`: Index for `phased_vcf`.
- `Array[String] contigs`: Contigs to process.
- `String? extra_args`: Additional arguments passed to WhatsHap.
- `File ref_fa`: From [references](references.md).
- `File ref_fai`: From [references](references.md).

Outputs:
- `haplotagged_bam`: Haplotagged BAM.
- `haplotagged_bai`: Index for the haplotagged BAM.
- `haplotag_lists`: Per-contig haplotag read assignments.
