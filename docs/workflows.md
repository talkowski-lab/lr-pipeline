# Workflows
This document describes each WDL workflow in the pipeline, including its purpose, inputs and outputs. Annotations, annotation utilities and tools are run directly and are registered in `.dockstore.yml`; the sub-workflows in the final section are imported building blocks that are never run on their own.

This file is generated from the `meta` and `parameter_meta` blocks of each workflow by [`generate_workflows_doc.py`](../.github/scripts/generate_workflows_doc.py). Edit those blocks rather than this document. Inputs described as `From references.` are the shared reference files listed in [references](references.md).


## Annotations


### [AnnotateAgeMetrics](../wdl/annotation/AnnotateAgeMetrics.wdl)
This workflow computes the age distribution of carriers for every variant in the input VCF. For each sample it derives an age from a date-of-birth table relative to a fixed reference date, then tabulates the number of heterozygous and homozygous carriers of each allele that fall into a set of user-defined age bins, along with overflow `smaller` and `larger` bins for ages outside the configured range. It emits a TSV of these per-allele age-bin counts.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File age_data`: CSV file with `person_id` and `date_of_birth` columns, used to derive each sample's age.
- `Array[Int] age_bins`: Age-bin edges, in years, into which carrier ages are binned.
- `String reference_date`: Reference date (`YYYY-MM-DD`) against which each sample's age is computed.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotations_tsv_age`: TSV of per-allele carrier counts across the age bins.

### [AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl)
This workflow ingests a callset VCF and two truth VCFs - one of SNVs & indels and one of SVs - and finds matching variants across them, annotating each matched callset variant with the truth callset's AC/AF/AN and genotype-count fields. This enables benchmarking annotations against an existing cohort (e.g. gnomAD) and surfacing variants that are outliers relative to it.

The workflow undergoes multiple rounds of variant matching in order to determine matched pairs: (1) Exact match across CHROM, POS, REF and ALT. (2) Truvari match with overlap percentages of 90%, 70% and 50%. (3) Matching based on `bedtools closest`, finetuned for SVs. Here the callset and truth variants are split by type and converted to a symbolic representation, after which separate `bedtools closest` passes are run - one tuned for deletions and duplications via reciprocal positional overlap, and one tuned for insertions via breakpoint proximity - so that each callset variant is paired with the nearest same-type truth variant above the per-callset minimum SV-length thresholds.

Note: When converting to symbolic representation, only canonical DUPs (allele_type = `DUP` exactly) are treated as DUP; other DUP subtypes (e.g., `dup_interspersed`, `inv_dup`) are treated as insertions.

Note: Callset DUPs are compared twice, because the two matching rules need different coordinates. Against truth DUPs they are repositioned to their `ORIGIN` coordinates and compared by reciprocal overlap; against truth insertions they are held at their insertion site and compared by breakpoint proximity and length ratio.

Note: The SV truth VCF is expected to be symbolic already. Set `convert_symbolic_truth_sv_vcf` when it instead carries sequence alleles in the same format as the callset, in which case it is converted with its DUPs repositioned onto their `ORIGIN` coordinates, matching how truth DUPs are positioned in a symbolic truth callset.

Both the exact-match and Truvari rounds can be sharded within a contig. Truvari shard boundaries are snapped forward to the next gap wider than the `min_shard_gap_truvari_match` input of `TruvariMatch`, which keeps results identical to an unsharded run because Truvari only groups records into a new comparison chunk once the next record clears the running end by more than its chunk size. Fixed-width bins alone would split colocated record pairs and silently lose matches.

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
- `Int? shard_bin_size_exact_match`: If set, shards the exact-match round into contig regions of roughly this many base pairs, run in parallel.
- `Int? shard_bin_size_truvari_match`: If set, shards the Truvari round into contig regions of at least this many base pairs, run in parallel. Each region is extended to the next safe gap, so a value of 1000000 or more is recommended.
- `Boolean convert_symbolic_truth_sv_vcf`: Whether the SV truth VCF represents alleles as sequence rather than symbolically. When true it is converted to a symbolic representation first, reading the same `type_field_vcf` and `length_field_vcf` INFO fields as the callset. (default `false`)
- `String type_field_vcf`: INFO field in the callset VCF giving each variant's allele type. (default `allele_type`)
- `String length_field_vcf`: INFO field in the callset VCF giving each variant's allele length. (default `allele_length`)
- `String source_tag_truth_snv_indel_vcf`: Label used to tag matches against the SNV & indel truth VCF. (default `SNV_indel`)
- `String source_tag_truth_sv_vcf`: Label used to tag matches against the SV truth VCF. (default `SV`)
- `String? args_string_vcf`: `bcftools view` arguments used to pre-subset the callset VCF.
- `String? args_string_truth_snv_indel_vcf`: `bcftools view` arguments used to pre-subset the SNV & indel truth VCF.
- `String? args_string_truth_sv_vcf`: `bcftools view` arguments used to pre-subset the SV truth VCF.
- `String? rename_id_string_vcf`: Expression used to rename variant IDs in the callset VCF prior to matching.
- `String? rename_id_string_truth_snv_indel_vcf`: Expression used to rename variant IDs in the SNV & indel truth VCF prior to matching.
- `String? rename_id_string_truth_sv_vcf`: Expression used to rename variant IDs in the SV truth VCF prior to matching.
- `Boolean? rename_id_strip_chr_vcf`: Whether to strip the `chr` prefix when renaming callset variant IDs.
- `Boolean? rename_id_strip_chr_truth_snv_indel_vcf`: Whether to strip the `chr` prefix when renaming SNV & indel truth variant IDs.
- `Boolean? rename_id_strip_chr_truth_sv_vcf`: Whether to strip the `chr` prefix when renaming SV truth variant IDs.
- `File? ref_fa`: From references. Only needed when either VCF represents alleles symbolically, since Truvari uses it solely to resolve those alleles to sequence.
- `File? ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String gatk_sv_lr_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (37).

Outputs:
- `File annotations_tsv_benchmark`: TSV mapping callset variants to their matched truth variants, match type, and the truth callset's AC/AF/AN and genotype-count fields.
- `File annotations_header_benchmark`: Header listing the extra annotation columns present in `annotations_tsv_benchmark`.

### [AnnotateDbSNP](../wdl/annotation/AnnotateDbSNP.wdl)
This workflow annotates each variant in the input VCF with its dbSNP reference SNP identifier (rsID). It matches variants against a per-contig dbSNP VCF on CHROM, POS, REF and ALT, emitting a TSV mapping each matched variant to its `dbSNP_ID`.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File dbsnp_vcf`: From references.
- `File dbsnp_vcf_idx`: From references.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? shard_bin_size`: If set, shards each contig into regions each containing roughly this many combined VCF + dbSNP records, run in parallel.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File annotations_tsv_dbsnp`: TSV mapping variants to their dbSNP identifiers.

### [AnnotateDbVaR](../wdl/annotation/AnnotateDbVaR.wdl)
This workflow annotates structural variants in the input VCF with matching records from dbVar. It restricts to variants at or above a minimum length, converts them to a symbolic representation, and matches deletions, duplications and insertions separately against a per-contig dbVar VCF using type-specific size-similarity, reciprocal-overlap and breakpoint-window thresholds. It emits a TSV linking matched variants to their dbVar records.

Note: When converting to symbolic representation, only canonical DUPs (allele_type = `DUP` exactly) are treated as DUP; other DUP subtypes (e.g., `dup_interspersed`, `inv_dup`) are treated as insertions.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File dbvar_vcf`: From references.
- `File dbvar_vcf_idx`: From references.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int min_length`: Minimum variant length to consider for matching.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching deletions. (default `500`)
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for matching deletions. (default `0.7`)
- `Float del_size_similarity`: Minimum size similarity for matching deletions. (default `0.7`)
- `Int dup_breakpoint_window`: Breakpoint window, in bp, for matching duplications. (default `500`)
- `Float dup_reciprocal_overlap`: Minimum reciprocal overlap for matching duplications. (default `0.7`)
- `Float dup_size_similarity`: Minimum size similarity for matching duplications. (default `0.7`)
- `Int ins_breakpoint_window`: Breakpoint window, in bp, for matching insertions. (default `200`)
- `Float ins_reciprocal_overlap`: Minimum reciprocal overlap for matching insertions. (default `0.0`)
- `Float ins_size_similarity`: Minimum size similarity for matching insertions. (default `0.5`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File annotations_tsv_dbvar`: TSV mapping variants to their matched dbVar records.

### [AnnotateIndelTRs](../wdl/annotation/AnnotateIndelTRs.wdl)
This workflow flags short insertions and deletions that represent tandem repeats. Using the str-analysis `filter_vcf_to_tandem_repeats` tool, it inspects each indel's sequence and marks it as a tandem repeat when it meets a minimum total repeat length, minimum number of repeats and minimum repeat-unit length, emitting a TSV of the flagged variants.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCF before tandem-repeat filtering. (default `-i 'INFO/allele_type!=\"trv\" && INFO/TR_ENVELOPED!=1'`)
- `Int min_tandem_repeat_length`: Minimum total tandem-repeat length for an indel to be flagged. (default `9`)
- `Int min_repeats`: Minimum number of repeats for an indel to be flagged. (default `3`)
- `Int min_repeat_unit_length`: Minimum repeat-unit length for an indel to be flagged. (default `1`)
- `String prefix`: Prefix for output file names.
- `String stranalysis_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotations_tsv_trs`: TSV of indels flagged as tandem repeats.

### [AnnotateInSilicoPredictors](../wdl/annotation/AnnotateInSilicoPredictors.wdl)
This workflow annotates SNVs and indels with precomputed in-silico predictor scores - CADD, Pangolin, PhyloP, REVEL and SpliceAI - drawn from the gnomAD V4 Hail Tables. It shards the VCF and uses a Hail-based script to look up each variant's scores, emitting a TSV of per-variant predictions.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String cadd_ht`: From references.
- `String pangolin_ht`: From references.
- `String phylop_ht`: From references.
- `String revel_ht`: From references.
- `String spliceai_ht`: From references.
- `String annotate_in_silico_predictors_script`: Path to the Hail script that performs the lookups (defaults to this repository's copy on `main`). (default `https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/annotation/annotate_insilico_predictors.py`)
- `String genome_build`: Genome build to annotate against. (default `GRCh38`)
- `String prefix`: Prefix for output file names.
- `String hail_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotations_tsv_insilico`: TSV of per-variant in-silico predictor scores.

### [AnnotateGnomADSTR](../wdl/annotation/AnnotateGnomADSTR.wdl)
This workflow annotates tandem-repeat variants in the input VCF with overlapping loci from the gnomAD V4 tandem-repeat catalog. It subsets to tandem-repeat calls and matches each against the catalog using a minimum reciprocal-overlap threshold, emitting a TSV linking calls to their gnomAD TR locus.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File gnomad_tr_json`: From references.
- `Float trv_reciprocal_overlap`: Minimum reciprocal overlap between a call and a catalog locus to be matched. (default `0.7`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File annotations_tsv_gnomad_str`: TSV mapping tandem-repeat calls to their gnomAD TR loci.

### [AnnotateGQMetrics](../wdl/annotation/AnnotateGQMetrics.wdl)
This workflow computes binned distributions of genotype-quality metrics across the carriers of each variant. For every configured FORMAT field it counts the genotypes whose value falls into each bin - optionally restricted to a variant filter and respecting whether larger or smaller values of that field are better - and can additionally bin allele-balance values. It emits a per-variant TSV of these distribution counts.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Array[String] gq_fields`: FORMAT fields whose values are binned, one per field to annotate.
- `Array[Array[Int]] gq_bins`: Bin edges for each field in `gq_fields`.
- `Array[String] gq_variant_filters`: Per-field expression restricting which variants the field is binned over.
- `Array[Boolean] gq_larger_field`: Per-field flag indicating whether larger values of the field are better.
- `Boolean ab_annotation`: Whether to additionally compute allele-balance distributions.
- `Array[Float] ab_bins`: Bin edges for allele-balance values.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File annotations_tsv_gq`: TSV of per-variant genotype-quality (and optional allele-balance) distribution counts.

### [AnnotateL1MEAID](../wdl/annotation/AnnotateL1MEAID.wdl)
This workflow first runs RepeatMasker on the insertions in an input VCF. It then uses its output to run L1ME-AID (https://github.com/Markloftus/L1ME-AID) and INTACT_MEI (https://github.com/xzhuo/INTACT_MEI) in order to identify, annotate and filter mobile element insertion (MEI) calls. It restricts to insertions at or above a minimum length and emits a TSV of the resulting MEI annotations.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int min_length`: Minimum insertion length to consider for MEI annotation.
- `String prefix`: Prefix for output file names.
- `String intact_mei_docker`, `String l1meaid_docker`, `String repeatmasker_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (9).

Outputs:
- `File annotations_tsv_l1meaid`: TSV of L1ME-AID and INTACT_MEI MEI annotations.

### [AnnotateMEDs](../wdl/annotation/AnnotateMEDs.wdl)
This workflow annotates mobile element deletions (MEDs) by intersecting the deletions in the input VCF against a catalog of known mobile-element loci. Deletions are extracted to BED form and matched to the catalog using size-similarity, reciprocal-overlap, breakpoint-window and sequence-similarity thresholds, producing a TSV of the deletions identified as MEDs.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int del_breakpoint_window`: Breakpoint window, in bp, for matching. (default `500`)
- `Float del_reciprocal_overlap`: Minimum reciprocal overlap for a deletion to match a catalog locus. (default `0.0`)
- `Float del_sequence_similarity`: Minimum sequence similarity for a deletion to match a catalog locus. (default `0.7`)
- `Float del_size_similarity`: Minimum size similarity for a deletion to match a catalog locus. (default `0.7`)
- `File mei_catalog`: From references.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File annotations_tsv_meds`: TSV of deletions identified as mobile element deletions.

### [AnnotateMEIs](../wdl/annotation/AnnotateMEIs.wdl)
This workflow consolidates the mobile element insertion calls produced by the `AnnotateL1MEAID`, `AnnotatePALMER` and `AnnotateSVAN` workflows into a single harmonized set. It reconciles the three per-tool annotation TSVs, using the SVAN annotation header for typing, to produce a final TSV of MEI calls.

Inputs:
- `File annotations_tsv_l1meaid`: MEI annotation TSV output by `AnnotateL1MEAID`.
- `File annotations_tsv_palmer`: MEI annotation TSV output by `AnnotatePALMER`.
- `File annotations_tsv_svan`: MEI annotation TSV output by `AnnotateSVAN`.
- `File annotations_header_svan`: Annotation header output by `AnnotateSVAN`, used to type the consolidated fields.
- `Array[String] contigs`: Contigs to annotate.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File annotations_tsv_meis`: Consolidated TSV of mobile element insertion calls.

### [AnnotatePALMER](../wdl/annotation/AnnotatePALMER.wdl)
This workflow leverages PALMER (https://github.com/WeichenZhou/PALMER) in order to annotate MEI calls for a cohort in a given cohort VCF. It retains the genotypes present in the VCF, simply adding an INFO field `ME_TYPE` to insertions whose characteristics match those of the PALMER calls. Matching is performed per MEI type using type-specific reciprocal-overlap, size-similarity, sequence-similarity, breakpoint-window and minimum-shared-sample thresholds.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File PALMER_vcf`: VCF of PALMER MEI calls to match against.
- `File PALMER_vcf_idx`: Index for `PALMER_vcf`.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Array[String] mei_types`: MEI types to run on - must be a subset of [`ALU`, `SVA`, `LINE` or `HERVK`].
- `Int min_length`: Minimum insertion length to consider for annotation.
- `File rm_out`: RepeatMasker output for the input VCF's insertions.
- `Int rm_buffer`: Padding, in bp, applied around RepeatMasker annotations when matching.
- `File ref_fai`: From references.
- `Int ins_breakpoint_window_alu`: Per-type breakpoint window, in bp, for matching. (default `200`)
- `Float ins_reciprocal_overlap_alu`: Per-type minimum reciprocal overlap for matching. (default `0.9`)
- `Float ins_sequence_similarity_alu`: Per-type minimum sequence similarity for matching. (default `0.9`)
- `Float ins_size_similarity_alu`: Per-type minimum size similarity for matching. (default `0.9`)
- `Int ins_min_shared_samples_alu`: Per-type minimum number of shared samples for matching. (default `0`)
- `Int ins_breakpoint_window_line`: Per-type breakpoint window, in bp, for matching. (default `200`)
- `Float ins_reciprocal_overlap_line`: Per-type minimum reciprocal overlap for matching. (default `0.9`)
- `Float ins_sequence_similarity_line`: Per-type minimum sequence similarity for matching. (default `0.9`)
- `Float ins_size_similarity_line`: Per-type minimum size similarity for matching. (default `0.9`)
- `Int ins_min_shared_samples_line`: Per-type minimum number of shared samples for matching. (default `0`)
- `Int ins_breakpoint_window_sva`: Per-type breakpoint window, in bp, for matching. (default `200`)
- `Float ins_reciprocal_overlap_sva`: Per-type minimum reciprocal overlap for matching. (default `0.9`)
- `Float ins_sequence_similarity_sva`: Per-type minimum sequence similarity for matching. (default `0.9`)
- `Float ins_size_similarity_sva`: Per-type minimum size similarity for matching. (default `0.9`)
- `Int ins_min_shared_samples_sva`: Per-type minimum number of shared samples for matching. (default `0`)
- `Int ins_breakpoint_window_hervk`: Per-type breakpoint window, in bp, for matching. (default `200`)
- `Float ins_reciprocal_overlap_hervk`: Per-type minimum reciprocal overlap for matching. (default `0.9`)
- `Float ins_sequence_similarity_hervk`: Per-type minimum sequence similarity for matching. (default `0.9`)
- `Float ins_size_similarity_hervk`: Per-type minimum size similarity for matching. (default `0.9`)
- `Int ins_min_shared_samples_hervk`: Per-type minimum number of shared samples for matching. (default `0`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotations_tsv_palmer`: TSV of insertions annotated with their PALMER `ME_TYPE`.

### [AnnotateRegion](../wdl/annotation/AnnotateRegion.wdl)
This workflow annotates each variant with the genomic region class it falls within - simple repeat (`SR`), segmental duplication (`SD`), RepeatMasker region (`RM`) or unique sequence (`US`) - by intersecting it against the corresponding BED panels. It emits a TSV of per-variant `REGION` assignments.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File simple_repeats_bed`: From references.
- `File seg_dup_bed`: From references.
- `File repeat_masked_bed`: From references.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotations_tsv_region`: TSV of per-variant genomic-region assignments.

### [AnnotateSQMetrics](../wdl/annotation/AnnotateSQMetrics.wdl)
This workflow recomputes site-level quality metrics for each variant directly from its genotype-level data. It clears any stale allele-specific INFO fields and recalculates Hardy-Weinberg equilibrium, the inbreeding coefficient, the maximum p(allele balance), and the allele-specific quality approximation, quality-by-depth and variant depth from the per-sample DP, PL and AD fields, emitting a per-variant TSV.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File annotations_tsv_sq`: TSV of recomputed site-level quality metrics.

### [AnnotateSVAN](../wdl/annotation/AnnotateSVAN.wdl)
This workflow leverages SVAN (https://github.com/REPBIO-LAB/SVAN) in order to annotate Mobile Element Insertions (MEIs), Mobile Element Deletions, Tandem Duplications, Dispersed Duplications and Nuclear Mitochondrial Segments (NUMT). It processes insertions and deletions separately, first running Tandem Repeat Finder (TRF) on the inserted or deleted sequence of each SV in the input VCF and then running SVAN over the result, before extracting and aligning the annotations into a single TSV. Before extraction, the `DUP_COORD` field produced by SVAN is reformatted: any `flank_`-prefixed relative coordinates are resolved to absolute genomic positions, preserving the original order of comma-separated values.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `String type_field`: INFO field giving each variant's allele type, used to select insertions/deletions for annotation. (default `allele_type`)
- `String type_ins`: Value of `type_field` identifying an insertion. (default `ins`)
- `String type_del`: Value of `type_field` identifying a deletion. (default `del`)
- `String length_field`: INFO field giving each variant's allele length. (default `allele_length`)
- `Int min_length`: Minimum insertion/deletion length to consider for annotation. (default `0`)
- `Boolean annotate_ins`: Whether to annotate insertions. (default `true`)
- `Boolean annotate_del`: Whether to annotate deletions. (default `true`)
- `File vntr_bed`: From references.
- `File exons_bed`: From references.
- `File repeats_bed`: From references.
- `File ref_fa`: From references.
- `Array[File] ref_fa_idx`: From references.
- `File mei_fa`: From references.
- `Array[File] mei_fa_idx`: From references.
- `String prefix`: Prefix for output file names.
- `String svan_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (19).

Outputs:
- `File annotations_tsv_svan`: TSV of SVAN annotations.
- `File annotations_header_svan`: Header lines describing the SVAN annotation fields.

### [AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl)
This workflow leverages SVAnnotate (https://gatk.broadinstitute.org/hc/en-us/articles/30332011989659-SVAnnotate) in order to annotate predicted functional effects for SVs. It conditionally only runs SVs through this workflow, ignoring all SNVs and indels, converting each SV to a symbolic representation before annotating it against coding and noncoding panels and extracting the resulting `PREDICTED_` annotations into a TSV.

Note: When converting to symbolic representation, all DUP allele types (including `dup_interspersed`, `inv_dup`, `complex_dup`, etc.) are treated as DUP.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `Int min_length`: Minimum length for a variant to be treated as an SV and annotated.
- `File coding_gtf`: From references.
- `File noncoding_bed`: From references.
- `String prefix`: Prefix for output file names.
- `String gatk_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (8).

Outputs:
- `File annotations_tsv_svannotate`: TSV of SVAnnotate functional-effect annotations.
- `File annotations_header_svannotate`: Header lines describing the SVAnnotate annotation fields.

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

### [AnnotateVEPHail](../wdl/annotation/AnnotateVEPHail.wdl)
This workflow leverages the Ensembl Variant Effect Predictor (VEP) (https://useast.ensembl.org/info/docs/tools/vep/index.html) in order to annotate predicted functional effects based on site-level information. It strips genotypes, scatters the VCF into shards, optionally normalizes and splits multiallelics around the VEP call, and uses Hail in order to run this annotation process in a more efficient and scalable manner before concatenating the per-shard annotations into a single TSV.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File ref_fa_gz`: bgzipped `ref_fa`, from references.
- `File ref_fai_gz`: Index for `ref_fa_gz`, from references.
- `File ref_vep_cache`: From references.
- `String? subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCF before annotation.
- `String split_vcf_hail_script`: Path to the Hail script used to scatter the VCF (defaults to this repository's copy on `main`). (default `https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/helper/split_vcf_hail.py`)
- `String vep_annotate_hail_python_script`: Path to the Hail script used to run VEP (defaults to this repository's copy on `main`). (default `https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/helper/vep_annotate_hail.py`)
- `String genome_build`: Genome build to annotate against. (default `GRCh38`)
- `String vep_json_schema`: Hail type schema describing the structure of VEP's JSON output. (default `Struct{allele_string:String,colocated_variants:Array[Struct{allele_string:String,clin_sig:Array[String],clin_sig_allele:String,end:Int32,id:String,phenotype_or_disease:Int32,pubmed:Array[Int32],somatic:Int32,start:Int32,strand:Int32}],context:String,end:Int32,id:String,input:String,intergenic_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],impact:String,minimised:Int32,variant_allele:String}],most_severe_consequence:String,motif_feature_consequences:Array[Struct{allele_num:Int32,consequence_terms:Array[String],high_inf_pos:String,impact:String,minimised:Int32,motif_feature_id:String,motif_name:String,motif_pos:Int32,motif_score_change:Float64,transcription_factors:Array[String],strand:Int32,variant_allele:String}],regulatory_feature_consequences:Array[Struct{allele_num:Int32,biotype:String,consequence_terms:Array[String],impact:String,minimised:Int32,regulatory_feature_id:String,variant_allele:String}],seq_region_name:String,start:Int32,strand:Int32,transcript_consequences:Array[Struct{allele_num:Int32,amino_acids:String,appris:String,biotype:String,canonical:Int32,ccds:String,cdna_start:Int32,cdna_end:Int32,cds_end:Int32,cds_start:Int32,codons:String,consequence_terms:Array[String],distance:Int32,domains:Array[Struct{db:String,name:String}],exon:String,flags:String,gene_id:String,gene_pheno:Int32,gene_symbol:String,gene_symbol_source:String,hgnc_id:String,hgvsc:String,hgvsp:String,hgvs_offset:Int32,impact:String,intron:String,lof:String,lof_flags:String,lof_filter:String,lof_info:String,mane_select:String,mane_plus_clinical:String,minimised:Int32,pick:Int32,mirna:Array[String],polyphen_prediction:String,polyphen_score:Float64,protein_end:Int32,protein_start:Int32,protein_id:String,sift_prediction:String,sift_score:Float64,source:String,strand:Int32,swissprot:String,transcript_id:String,trembl:String,tsl:Int32,uniparc:String,uniprot_isoform:Array[String],variant_allele:String}],variant_class:String}`)
- `String normalize_check_ref`: `bcftools norm` `--check-ref` mode used when normalizing. (default `w`)
- `Boolean normalize_vcf`: Whether to normalize and split multiallelics around the VEP call. (default `false`)
- `Boolean localize_vcf`: Whether to localize the VCF before annotation. (default `true`)
- `Boolean has_index`: Whether the input VCF is indexed. (default `true`)
- `Boolean get_chromosome_sizes`: Whether to compute chromosome sizes. (default `false`)
- `Boolean split_by_chromosome`: Whether to scatter the VCF by chromosome. (default `false`)
- `Boolean split_into_shards`: Whether to scatter the VCF into fixed-size shards. (default `false`)
- `String prefix`: Prefix for output file names.
- `String hail_docker`, `String sv_base_mini_docker`, `String utils_docker`, `String vep_hail_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (11).

Outputs:
- `File annotations_tsv_vep`: TSV of VEP functional-effect annotations.

### [AnnotateVRS](../wdl/annotation/AnnotateVRS.wdl)
This workflow annotates each variant with its GA4GH Variant Representation Specification (VRS) attributes using a seqrepo sequence repository. It runs `vrs-annotate` per contig to add the VRS INFO fields, then extracts them into an annotation TSV of five locating columns (CHROM, POS, REF, ALT, ID) followed by a column for each VRS field.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation.
- `File seqrepo_tar`: From references.
- `String prefix`: Prefix for output file names.
- `String utils_docker`, `String vrs_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File annotations_tsv_vrs`: TSV of per-variant VRS attributes (`VRS_Allele_IDs`, `VRS_Error`, `VRS_Starts`, `VRS_Ends`, `VRS_States`, `VRS_Lengths`, `VRS_RepeatSubunitLengths`).


## Annotation Utilities


### [AnnotateTREndTags](../wdl/annotation_utils/AnnotateTREndTags.wdl)
This utility adds an `END` INFO tag to the tandem-repeat records of a VCF, computed per contig, so that downstream tools correctly interpret the span of each TR call. It outputs the updated VCF.

Inputs:
- `File vcf`: VCF to update.
- `File vcf_idx`: Index for VCF to update.
- `Array[String] contigs`: Contigs to process within the input VCF.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File vcf_with_end`: VCF with `END` tags added to tandem-repeat records.
- `File vcf_with_end_idx`: Index for the updated VCF.

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
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File allele_type_annotated_vcf`: VCF annotated with `allele_type`.
- `File allele_type_annotated_vcf_idx`: Index for the annotated VCF.

### [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)
This utility applies a list of annotation TSVs to a VCF as new INFO fields, adding the specified field names, descriptions, types and numbers for each TSV in turn. Each annotation source can optionally have its TSV sorted, the VCF pre-subset and its TSV rows filtered beforehand. It outputs the annotated VCF.

Inputs:
- `File vcf`: VCF to annotate.
- `File vcf_idx`: Index for VCF to annotate.
- `Array[File] annotations_tsvs`: Annotation TSVs to apply, each as a set of INFO fields.
- `Array[String] contigs`: Contigs to annotate within the input VCF.
- `Int? records_per_shard`: Number of variants to keep within a single shard during annotation. When set, each contig VCF is sharded by record count, annotated in parallel and concatenated.
- `Array[Boolean] sort_tsvs`: Per-TSV flag indicating whether to sort the TSV before annotation. (default `[]`)
- `Array[Boolean] strip_info_fields_per_tsv`: Per-TSV flags selecting which annotations have their pre-existing INFO fields stripped before the annotation is written. Must be the same length as the annotation TSV list, and takes precedence over `strip_info_fields_by_name`. (default `[]`)
- `Array[String] strip_info_fields_by_name`: Names of INFO fields to strip from the input VCF before annotation. Ignored when `strip_info_fields_per_tsv` is set. (default `[]`)
- `Array[String] subset_vcf_strings`: Per-TSV `bcftools view` arguments used to pre-subset the VCF. (default `[]`)
- `Array[String] awk_tsv_conditions`: Per-TSV `awk` condition used to filter the TSV rows applied. (default `[]`)
- `Array[Array[Int]] subset_tsv_columns`: Per-TSV 1-based column indices (beyond the always-kept CHROM/POS/REF/ALT/ID columns 1-5) to keep from the TSV before annotation; an empty list for a TSV keeps all columns. (default `[]`)
- `Array[Array[String]] info_names`: INFO field names added by each annotation TSV.
- `Array[Array[String]] info_descriptions`: INFO field header descriptions for each annotation TSV.
- `Array[Array[String]] info_types`: INFO field types for each annotation TSV.
- `Array[Array[String]] info_numbers`: INFO field `Number` values for each annotation TSV.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File annotated_vcf`: Annotated VCF.
- `File annotated_vcf_idx`: Index for the annotated VCF.

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

### [CombineTRs](../wdl/annotation_utils/CombineTRs.wdl)
This utility combines tandem-repeat VCFs from multiple callers for one sample or a cohort into one VCF. It checks sample consistency, sets missing filters to pass, tags each caller's calls, assigns TR identifiers, deduplicates overlapping variants and priority-merges the callers per contig. It outputs the combined TR VCF.

Inputs:
- `Array[File] tr_vcfs`: Tandem-repeat VCFs to combine, one per caller.
- `Array[File] tr_vcf_idxs`: Indexes for `tr_vcfs`.
- `Array[String] tr_callers`: Caller name for each VCF in `tr_vcfs`, used for tagging and merge priority.
- `Array[String] contigs`: Contigs to process.
- `Array[String] sample_ids`: Sample IDs expected across all input VCFs.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (8).

Outputs:
- `File combined_tr_vcf`: Combined tandem-repeat VCF.
- `File combined_tr_vcf_idx`: Index for the combined VCF.

### [ConcatenateMosDepth](../wdl/annotation_utils/ConcatenateMosDepth.wdl)
This utility concatenates a sample's per-contig `MosDepth` per-base coverage BED files into a single indexed BED. It outputs the combined per-base coverage BED and its index.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-contig mosdepth per-base coverage BED files to concatenate.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File mosdepth_per_base_combined`: Combined per-base coverage BED.
- `File mosdepth_per_base_combined_idx`: Index for the combined BED.

### [ConcatenateVcfsAcrossContigs](../wdl/annotation_utils/ConcatenateVcfsAcrossContigs.wdl)
This utility concatenates exactly one VCF per contig into a single VCF, optionally dropping genotypes before concatenation. It validates that the VCF, index, and contig arrays are aligned and that no contig is duplicated.

Inputs:
- `Array[File] vcfs`: Per-contig VCFs to concatenate.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs corresponding to `vcfs`.
- `Boolean drop_genotypes`: Whether to strip genotypes before concatenation. (default `false`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File concat_vcf`: Combined VCF.
- `File concat_vcf_idx`: Index for the combined VCF.

### [ExtractDisparateTRLoci](../wdl/annotation_utils/ExtractDisparateTRLoci.wdl)
This utility subsets two VCFs to tandem-repeat variants (`INFO/allele_type=trv`) on one contig, then compares their loci. It produces one TSV for identities present in only one VCF, where identity is `CHROM`, `POS` and `len(REF)`, and another TSV for positive-base overlaps with distinct identities. Overlaps are identified with `bedtools intersect`; the overlapping TSV includes the `INFO/TRID` value from both VCFs.

Inputs:
- `File vcf_a`: First VCF to compare.
- `File vcf_a_idx`: Index for `vcf_a`.
- `File vcf_b`: Second VCF to compare.
- `File vcf_b_idx`: Index for `vcf_b`.
- `String contig`: Contig to compare within both VCFs.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File missing_variants_tsv`: Locus identities present in one VCF but missing from the other.
- `File overlapping_variants_tsv`: Overlapping locus pairs with distinct identities and their `TRID` values.

### [SummarizeAnnotations](../wdl/annotation_utils/SummarizeAnnotations.wdl)
This utility tallies annotation values across one or more VCFs to produce summary count tables, size-binned by allele class (SNV/DEL/INS/DUP/TRV). It always counts at the site level and can optionally count per sample, per allele, per functional gene consequence (from VEP/SVAnnotate `PREDICTED_*` fields), as raw per-variant value lists, and, when `create_plotting` is enabled, produce a separate set of AF-binned, region-aware Parquet tables for plotting (including a de novo transmission breakdown when a PED file is supplied and trios are found).

Inputs:
- `Array[File] vcfs`: VCFs whose annotations are counted.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[Int] length_bins_summary`: Size-bin edges used for the summary count tables. (default `[0, 1, 50, 500]`)
- `Array[Int] length_bins_plotting`: Size-bin edges used for the plotting tables. (default `[0, 1, 50, 100, 500, 5000, 50000]`)
- `Array[Float] af_bins_plotting`: Allele-frequency bin edges used for the plotting tables. (default `[0.0, 0.01, 0.05, 0.1, 0.5]`)
- `Boolean create_per_sample`: Whether to additionally produce per-sample counts. (default `false`)
- `Boolean create_per_allele`: Whether to additionally produce per-allele counts. (default `false`)
- `Boolean create_list`: Whether to additionally produce raw per-variant value-list tables. (default `false`)
- `Boolean create_functional`: Whether to additionally produce per-gene functional counts. (default `false`)
- `Boolean create_plotting`: Whether to additionally produce AF-binned Parquet tables for plotting. (default `false`)
- `Boolean use_ssd`: Whether to use SSD-backed local disks. (default `false`)
- `Boolean split_by_region`: Whether to split each VCF by genomic region before counting. (default `false`)
- `String subset_vcf_string`: `bcftools view` arguments used to pre-subset the VCFs. (default empty)
- `Int max_length`: Maximum variant length to count, or `-1` for no maximum. (default `-1`)
- `Int min_length`: Minimum variant length to count, or `-1` for no minimum. (default `-1`)
- `File? ped`: PED file used to identify trios for the de novo transmission breakdown (only used when `create_plotting` is enabled).
- `Int? records_per_shard`: Number of variants to keep within a single shard.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File summary_sites_tsv`: Site-level annotation counts.
- `File? summary_samples_tsv`: Per-sample counts (when `create_per_sample`).
- `File? summary_alleles_tsv`: Per-allele counts (when `create_per_allele`).
- `File? summary_list_tsv`: Raw per-variant value lists (when `create_list`).
- `File? summary_functional_tsv`: Per-gene functional counts (when `create_functional`).
- `File? summary_functional_samples_tsv`: Per-gene per-sample functional counts (when `create_functional` and `create_per_sample`).
- `File? summary_functional_alleles_tsv`: Per-gene per-allele functional counts (when `create_functional` and `create_per_allele`).
- `File? plotting_sites_parquet`: Site-level AF/size-binned counts, as Parquet (when `create_plotting`).
- `File? plotting_samples_parquet`: Per-sample AF/size-binned counts, as Parquet (when `create_plotting`).
- `File? plotting_alleles_parquet`: Per-allele AF/size-binned counts, as Parquet (when `create_plotting`).
- `File? plotting_denovo_parquet`: Per-proband de novo transmission counts, as Parquet (when `create_plotting` and trios are found via `ped`).
- `File? plotting_variant_list_parquet`: Raw per-variant genotype-count list, as Parquet (when `create_plotting`).

### [CreateCohortMethylationFile](../wdl/annotation_utils/CreateCohortMethylationFile.wdl)
This utility builds cohort-level CpG methylation matrices from per-sample `MethylationProfiling` BED outputs. For each contig, it merges every sample's combined and per-haplotype modification-score BEDs into a wide site-by-sample(/haplotype) matrix, filling `.` for sites missing in a given sample or haplotype. Samples can optionally be processed in shards (merged independently, then joined column-wise) to bound how many sample files are localized onto a single task at once.

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
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `Array[File] combined_methylation_beds`: Per-contig site-by-sample modification-score matrix BEDs.
- `Array[File] haplotype_methylation_beds`: Per-contig site-by-haplotype modification-score matrix BEDs.

### [CreateCohortCoverageSummary](../wdl/annotation_utils/CreateCohortCoverageSummary.wdl)
This utility builds a binned coverage matrix across a cohort from per-sample mosdepth BED outputs. It tiles the genome into windows, computes the mean coverage and threshold-crossing counts within each bin for every sample, and concatenates the results into a single coverage TSV.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-sample mosdepth coverage BED files.
- `Array[File] mosdepth_bed_idx`: Indexes for `mosdepth_bed_files`.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs over which to compute coverage.
- `Int window_size`: Size, in bp, of each genomic window.
- `Int bin_size`: Size, in bp, of each coverage bin within a window.
- `Array[Int] thresholds`: Coverage thresholds at which to count bins as covered.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File binned_coverage_tsv`: Binned coverage matrix across the cohort.

### [CreateDepthIntervals](../wdl/annotation_utils/CreateDepthIntervals.wdl)
This utility creates an interval file matching the fixed-width bins emitted by `CreateSampleReadCounts`. It writes 1-based, inclusive `contig:start-end` intervals in the requested contig order and omits each contig's trailing partial bin.

Inputs:
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs for which to create intervals, in output order.
- `Int bin_size`: Size, in bp, of each interval. Use the same value supplied to `CreateSampleReadCounts`.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File intervals`: Fixed-width interval file.

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

### [CreateCohortMetadata](../wdl/annotation_utils/CreateCohortMetadata.wdl)
This utility builds a cohort metadata file by combining a pedigree file with an ancestry-assignment file. It outputs the merged metadata file.

Inputs:
- `File ped_file`: Cohort pedigree file.
- `File ancestry_file`: Two-column file of sample IDs and ancestry labels.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File metadata`: Merged cohort metadata file.

### [CreateCohortPedigreeAncestryFilesAoUPhase2](../wdl/annotation_utils/CreateCohortPedigreeAncestryFilesAoUPhase2.wdl)
This utility builds the pedigree and two-column ancestry files consumed by `CreateCohortMetadata` from the All of Us Phase 2 ancestry-prediction table, keeping the `ancestry_pred` label for each requested sample. Every sample becomes a singleton PED row with unknown parents, sex and phenotype, since All of Us releases no pedigree or sex calls. The All of Us `eur` label is renamed to `nfe`, since that is the gnomAD label `compute_AFs.py` accepts; the remaining All of Us labels pass through unchanged. Samples the predictions do not cover are labelled `.`, the sentinel `compute_AFs.py` drops from its population list, so they contribute to the global `AC`/`AF`/`AN` but to no population-specific field. Both directions of mismatch between the cohort and the predictions are reported separately.

Inputs:
- `File ancestry_predictions`: All of Us Phase 2 ancestry predictions, keyed by `research_id`.
- `Array[String] sample_ids`: Sample IDs making up the cohort.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File ped`: Pedigree file covering every sample in `sample_ids`.
- `File ancestry`: Two-column file of sample IDs and ancestry labels.
- `File missing_samples`: Samples present in `ancestry_predictions` but not in `sample_ids`.
- `File samples_missing_ancestry`: Samples present in `sample_ids` but not in `ancestry_predictions`.

### [CreateSampleReadCounts](../wdl/annotation_utils/CreateSampleReadCounts.wdl)
This utility produces a binned read-counts file for a single sample from its per-contig mosdepth BED outputs, binning counts at a fixed resolution and merging across contigs. It outputs the binned read-counts file.

Inputs:
- `Array[File] mosdepth_bed_files`: Per-contig mosdepth coverage BED files for the sample.
- `Array[File] mosdepth_bed_idx`: Indexes for `mosdepth_bed_files`.
- `File ref_dict`: From references.
- `Array[String] contigs`: Contigs over which to bin read counts.
- `Int bin_size`: Size, in bp, of each read-count bin.
- `String sample_id`: ID of the sample being processed.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File binned_read_counts`: Binned read-counts file for the sample.

### [SummarizeSingletonCalls](../wdl/annotation_utils/SummarizeSingletonCalls.wdl)
This utility counts each sample's called genotypes across variant type, allele-length range, genomic region, evidence source, and sample-level alternate-allele count. It classifies calls supported only by `hapdiff` and/or `dipcall` as assemblies, calls with any other `EV` value as alignments, and calls without either kind of evidence as other. Output columns use the format `variant_type - size_range - region - count_type - singleton_type`.

Inputs:
- `Array[File] vcfs`: VCFs whose sample calls are counted.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Int min_length`: Minimum absolute `INFO/allele_length` to count. (default `50`)
- `Int? records_per_shard`: Number of variants to keep within a single shard.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File singleton_counts_tsv`: Wide per-sample count table.

### [SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)
This utility splits a VCF into one VCF per requested contig. It can add missing INFO-header lines, create genotype-free copies, rewrite SNV IDs, and rename source contigs to dbSNP or dbVar naming.

Inputs:
- `File vcf`: VCF to split.
- `File vcf_idx`: Index for `vcf`.
- `Array[String] contigs`: Contigs to extract.
- `Boolean create_no_geno`: Whether to also produce genotype-free VCFs. (default `false`)
- `Boolean modify_snv_ids`: Whether to rewrite SNV IDs. (default `false`)
- `Boolean rename_dbsnp_contigs`: Whether to rename source contigs from dbSNP naming. (default `false`)
- `Boolean rename_dbvar_contigs`: Whether to rename source contigs from dbVar naming. (default `false`)
- `Array[String]? missing_info_header_fields`: INFO-header lines to add if missing.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `Array[File] contig_vcfs`: Per-contig VCFs.
- `Array[File] contig_vcf_idxs`: Indexes for `contig_vcfs`.
- `Array[File] contig_no_geno_vcfs`: Per-contig genotype-free VCFs when requested.
- `Array[File] contig_no_geno_vcf_idxs`: Indexes for `contig_no_geno_vcfs` when requested.

### [StripGenotypes](../wdl/annotation_utils/StripGenotypes.wdl)
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

### [ExtractSampleVcfs](../wdl/annotation_utils/ExtractSampleVcfs.wdl)
This utility extracts per-sample VCFs from a cohort VCF, splitting each sample's variants into a SNV/indel VCF and an SV VCF based on a minimum SV length. It outputs the per-sample SNV/indel and SV VCFs.

Inputs:
- `File cohort_vcf`: Cohort VCF to extract from.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[String] contigs`: Contigs to process within the cohort VCF.
- `Array[String] sample_ids`: Samples to extract.
- `Int min_sv_length`: Minimum length at which a variant is routed to the SV VCF rather than the SNV/indel VCF.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `Array[File] snv_indel_vcfs`: Per-sample SNV/indel VCFs.
- `Array[File] snv_indel_vcf_idxs`: Indexes for the SNV/indel VCFs.
- `Array[File] sv_vcfs`: Per-sample SV VCFs.
- `Array[File] sv_vcf_idxs`: Indexes for the SV VCFs.

### [IdentifyLowCoverageRegions](../wdl/annotation_utils/IdentifyLowCoverageRegions.wdl)
This utility finds recurrent low-coverage regions from cohort mosdepth per-base BED files. It streams each sample independently, divides each chromosome into fixed bins anchored at position 0, and calculates each bin's base-weighted median coverage. Each sample's median binned coverage is rounded down before its regular low-coverage cutoff is calculated as `floor(median_coverage * median_coverage_cutoff)`; bins at or below this inclusive cutoff are flagged. A cohort bin fails when its low-coverage sample proportion is at or above the inclusive `sample_proportion_cutoff`.

Sample sex is read from a six-column PED (`1` for male and `2` for female). Male chrX and chrY use half the floored regular cutoff. Female chrY is excluded from binning, the sample median, sample histograms, and low-coverage calls. Samples missing from the PED or carrying unsupported sex codes fail explicitly. The workflow scatters once across samples, then aggregates flagged bins in one cohort task to avoid a nested scatter.

Inputs:
- `Array[File] mosdepth_bed_files`: One sorted, contiguous mosdepth per-base BED or BED.GZ per sample. Files must use four columns (`chrom`, `start`, `end`, `coverage`), cover each included chromosome from position 0, and have the same chromosome coordinate system.
- `Array[String] sample_ids`: Sample IDs corresponding by array index to `mosdepth_bed_files`.
- `File ped`: Six-column PED containing every input sample and its sex.
- `Int bin_size`: Fixed genomic bin size in bp.
- `Float median_coverage_cutoff`: Fraction of the floored sample median coverage defining the inclusive regular low-coverage cutoff. The operational cutoff is also rounded down: `floor(floor(median_coverage) * median_coverage_cutoff)`. For example, `0.2` gives a `12x` cutoff for a sample with median coverage `60x`. (default `0.2`)
- `Float sample_proportion_cutoff`: Inclusive minimum proportion of eligible samples with low coverage required to fail a cohort bin. (default `0.5`)
- `Float chrY_coverage_cutoff`: Fraction of a sample's autosomal median coverage that its chrY median coverage must reach for the sample to be inferred male. (default `0.2`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File sample_histograms_tar`: Tarball of per-sample coverage histograms with weighted, data-driven bins. The displayed range extends through at least the 95th percentile and the `Q3 + 3 * IQR` upper fence; omitted extreme bins are counted, actual low-coverage bins are orange, and the sample median plus regular/sex-chromosome cutoffs are annotated.
- `File chromosome_low_coverage_tar`: Tarball of per-chromosome plots showing number of samples flagged in each genomic bin.
- `File cohort_coverage_counts_tsv`: Raw, uncompressed TSV underlying the chromosome plots, with low-coverage and eligible-sample counts for every eligible genomic bin, including bins with zero low-coverage samples. chrY includes only male samples in its eligibility denominator and is omitted when no male samples are present.
- `File failed_bins_bed`: Raw, uncompressed, naturally chromosome-sorted BED of cohort bins whose low-coverage sample proportion is greater than or equal to `sample_proportion_cutoff`.
- `File sample_cutoffs_tsv`: TSV with `sample_id`, floored regular `cutoff`, and floored `median_coverage` for every input sample. Both values match those used for low-coverage calls and plots.

### [FilterLowCoverageRegions](../wdl/annotation_utils/FilterLowCoverageRegions.wdl)
This utility adds the `LOW_COVERAGE_REGION` FILTER to a VCF record when at least `min_region_coverage_cutoff` of its entire REF span overlaps a supplied low-coverage BED. The span is `POS-1` through `POS-1 + len(REF)` in 0-based half-open BED coordinates; therefore SNVs and insertions both have a one-base REF span at their VCF position. The decision is per record, independent of ALT count or content and `INFO/allele_type`. It adds the FILTER definition to the VCF header, preserves existing filters, and emits an indexed filtered VCF.

Inputs:
- `File vcf`: VCF to filter.
- `File vcf_idx`: Index for VCF to filter.
- `File low_coverage_regions_bed`: BED of low-coverage regions. It may include multiple contigs; the workflow subsets it to the requested `contigs`.
- `Array[String] contigs`: Contigs to process within the input VCF.
- `Float min_region_coverage_cutoff`: Inclusive minimum fraction of a variant's reference span that must overlap a low-coverage region to add the filter.
- `Int? records_per_shard`: Number of variants per shard. When set, variants are processed in parallel shards and concatenated.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File low_coverage_region_filtered_vcf`: VCF with low-coverage-region filters added.
- `File low_coverage_region_filtered_vcf_idx`: Index for the filtered VCF.

### [FilterLowCoverageGenotypes](../wdl/annotation_utils/FilterLowCoverageGenotypes.wdl)
This utility sets selected individual genotype calls to missing (`./.`) when `FORMAT/DP` is present and at or below that sample's low-coverage cutoff, while preserving every other FORMAT field. Male chrX/chrY calls use half the sample's cutoff, with sex read from a six-column PED. Every called genotype is eligible, whether homozygous reference or carrying an alternate allele, including partially called genotypes such as `./1`. Calls missing every GT allele or DP are left unchanged. Filtering can be restricted to variants whose INFO field matches a given value. It optionally shards by record count and outputs the filtered VCF plus a report of affected variants.

Inputs:
- `File vcf`: Cohort VCF to filter.
- `File vcf_idx`: Index for the cohort VCF.
- `File sample_cutoffs_tsv`: `sample_cutoffs_tsv` output from `IdentifyLowCoverageRegions`, containing `sample_id` and `cutoff` columns for every VCF sample.
- `File ped`: Six-column PED containing every VCF sample and its sex.
- `String? subset_unfilled_vcf_field`: INFO field used to limit which variants are filtered. Requires `subset_unfilled_vcf_value`.
- `String? subset_unfilled_vcf_value`: Value that `subset_unfilled_vcf_field` must equal for a variant to be filtered. Variants that don't match are left unfiltered.
- `Int? records_per_shard`: Number of variants per shard. When set, variants are processed in parallel shards and concatenated.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File filtered_vcf`: VCF with low-coverage genotypes set to missing.
- `File filtered_vcf_idx`: Index for `filtered_vcf`.
- `File filtered_genotypes_tsv`: TSV with one row per affected variant: `CHROM`, `POS`, `REF`, `ALT`, `ID`, pre- and post-filter allele counts, number of filtered samples, and comma-separated filtered sample IDs.

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

### [FillPhasedGenotypes](../wdl/annotation_utils/FillPhasedGenotypes.wdl)
This utility transfers phasing information from a phased VCF onto the genotypes of an unphased VCF over matching sites, optionally sharding each contig by region. It outputs the phased VCF.

Inputs:
- `File phased_vcf`: VCF providing the phasing information.
- `File phased_vcf_idx`: Index for `phased_vcf`.
- `File unphased_vcf`: VCF whose genotypes are phased.
- `File unphased_vcf_idx`: Index for `unphased_vcf`.
- `Array[String] contigs`: Contigs to process.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding each contig.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File hiphase_phased_vcf`: Phased VCF.
- `File hiphase_phased_vcf_idx`: Index for the phased VCF.

### [CreateTRGTHistograms](../wdl/annotation_utils/CreateTRGTHistograms.wdl)
This utility generates per-locus tandem repeat allele-frequency histograms, stratified by population and sex, from a multisample LPS (longest polymer sequence) table for use in the TR browser. It outputs a single combined histograms TSV.

Inputs:
- `File lps_tsv`: Multisample LPS table.
- `File metadata_tsv`: Sample metadata (population, sex) used to stratify the histograms.
- `Array[File] vcf_trid_metadata_tsvs`: Per-contig TRID metadata from `TRGTLPS.vcf_trid_metadata_tsvs`, index-aligned with the `contigs` input array. Required for any callset genotyped against a catalog that contains variation clusters, as it allows TRIDs that include several comma-separated LocusIds to be processed correctly. Defaults to empty, which keeps the previous behavior for catalogs of isolated repeats only.
- `Array[String] contigs`: Contigs to process within the LPS table.
- `Array[Array[String]] filter_trid_motif_pairs`: `(TRID, motif)` pairs to drop from `lps_tsv` before the histograms are computed, given as two-element arrays of the LPS table's `trid` and `motif` column values - e.g. `[['X-149631602-149631617-TCC,X-149631685-149631694-GCT,X-149631723-149631735-CGCCGT', 'CGC']]`. Use it for a row trgt-lps emitted from a spurious `INFO/MOTIFS` value, which cannot be resolved against `vcf_trid_metadata_tsvs` because no LocusId in the TRID carries that motif. Every pair must match at least one LPS row or the task fails. Pass an empty array to filter nothing.
- `String prefix`: Prefix for output file names.
- `String stranalysis_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File trgt_histograms_tsv`: Combined per-locus allele-frequency histograms TSV.

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
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (11).

Outputs:
- `File tr_annotated_vcf`: Base VCF annotated with integrated TR calls.
- `File tr_annotated_vcf_idx`: Index for the annotated VCF.

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

### [PreprocessVcfs](../wdl/annotation_utils/PreprocessVcfs.wdl)
This utility preprocesses and integrates one or more cohort VCFs into a single VCF. It first optionally converts symbolic alleles to sequence alleles, then applies any per-VCF sample-ID swaps, optionally subsets every VCF to the requested samples, and validates that the resulting sample sets are identical. Each VCF is then optionally normalized, annotated with core variant attributes and an optional source label, and length-filtered. Per-VCF controls are required arrays: an empty array disables that control for every VCF; a non-empty array must align with `vcfs`.

Inputs:
- `Array[File] vcfs`: Cohort VCFs to preprocess and merge.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[Boolean] normalize_vcfs`: Per-VCF normalization settings. `[]` disables normalization; otherwise aligned with `vcfs`.
- `Array[Boolean] convert_symbolic_to_sequence`: Per-VCF symbolic-allele conversion settings. `[]` disables conversion; otherwise aligned with `vcfs`. For enabled VCFs, `<DEL>` becomes a reference-anchored deletion, and `<DUP>` becomes a reference-anchored insertion using `SVLEN` or `END` as a fallback to determine the inserted-reference length. `<INS>` remains unchanged; existing `INFO/allele_length` and `INFO/allele_type` values are preserved, while missing values receive an absolute length from `SVLEN` or `END` as a fallback and `allele_type=ins`. `<INV>` remains symbolic but receives `INFO/allele_type=inv` and an absolute `INFO/allele_length` from `SVLEN` or `END` as a fallback. Any other angle-bracket symbolic ALT fails the workflow.
- `Array[String] source_tags`: Per-VCF `SOURCE` values. `[]` disables source tagging; otherwise aligned with `vcfs`. Required when at least one length cutoff is enabled.
- `Array[File] swap_sample_lists`: Per-VCF sample-ID swap maps, applied before sample subsetting. `[]` disables swapping; otherwise aligned with `vcfs`. A zero-byte map means no swap for that VCF.
- `Array[Int] min_length_cutoffs`: Per-VCF minimum absolute allele lengths. `[]` disables minimum-length filtering; otherwise aligned with `vcfs`. A value of `-1` disables this filter for that VCF. Calls with `abs(allele_length)` strictly below an enabled cutoff receive `SMALL_{source_tags[i]}`.
- `Array[Int] max_length_cutoffs`: Per-VCF maximum absolute allele lengths. `[]` disables maximum-length filtering; otherwise aligned with `vcfs`. A value of `-1` disables this filter for that VCF. Calls with `abs(allele_length)` strictly above an enabled cutoff receive `LARGE_{source_tags[i]}`.
- `Int? records_per_shard`: Number of variants to keep within a single shard during processing.
- `Array[String] sample_ids`: Samples to retain in every VCF. `[]` skips sample subsetting. After swaps and any subsetting, all input VCFs must contain identical sample sets.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (14).

Outputs:
- `File preprocessed_vcf`: Preprocessed and merged cohort VCF.
- `File preprocessed_vcf_idx`: Index for the preprocessed VCF.

### [FillBackbonePhasedGenotypes](../wdl/annotation_utils/FillBackbonePhasedGenotypes.wdl)
This utility merges a backbone-phased VCF with its no-TRGT counterpart (the same backbone-phasing run without TRGT calls included): for each still-unphased heterozygous genotype in `backbone_phased_vcf`, if a matching variant exists in `backbone_phased_notrgt_vcf` with a phased genotype, that phased `GT` (and `PS`) is pulled into the output. Region sharding is optional. It outputs the merged VCF and a per-sample TSV of heterozygous/unphased/pulled genotype counts.

Inputs:
- `File backbone_phased_vcf`: Backbone-phased VCF whose remaining unphased het genotypes are filled.
- `File backbone_phased_vcf_idx`: Index for `backbone_phased_vcf`.
- `File backbone_phased_notrgt_vcf`: Backbone-phased VCF (without TRGT calls) providing phased genotypes to pull from.
- `File backbone_phased_notrgt_vcf_idx`: Index for `backbone_phased_notrgt_vcf`.
- `String contig`: Contig to process.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the contig.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File backbone_merged_vcf`: Merged VCF with phased genotypes pulled in where available.
- `File backbone_merged_vcf_idx`: Index for the merged VCF.
- `File backbone_merged_tsv`: Per-sample TSV of heterozygous, unphased, and post-pull unphased genotype counts.

### [EvaluateBackbonePhasing](../wdl/annotation_utils/EvaluateBackbonePhasing.wdl)
This utility evaluates backbone-phasing accuracy by comparing backbone-phased VCFs against base (truth) VCFs. It assigns samples to their base VCFs, compares phased genotypes per contig and aggregates the results into tables broken down by variants outside tandem repeats, TR-enveloped variants and TR variants. It outputs these summary tables plus per-VCF status tables.

Inputs:
- `Array[File] backbone_phased_vcfs`: Backbone-phased VCFs to evaluate.
- `Array[File] backbone_phased_vcf_idxs`: Indexes for `backbone_phased_vcfs`.
- `Array[File] base_vcfs`: Base (truth) VCFs to compare against.
- `Array[File] base_vcf_idxs`: Indexes for `base_vcfs`.
- `Array[String] contigs`: Contigs to process.
- `Int max_variants`: Maximum number of variants to evaluate, or `-1` for no limit. (default `-1`)
- `Array[String]? subset_samples`: Samples to restrict the evaluation to.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (9).

Outputs:
- `File outside_tr_table`: Phasing-accuracy table for variants outside tandem repeats.
- `File tr_enveloped_table`: Phasing-accuracy table for TR-enveloped variants.
- `File trv_table`: Phasing-accuracy table for tandem-repeat variants.
- `File missing_samples`: Samples with no matching base VCF.
- `Array[File] vcf_tables`: Per-VCF variant-status tables.

### [EvaluateOverlappingTRLoci](../wdl/annotation_utils/EvaluateOverlappingTRLoci.wdl)
This utility evaluates how consistently TRGT genotypes the same reference bases when the repeat catalog defines overlapping loci. It takes a single-sample TRGT VCF, pairs every two records whose `POS`-to-`INFO/END` reference spans intersect and, for each pair, projects both haplotype sequences onto the shared reference interval and scores their agreement. Only pairs where at least one locus carries a non-reference call are evaluated, since two reference calls agree on the shared bases by construction. Each TRGT record is a full-locus replacement - `REF` is the reference sequence spanning the locus and each `ALT` is a complete haplotype sequence - so no reference FASTA or prior `bcftools norm` is needed. The projection aligns each haplotype to its own `REF` with `edlib` and assigns inserted bases to the reference base they follow, which means the projected sequence of a haplotype carrying a length change inherits the aligner's placement of that change; when a length change is placed at the edge of the shared interval its whole length is charged to that interval, so `max_similarity` can fall below `0`. Genotypes are unphased, so both haplotype assignments are scored and the one with the lower summed edit distance is reported. Haploid records are scored on their single haplotype and report `.` for the absent haplotype, and records with a missing genotype are skipped.

Inputs:
- `File vcf`: Single-sample TRGT VCF to evaluate.
- `File vcf_idx`: Index for `vcf`.
- `Array[String] contigs`: Contigs to process.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File overlapping_loci_tsv`: One row per evaluated locus pair, with the `INFO/TRID` and `INFO/MOTIFS` of both records reproduced verbatim, the shared reference interval and its length, the projected haplotype sequences of both records (`-` where the interval is deleted on that haplotype, `.` where the record is haploid), `min_edit_distance` - the edit distance summed over both haplotype pairs under the better of the two assignments - and `max_similarity`, that distance divided by the shared interval length times the number of compared haplotypes and subtracted from `1`, so `1` means the records agree exactly.

### [PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)
This utility bundles every genotype-update and post-processing step applied to a near-final callset into one workflow, with a required `run_` Boolean guarding each step so that the input VCF is left untouched when all are set to `false`. The per-record steps are applied in a single pass over the VCF: each variant is first matched against `transfer_vcf` and has its genotypes transferred (when `run_transfer_genotypes` is set) using its unmodified properties, after which the remaining steps - unphasing, ploidy normalization, TR-ID decrementing, MEI pruning, homopolymer flagging, singleton filtering and same-coordinate sorting - run in order. Some steps require an accompanying field - `run_transfer_genotypes` needs `transfer_vcf`, `run_unphase_samples` needs `unphase_samples`, and `run_normalize_ploidy` needs `ped`. The per-record pass can optionally be region-sharded via `shard_bin_size`.

Inputs:
- `File vcf`: VCF to post-process.
- `File vcf_idx`: Index for VCF to post-process.
- `Array[String] contigs`: Contigs to process within the input VCF.
- `Int? shard_bin_size`: Region-bin size, in bp, used when sharding the per-record pass.
- `Boolean run_clean_vcf_header`: Whether to run the header-cleaning step.
- `Boolean run_decrement_trv_ids`: Whether to decrement tandem-repeat variant IDs.
- `Boolean run_drop_filters`: Whether to drop the FILTER values listed in `drop_filters`.
- `Boolean run_filter_assembly_only_singletons`: Whether to apply the `ASSEMBLY_ONLY_SINGLETON` filter and emit a matching TSV.
- `Boolean run_filter_single_read_singletons`: Whether to apply the `SINGLE_READ_SUPPORT` filter to singleton calls.
- `Boolean run_flag_homopolymer_trvs`: Whether to flag tandem repeats with a length-1 shortest motif as `HOMOPOLYMER_TRV`.
- `Boolean run_normalize_ploidy`: Whether to normalize ploidy by sex - clearing chrY female calls, making chrX/chrY male calls hemizygous, enforcing diploidy and right-aligning unphased calls (requires `ped`).
- `Boolean run_prune_meis`: Whether to reclassify mobile elements whose length falls outside the expected bounds back to plain insertions/deletions.
- `Boolean run_reassign_suffixes`: Whether to run the variant-ID suffix reassignment step.
- `Boolean run_sorting`: Whether to sort records sharing a coordinate by absolute allele length and variant ID.
- `Boolean run_transfer_genotypes`: Whether to transfer genotypes from `transfer_vcf` onto heterozygous calls (run first; requires `transfer_vcf`).
- `Boolean run_unphase_samples`: Whether to unphase the samples in `unphase_samples` (requires `unphase_samples`).
- `Array[String] unphase_samples`: Samples to unphase when `run_unphase_samples` is set (defaults to empty). (default `[]`)
- `Array[String] drop_filters`: FILTER values to drop when `run_drop_filters` is set. (default `[]`)
- `File? transfer_vcf`: VCF whose genotypes are transferred when `run_transfer_genotypes` is set.
- `File? transfer_vcf_idx`: Index for `transfer_vcf`.
- `File? ped`: Cohort pedigree file, used for ploidy normalization when `run_normalize_ploidy` is set.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File post_processed_vcf`: Post-processed VCF.
- `File post_processed_vcf_idx`: Index for the post-processed VCF.
- `File? assembly_only_singletons_tsv`: Optional TSV containing one row per assembly-only singleton ALT allele; present only when `run_filter_assembly_only_singletons` is true.

### [ResolveHaplotypeOverlaps](../wdl/annotation_utils/ResolveHaplotypeOverlaps.wdl)
This utility detects and resolves haplotype-level overlaps among non-TR, non-TR-enveloped variants in a phased cohort VCF. For each sample, it extracts the sample's non-ref calls (excluding `allele_type='trv'` and `INFO/TR_ENVELOPED` variants), then sweeps each haplotype's variant intervals to find all overlapping pairs. Overlapping pairs are resolved by keeping the variant that spans more reference sequence (larger `len(REF)`) - which always favors DELs over INS or SNVs. When two variants span the same reference length, the higher-GQ call wins; remaining ties are broken by `INFO/allele_length`, then type rank (DEL > INS > SNV), then QUAL, then input-file order. The loser's FORMAT fields (`GT`, `GQ`, `DP`, `EV`, `BEV`, `AD`, `PL`) are cleared in the output VCF. The workflow scatters per-sample detection across all samples, then applies the collected clears to the given contig (with optional record-count sharding) to produce the resolved VCF.

Inputs:
- `File vcf`: Phased cohort VCF to resolve.
- `File vcf_idx`: Index for `vcf`.
- `String contig`: Contig to process.
- `Int? records_per_shard`: When set, shards the contig into chunks of this many records for the clearing step.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File overlap_resolved_vcf`: VCF with overlapping loser genotypes cleared.
- `File overlap_resolved_vcf_idx`: Index for `overlap_resolved_vcf`.
- `File overlap_tsv`: TSV of all detected overlap pairs, with columns `sample`, `haplotype`, `variant_id_retained`, `var_type_retained`, `size_bin_retained`, `variant_id_cleared`, `var_type_cleared`, `size_bin_cleared`.

### [SubsetTsvToColumns](../wdl/annotation_utils/SubsetTsvToColumns.wdl)
This utility subsets an annotation TSV to a chosen set of columns, optionally filtering rows to those whose columns match specified values. It outputs the subset TSV.

Inputs:
- `File annotations_tsv`: Annotation TSV to subset.
- `File annotations_header`: Header describing the TSV columns.
- `Array[String] subset_columns`: Columns to retain.
- `Array[Array[String]]? subset_column_values`: Per-column values to filter rows by.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File subset_tsv`: Column-subset TSV.

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

### [ConvertVcfToBed](../wdl/annotation_utils/ConvertVcfToBed.wdl)
This utility converts per-contig VCFs to one BED-like table with `svtk vcf2bed`. It can filter by variant length, rewrite selected INFO fields, convert records to insertion/deletion classes, shard large inputs, and control INFO, sample, filter, BND, CPX, compression, and output-extension behavior.

Inputs:
- `Array[File] vcfs`: Per-contig VCFs to convert.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs corresponding to `vcfs`.
- `Int? records_per_shard`: Optional record count per conversion shard.
- `Int? min_length`: Optional minimum absolute variant length.
- `String length_field`: INFO field containing variant length. (default `allele_length`)
- `Boolean convert_to_ins_del`: Whether to reduce variant classes to insertions/deletions. (default `false`)
- `Array[Array[String]]? switch_info_fields`: INFO-field/value rewrites applied before conversion.
- `Array[String] info_columns`: INFO columns included in output. (default `["ALL"]`)
- `Boolean include_samples`: Whether to include sample columns. (default `true`)
- `Boolean include_filters`: Whether to include FILTER. (default `true`)
- `Boolean split_bnd`: Whether to split BND records. (default `false`)
- `Boolean split_cpx`: Whether to split complex records. (default `false`)
- `Boolean output_gz`: Whether to gzip output. (default `false`)
- `Boolean output_bed`: Whether to use a `.bed` extension instead of `.tsv`. (default `false`)
- `String prefix`: Prefix for output file names.
- `String gatk_sv_lr_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File bed`: Combined converted BED/TSV artifact.

### [ExtractRandomCalls](../wdl/annotation_utils/ExtractRandomCalls.wdl)
This utility draws a random sample of variant/sample pairs from a set of VCF shards, for manual review or IGV curation. Each shard is sampled independently under the supplied filters, then the per-shard draws are pooled and down-sampled to the requested count using a fixed seed, so a given seed always yields the same selection.

Inputs:
- `Array[File] vcfs`: VCF shards to sample from.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Int count`: Number of variant/sample pairs to return.
- `Int random_seed`: Seed for the draw, so the selection is reproducible. (default `42`)
- `Float? min_af`: Restrict to variants within an allele-frequency range.
- `Float? max_af`: Restrict to variants within an allele-frequency range.
- `Int? min_ac`: Restrict to variants within an allele-count range.
- `Int? max_ac`: Restrict to variants within an allele-count range.
- `Boolean singleton`: Restrict to singletons. (default `false`)
- `Array[String] filters`: Restrict to variants carrying these FILTER values. (default `[]`)
- `Array[String] allele_types`: Restrict to these allele types. (default `[]`)
- `Int? min_allele_length`: Restrict to variants within an allele-length range.
- `Int? max_allele_length`: Restrict to variants within an allele-length range.
- `Array[String] include_samples`: Restrict the draw to these samples. (default `[]`)
- `Array[String] exclude_samples`: Exclude these samples from the draw. (default `[]`)
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File variant_sample_pairs`: TSV of the selected variant/sample pairs.
- `File candidate_summary`: TSV summarizing how many candidates each shard contributed.
- `File variant_vcf`: VCF containing just the selected variants.
- `File variant_vcf_idx`: Index for `variant_vcf`.


## Tools


### [Automop](../wdl/tools/Automop.wdl)
This tool runs `mop` (via FISS) to clean up unreferenced intermediate files in a Terra workspace, freeing storage. A dry-run mode reports what would be deleted without removing anything.

Inputs:
- `String workspace_namespace`: Terra workspace namespace to clean.
- `String workspace_name`: Terra workspace name to clean.
- `String user`: User running the cleanup.
- `Boolean dry_run`: Whether to report rather than perform deletions.
- `String prefix`: Prefix for output file names.
- `String automop_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File fissfc_log`: Log of the cleanup run.

### [BackbonePhase](../wdl/tools/BackbonePhase.wdl)
This tool transfers ('backbone') phasing from a set of base VCFs onto a target VCF for a single contig. It assigns each sample to its base VCF, computes the phase-flip orientation needed to make the target consistent with the backbone, and applies those flips. It outputs the phase-transferred VCF and a list of samples with no matching base VCF.

Inputs:
- `File vcf`: Target VCF to phase.
- `File vcf_idx`: Index for `vcf`.
- `Array[File] base_vcfs`: Base VCFs providing the backbone phasing.
- `Array[File] base_vcf_idxs`: Indexes for `base_vcfs`.
- `String contig`: Contig to process.
- `File? swap_samples_base`: Sample-ID swap map applied to the base VCFs.
- `Boolean allow_unphased_match_phase`: Whether to allow unphased genotypes to set the phase orientation. (default `false`)
- `String prefix`: Prefix for output file names.
- `String docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File transferred_vcf`: Phase-transferred VCF.
- `File transferred_vcf_idx`: Index for the transferred VCF.
- `File missing_samples`: Samples with no matching base VCF.

### [DeepVariant](../wdl/tools/DeepVariant.wdl)
This tool follows the DeepVariant path of the linked source workflow: it reads the source reference-bundle and small-variant-options JSON files, subsets the input BAM into its size-balanced shards, calls CPU or GPU DeepVariant, then merges the VCFs and gVCFs. The same region list controls BAM subsetting and DeepVariant `--regions`, preventing duplicate off-shard zero-depth gVCF blocks. It performs no explicit GCS-directory copy; map its declared outputs directly to attributes in the main Terra entity table.

Inputs:
- `File bam`: Aligned whole-genome BAM and index.
- `File bai`: Aligned whole-genome BAM and index.
- `String sex`: Biological sex; `M` enables the source reference bundle's haploid-contig and PAR settings.
- `String model_for_dv_andor_pepper`: DeepVariant model, such as `PACBIO` or `ONT_R104`.
- `File ref_bundle_json_file`: Source-compatible reference bundle JSON. Its size-balanced shard manifests are required.
- `File small_variant_calling_options_json`: Source-compatible small-variant options JSON, providing DeepVariant threads, memory, GPU use, and haploid contigs.
- `Array[String] gcp_zones`: Placement zones. (default `["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]`)
- `String prefix`: Prefix for output file names.
- `String deepvariant_docker`, `String deepvariant_gpu_docker`, `String utils_docker`, `String resource_visualization_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File gvcf`: Merged gVCF.
- `File gvcf_idx`: Index for gvcf.
- `File vcf`: Merged VCF.
- `File vcf_idx`: Index for vcf.
- `Array[File] resource_usage_logs`: Per-shard diagnostic outputs.
- `Array[File] resource_usage_visualizations`: Per-shard diagnostic outputs.
- `Array[File] visual_reports`: Per-shard diagnostic outputs.

### [Hifiasm](../wdl/tools/Hifiasm.wdl)
This tool assembles a sample's long reads into a haplotype-resolved de novo assembly using hifiasm (https://github.com/chhylp123/hifiasm). Reads are converted to FASTQ, assembled in bubble-phasing mode and the resulting assembly graphs are converted to bgzipped FASTA.

Without parental or Hi-C data the two haplotype assignments are arbitrary and switch between bubbles, so 'hap1' and 'hap2' do not correspond to the maternal and paternal haplotypes. Downstream callers that assume parental phase should not rely on which output a contig came from.

Inputs:
- `Array[File] bams`: Unaligned BAMs for the sample, one per SMRT cell.
- `String prefix`: Prefix for output file names.
- `String hifiasm_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File hifiasm_hap1_fa`: Bgzipped FASTA of the first haplotype assembly.
- `File hifiasm_hap2_fa`: Bgzipped FASTA of the second haplotype assembly.
- `File hifiasm_primary_fa`: Bgzipped FASTA of the primary contig assembly.
- `File hifiasm_hap1_gfa`: Assembly graph for the first haplotype assembly.
- `File hifiasm_hap2_gfa`: Assembly graph for the second haplotype assembly.
- `File hifiasm_primary_gfa`: Assembly graph for the primary contig assembly.
- `File hifiasm_log`: Console log from the hifiasm run, including the inferred coverage histogram.

### [HiFiCNV](../wdl/tools/HiFiCNV.wdl)
This tool runs PacBio HiFiCNV (https://github.com/PacificBiosciences/HiFiCNV) on a sample's aligned HiFi BAM to call copy number variants from read depth. It outputs the CNV VCF, a copy-number bedgraph, a depth BigWig track and the tool's log.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `String sex`: Sex of sample (one of `M` or `F`), used to select the matching expected-CN file.
- `File ref_fa`: Reference sequences FASTA file.
- `File ref_fai`: Index for `ref_fa`.
- `File exclude_bed`: Regions to exclude from CNV calling (e.g. centromeres).
- `File exclude_bed_idx`: Index for `exclude_bed`.
- `File expected_cn_male`: PAR regions and expected copy numbers for sex chromosomes, male.
- `File expected_cn_female`: PAR regions and expected copy numbers for sex chromosomes, female.
- `File? maf`: Optional minor-allele-frequency track passed to HiFiCNV as `--maf`.
- `String? cov_regex`: Optional regular expression passed to HiFiCNV as `--cov-regex`, selecting the contigs used to estimate expected coverage.
- `Boolean disable_vcf_filters`: Whether to pass `--disable-vcf-filters`, emitting every call rather than only those HiFiCNV would keep. (default `false`)
- `String prefix`: Prefix for output file names.
- `String hificnv_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File hificnv_vcf`: CNV calls VCF.
- `File hificnv_vcf_idx`: Index for the CNV calls VCF.
- `File hificnv_bedgraph`: Per-window copy number bedgraph.
- `File hificnv_depth_bw`: Depth BigWig track.
- `File hificnv_log`: HiFiCNV log file.

### [HiPhase](../wdl/tools/HiPhase.wdl)
This tool runs PacBio HiPhase (https://github.com/PacificBiosciences/HiPhase) to jointly phase a sample's small-variant, SV and (optionally) TRGT VCFs against its aligned reads. It preprocesses and synchronizes the input VCFs per contig, phases them together and optionally haplotags the BAM. It outputs the phased VCF, per-contig phasing statistics and an optional haplotagged BAM.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `File small_vcf`: Small-variant (SNV/indel) VCF to phase.
- `File small_vcf_idx`: Index for `small_vcf`.
- `File sv_vcf`: SV VCF to phase.
- `File sv_vcf_idx`: Index for `sv_vcf`.
- `File? trgt_vcf`: TRGT tandem-repeat VCF to additionally phase.
- `File? trgt_vcf_idx`: Index for `trgt_vcf`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to phase.
- `Int? trgt_min_repeat_unit`: Minimum repeat-unit length retained when filtering the TRGT VCF.
- `Boolean? trgt_normalize`: Whether to normalize the TRGT VCF before phasing.
- `Int? trgt_min_length_diff`: Minimum length difference retained when filtering the TRGT VCF.
- `Int? trgt_max_catalog_length`: Maximum catalog length retained when filtering the TRGT VCF.
- `String? hiphase_extra_args`: Additional arguments passed to HiPhase.
- `Boolean run_haplotagging`: Whether to also haplotag the BAM. (default `false`)
- `String prefix`: Prefix for output file names.
- `String hiphase_docker`, `String hiphase_preprocess_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (10).

Outputs:
- `File hiphase_vcf`: Phased VCF.
- `File hiphase_vcf_idx`: Index for the phased VCF.
- `Array[File] hiphase_haplotag_files`: Per-contig haplotag read assignments.
- `Array[File] hiphase_stats`: Per-contig phasing statistics.
- `Array[File] hiphase_blocks`: Per-contig phase blocks.
- `Array[File] hiphase_summary`: Per-contig phasing summaries.
- `File? hiphase_haplotagged_bam`: Haplotagged BAM (only when `run_haplotagging`).
- `File? hiphase_haplotagged_bam_idx`: Index for the haplotagged BAM (only when `run_haplotagging`).

### [MergeHiPhaseCallsets](../wdl/tools/MergeHiPhaseCallsets.wdl)
This tool merges per-sample HiPhase-phased VCFs into a cohort VCF on a per-contig basis, optionally also merging the TRGT tandem-repeat calls separately - fixing TRGT `END`/`AL` headers and propagating phase-set tags. It outputs the merged integrated VCF and an optional merged TRGT VCF.

Inputs:
- `Array[File] phased_vcfs`: Per-sample HiPhase-phased VCFs to merge.
- `Array[File] phased_vcf_idxs`: Indexes for `phased_vcfs`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to process.
- `Boolean merge_trgt`: Whether to additionally merge the TRGT tandem-repeat calls separately.
- `String merge_args`: Arguments controlling the VCF merge. (default `--merge id`)
- `String prefix`: Prefix for output file names.
- `String trgt_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (11).

Outputs:
- `File hiphase_merged_integrated_vcf`: Merged integrated cohort VCF.
- `File hiphase_merged_integrated_vcf_idx`: Index for the merged integrated VCF.
- `File? hiphase_merged_trgt_vcf`: Merged TRGT VCF (only when `merge_trgt`).
- `File? hiphase_merged_trgt_vcf_idx`: Index for the merged TRGT VCF (only when `merge_trgt`).

### [Kanpig](../wdl/tools/Kanpig.wdl)
This tool regenotypes a cohort SV VCF against each sample's aligned reads using Kanpig (https://github.com/ACEnglish/kanpig). It subsets the cohort to the target samples, runs Kanpig per sample with sex-aware ploidy beds, and merges the per-sample genotypes back into both a raw and a processed cohort VCF. It outputs the regenotyped (processed) and raw Kanpig VCFs.

Inputs:
- `File cohort_vcf`: Cohort SV VCF to regenotype.
- `File cohort_vcf_idx`: Index for the cohort VCF.
- `Array[File] bams`: Aligned reads, one per sample.
- `Array[File] bais`: Indexes for `bams`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File ploidy_bed_male`: From references.
- `File ploidy_bed_female`: From references.
- `Array[String] sample_ids`: Samples to regenotype.
- `Array[String] sexes`: Sex of each sample in `sample_ids`.
- `File? swap_samples`: Sample-ID swap map applied to the cohort VCF.
- `String merge_args`: Arguments controlling the per-sample genotype merge. (default `--merge id`)
- `String kanpig_params`: Parameters passed to Kanpig. (default `--neighdist 500 --gpenalty 0.04 --hapsim 0.97`)
- `String prefix`: Prefix for output file names.
- `String kanpig_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File sv_kanpig_vcf`: Regenotyped (processed) cohort VCF.
- `File sv_kanpig_vcf_idx`: Index for the processed VCF.
- `File sv_kanpig_raw_vcf`: Raw Kanpig cohort VCF.
- `File sv_kanpig_raw_vcf_idx`: Index for the raw VCF.

### [GLNexus](../wdl/tools/GLNexus.wdl)
This tool joint-calls per-sample gVCFs into a cohort VCF using GLnexus, then converts the result to a Hail MatrixTable. Calling is sharded over genomic ranges derived from the input gVCF names, and the per-range BCFs are concatenated back into a single VCF.

Inputs:
- `Array[File] gvcfs`: Per-sample gVCFs to joint-call.
- `Array[File] gvcf_idxs`: Indexes for `gvcfs`.
- `File ref_map_file`: Reference map describing the genome build.
- `Array[Array[File]]? background_sample_gvcfs`: Additional background-sample gVCFs to joint-call alongside the cohort.
- `Array[Array[File]]? background_sample_gvcf_idxs`: Indexes for `background_sample_gvcfs`.
- `Boolean force_add_missing_dp`: Add a `DP` FORMAT field to gVCFs that lack one before calling. (default `false`)
- `Boolean remove_duplicate_zero_depth_reference_blocks`: Drop duplicate zero-depth reference blocks before calling. (default `false`)
- `File? bed`: Restrict calling to these regions.
- `String config`: GLnexus preset configuration. (default `DeepVariantWGS`)
- `File? config_file`: Custom GLnexus configuration, used in place of `config`.
- `Boolean more_PL`: Emit additional PL values. (default `false`)
- `Boolean squeeze`: Squeeze the output representation. (default `false`)
- `Boolean trim_uncalled_alleles`: Remove alleles that no sample carries. (default `false`)
- `Int? num_cpus`: CPU count for the calling task; derived from the input count when unset.
- `Int max_cpus`: Upper bound on the derived CPU count. (default `64`)
- `String reference`: Reference genome build. (default `GRCh38`)
- `String? ref_fa`: Reference FASTA, used when registering a custom reference with Hail.
- `String? ref_fai`: Index for `ref_fa`.
- `String prefix`: Prefix for output file names.
- `String glnexus_docker`, `String hail_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `File joint_vcf`: Joint-called cohort VCF.
- `File joint_vcf_idx`: Index for `joint_vcf`.
- `File joint_mt`: Tarred Hail MatrixTable of the joint callset.

### [LongReadCNVs](../wdl/tools/LongReadCNVs.wdl)
This workflow calls cohort CNVs from long-read depth profiles with GATK gCNV, then converts, clusters and genotypes the depth calls. It outputs merged CNV calls, ploidy, and genotyped depth VCFs.

Inputs:
- `File intervals`: Interval list over which CNVs are called.
- `Array[String]+ sample_ids`: Sample IDs in the cohort.
- `Array[File]+ depth_profiles`: Per-sample read-depth profiles, aligned to `sample_ids`.
- `Array[String] contigs`: Contigs to process, given in reference dictionary order. `intervals` is subset to these before gCNV runs, and depth preprocessing, clustering and genotyping are likewise restricted, so genome-wide `intervals`, `merged_bincov`, `training_intervals` and `contig_ploidy_priors` are accepted.
- `Boolean sort_depth_profiles`: Whether to sort each depth profile by contig and position first, for profiles that are not already coordinate sorted.
- `String batch_id`: Identifier for the cohort batch.
- `File contig_ploidy_priors`: Contig ploidy priors used to determine per-sample contig ploidy.
- `File merged_bincov`: Merged read-depth evidence and its tabix index for depth genotyping.
- `File merged_bincov_idx`: Index for merged_bincov.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File ref_dict`: From references.
- `File pedigree`: Cohort pedigree used by depth genotyping.
- `File training_intervals`: Intervals used to train the depth genotyping model.
- `File median_coverage`: Per-sample median coverage table used by depth genotyping.
- `String variant_prefix`: Prefix used for generated variant IDs.
- `Int gcnv_qs_cutoff`: Minimum gCNV quality score for a segment to be kept. (default `30`)
- `Int num_intervals_per_scatter`: Number of intervals processed per scatter shard. (default `1500`)
- `String chr_x`: Name of the X contig in the reference. (default `chrX`)
- `String chr_y`: Name of the Y contig in the reference. (default `chrY`)
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
- `Float? ploidy_mean_bias_standard_deviation`: DetermineGermlineContigPloidy --mean-bias-standard-deviation.
- `Float? ploidy_mapping_error_rate`: DetermineGermlineContigPloidy --mapping-error-rate.
- `Float? ploidy_global_psi_scale`: DetermineGermlineContigPloidy --global-psi-scale.
- `Float? ploidy_sample_psi_scale`: DetermineGermlineContigPloidy --sample-psi-scale.
- `Float? gcnv_p_alt`: GermlineCNVCaller --p-alt.
- `Float? gcnv_p_active`: GermlineCNVCaller --p-active.
- `Float? gcnv_cnv_coherence_length`: GermlineCNVCaller --cnv-coherence-length.
- `Float? gcnv_class_coherence_length`: GermlineCNVCaller --class-coherence-length.
- `Int? gcnv_max_copy_number`: GermlineCNVCaller --max-copy-number.
- `Int? gcnv_max_bias_factors`: GermlineCNVCaller --max-bias-factors.
- `Float? gcnv_mapping_error_rate`: GermlineCNVCaller --mapping-error-rate.
- `Float? gcnv_interval_psi_scale`: GermlineCNVCaller --interval-psi-scale.
- `Float? gcnv_sample_psi_scale`: GermlineCNVCaller --sample-psi-scale.
- `Float? gcnv_depth_correction_tau`: GermlineCNVCaller --depth-correction-tau.
- `Float? gcnv_log_mean_bias_standard_deviation`: GermlineCNVCaller --log-mean-bias-standard-deviation.
- `Float? gcnv_init_ard_rel_unexplained_variance`: GermlineCNVCaller --init-ard-rel-unexplained-variance.
- `Int? gcnv_num_gc_bins`: GermlineCNVCaller --num-gc-bins.
- `Float? gcnv_gc_curve_standard_deviation`: GermlineCNVCaller --gc-curve-standard-deviation.
- `String? gcnv_copy_number_posterior_expectation_mode`: GermlineCNVCaller --copy-number-posterior-expectation-mode.
- `Boolean? gcnv_enable_bias_factors`: GermlineCNVCaller --enable-bias-factors.
- `Int? gcnv_active_class_padding_hybrid_mode`: GermlineCNVCaller --active-class-padding-hybrid-mode.
- `Float? gcnv_learning_rate`: GermlineCNVCaller --learning-rate.
- `Float? gcnv_adamax_beta_1`: GermlineCNVCaller --adamax-beta-1.
- `Float? gcnv_adamax_beta_2`: GermlineCNVCaller --adamax-beta-2.
- `Int? gcnv_log_emission_samples_per_round`: GermlineCNVCaller --log-emission-samples-per-round.
- `Float? gcnv_log_emission_sampling_median_rel_error`: GermlineCNVCaller --log-emission-sampling-median-rel-error.
- `Int? gcnv_log_emission_sampling_rounds`: GermlineCNVCaller --log-emission-sampling-rounds.
- `Int? gcnv_max_advi_iter_first_epoch`: GermlineCNVCaller --max-advi-iter-first-epoch.
- `Int? gcnv_max_advi_iter_subsequent_epochs`: GermlineCNVCaller --max-advi-iter-subsequent-epochs.
- `Int? gcnv_min_training_epochs`: GermlineCNVCaller --min-training-epochs.
- `Int? gcnv_max_training_epochs`: GermlineCNVCaller --max-training-epochs.
- `Float? gcnv_initial_temperature`: GermlineCNVCaller --initial-temperature.
- `Int? gcnv_num_thermal_advi_iters`: GermlineCNVCaller --num-thermal-advi-iters.
- `Int? gcnv_convergence_snr_averaging_window`: GermlineCNVCaller --convergence-snr-averaging-window.
- `Float? gcnv_convergence_snr_trigger_threshold`: GermlineCNVCaller --convergence-snr-trigger-threshold.
- `Int? gcnv_convergence_snr_countdown_window`: GermlineCNVCaller --convergence-snr-countdown-window.
- `Int? gcnv_max_calling_iters`: GermlineCNVCaller --max-calling-iters.
- `Float? gcnv_caller_update_convergence_threshold`: GermlineCNVCaller --caller-update-convergence-threshold.
- `Float? gcnv_caller_internal_admixing_rate`: GermlineCNVCaller --caller-internal-admixing-rate.
- `Float? gcnv_caller_external_admixing_rate`: GermlineCNVCaller --caller-external-admixing-rate.
- `Boolean? gcnv_disable_annealing`: GermlineCNVCaller --disable-annealing.
- `Int ref_copy_number_autosomal_contigs`: Reference copy number for autosomes. (default `2`)
- `Array[String]? allosomal_contigs`: Contigs treated as allosomal.
- `Int maximum_number_events_per_sample`: Maximum number of events permitted per sample. (default `1000`)
- `Float? defragment_max_dist`: Maximum gap, as a fraction of call length, across which adjacent calls are defragmented.
- `Boolean fast_mode`: Use SVCluster fast mode. (default `true`)
- `String clustering_algorithm`: SVCluster algorithm. (default `SINGLE_LINKAGE`)
- `Boolean? enable_cnv`: SVCluster behavior flags.
- `Boolean? default_no_call`: SVCluster behavior flags.
- `Boolean? omit_members`: SVCluster behavior flags.
- `String? breakpoint_summary_strategy`: SVCluster behavior flags.
- `Float? defrag_padding_fraction`: Defragmentation thresholds.
- `Float? defrag_sample_overlap`: Defragmentation thresholds.
- `Float depth_sample_overlap`: Required sample overlap for depth clustering. (default `0`)
- `Float depth_interval_overlap`: Required reciprocal interval overlap. (default `0.8`)
- `Float? depth_size_similarity`: Required size similarity.
- `Int depth_breakend_window`: Breakend join window in base pairs. (default `10000000`)
- `File? exclude_intervals`: Intervals whose overlapping calls are dropped.
- `Float exclude_overlap_fraction`: Overlap fraction at which a call is excluded. (default `0.5`)
- `File? gatk_to_svtk_script`: Override for the GATK-to-svtk conversion script.
- `Boolean svtk_set_pass`: Set FILTER to PASS during conversion. (default `false`)
- `String prefix`: Prefix for output file names.
- `String gatk_docker`, `String sv_base_mini_docker`, `String sv_pipeline_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (23).

Outputs:
- `File merged_cnvs_vcf`: Cohort CNV VCF after depth preprocessing.
- `File merged_cnvs_vcf_idx`: Index for `merged_cnvs_vcf`.
- `File ploidy_table`: Per-sample ploidy table.
- `File genotyped_depth_vcf`: Clustered CNV VCF genotyped from read depth.
- `File genotyped_depth_vcf_idx`: Index for `genotyped_depth_vcf`.
- `File genotyping_rd_table`: Read-depth evidence used for genotyping.

### [MethylationProfiling](../wdl/tools/MethylationProfiling.wdl)
This tool generates CpG methylation pileups from a haplotagged BAM using pb-CpG-tools (https://github.com/PacificBiosciences/pb-CpG-tools), producing combined and per-haplotype methylation BED tracks.

Inputs:
- `File bam`: Haplotagged aligned reads.
- `File bai`: Index for `bam`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String cpg_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File cpg_combined_bed`: Combined methylation pileup BED.
- `File cpg_combined_bed_idx`: Index for the combined BED.
- `File cpg_hap1_bed`: Haplotype 1 methylation pileup BED.
- `File cpg_hap1_bed_idx`: Index for the haplotype 1 BED.
- `File cpg_hap2_bed`: Haplotype 2 methylation pileup BED.
- `File cpg_hap2_bed_idx`: Index for the haplotype 2 BED.

### [MinimapAlignment](../wdl/tools/MinimapAlignment.wdl)
This workflow leverages Minimap2 (https://github.com/lh3/minimap2) in order to align a sample's maternal and paternal assemblies to a reference.

Inputs:
- `File assembly_mat`: Maternal assembly.
- `File assembly_pat`: Paternal assembly.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String sample_id`: ID of the sample being aligned.
- `String minimap_flags`: Parameters to use when running Minimap2. (default `-a -x asm20 --cs --eqx`)
- `Int minimap_threads`: Number of alignment threads. (default `32`)
- `String prefix`: Prefix for output file names.
- `String minimap_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File minimap_assembled_bam_mat`: Aligned maternal-assembly BAM.
- `File minimap_assembled_bai_mat`: Index for the maternal BAM.
- `File minimap_assembled_paf_mat`: Maternal-assembly PAF alignment.
- `File minimap_assembled_bam_pat`: Aligned paternal-assembly BAM.
- `File minimap_assembled_bai_pat`: Index for the paternal BAM.
- `File minimap_assembled_paf_pat`: Paternal-assembly PAF alignment.

### [MinimapReadAlignment](../wdl/tools/MinimapReadAlignment.wdl)
This tool aligns a sample's unaligned long reads to a reference using Minimap2 (https://github.com/lh3/minimap2). Every unaligned BAM for the sample is converted to FASTQ, streamed through Minimap2 in a single pass and coordinate-sorted into one indexed BAM.

Base modification tags are carried across from the unaligned BAM, since `samtools fastq` drops all tags by default and downstream methylation profiling needs them. Assemblies are aligned by `MinimapAlignment` instead.

Inputs:
- `Array[File] bams`: Unaligned BAMs for the sample, one per SMRT cell.
- `String sample_id`: ID of the sample being aligned, used for the read group ID and sample name.
- `String map_preset`: Minimap2 preset passed to '-x'. (default `map-hifi`)
- `Array[String] tags_to_preserve`: SAM tags carried over from the unaligned BAMs into the aligned BAM. (default `["MM", "ML"]`)
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String minimap2_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File aligned_bam`: Coordinate-sorted aligned reads.
- `File aligned_bai`: Index for the aligned reads.

### [MosDepth](../wdl/tools/MosDepth.wdl)
This tool runs mosdepth (https://github.com/brentp/mosdepth) to compute sequencing depth over a sample's BAM per contig. By default it emits per-base coverage; when `bin_size` is set, it instead windows depth into fixed-size bins (`--by`, `--no-per-base`) and emits per-region coverage.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for `bam`.
- `Array[String] contigs`: Contigs over which to compute depth.
- `Boolean single_contig`: Whether to run mosdepth once across all contigs instead of once per contig.
- `Boolean stream_mode`: Whether each per-contig run streams its region straight from the BAM rather than splitting the BAM by contig first. Ignored when `single_contig` is set.
- `Boolean fast_mode`: Use SVCluster fast mode.
- `Int? bin_size`: If set, windows depth into bins of this size (bp) and disables per-base output.
- `File? ref_fa`: Reference FASTA, required for CRAM input.
- `File? ref_fai`: Index for `ref_fa`.
- `String prefix`: Prefix for output file names.
- `String mosdepth_docker`, `String mosdepthstream_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `Array[File] mosdepth_dist`: Per-contig cumulative coverage distributions.
- `Array[File] mosdepth_summary`: Per-contig coverage summaries.
- `Array[File] mosdepth_per_base`: Per-contig per-base coverage (when `bin_size` is unset).
- `Array[File] mosdepth_per_base_csi`: Indexes for the per-base coverage.
- `Array[File] mosdepth_regions_bed`: Per-contig windowed coverage BEDs (when `bin_size` is set).
- `Array[File] mosdepth_regions_bed_csi`: Indexes for the windowed coverage BEDs.

### [PALMERAssembly](../wdl/tools/PALMERAssembly.wdl)
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
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to run PALMER on.
- `String sample`: ID of the sample being processed.
- `String mode`: PALMER run mode.
- `Array[String] mei_types`: MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `Array[String]? truvari_collapse_params`: Per-MEI-type Truvari parameters used when merging calls across haplotypes.
- `String prefix`: Prefix for output file names.
- `String palmer_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `Array[File] palmer_pat_calls`: Raw PALMER calls for the paternal haplotype, per MEI type.
- `Array[File] palmer_pat_tsd_reads`: PALMER TSD reads for the paternal haplotype, per MEI type.
- `Array[File] palmer_pat_vcfs`: Paternal-haplotype PALMER VCFs, per MEI type.
- `Array[File] palmer_pat_vcf_idxs`: Indexes for the paternal-haplotype VCFs.
- `Array[File] palmer_mat_calls`: Raw PALMER calls for the maternal haplotype, per MEI type.
- `Array[File] palmer_mat_tsd_reads`: PALMER TSD reads for the maternal haplotype, per MEI type.
- `Array[File] palmer_mat_vcfs`: Maternal-haplotype PALMER VCFs, per MEI type.
- `Array[File] palmer_mat_vcf_idxs`: Indexes for the maternal-haplotype VCFs.
- `Array[File] palmer_diploid_vcfs`: Diploid PALMER VCFs merged across haplotypes, per MEI type.
- `Array[File] palmer_diploid_vcf_idxs`: Indexes for the diploid VCFs.
- `File palmer_combined_vcf`: Final VCF combining all MEI types.
- `File palmer_combined_vcf_idx`: Index for the combined VCF.

### [PALMERDiploid](../wdl/tools/PALMERDiploid.wdl)
This tool runs PALMER (https://github.com/WeichenZhou/PALMER) on a single sample to generate mobile element insertion calls and convert them to a VCF. It shards the input BAM, runs PALMER per MEI type, merges the shard outputs and converts the raw calls into a per-type VCF, optionally bypassing execution when PALMER calls are supplied directly. It outputs the raw PALMER call and TSD files, per-type VCFs and a combined VCF.

Inputs:
- `File? bam`: Aligned reads to run PALMER on.
- `File? bai`: Index for `bam`.
- `Array[File]? override_palmer_calls`: Optional precomputed PALMER calls, causing the workflow to bypass execution.
- `Array[File]? override_palmer_tsd_files`: Optional precomputed PALMER TSD files, causing the workflow to bypass execution.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `Array[String] contigs`: Contigs to run PALMER on.
- `String sample`: ID of the sample being processed.
- `String mode`: PALMER run mode.
- `Array[String] mei_types`: MEI modes to run PALMER in - a subset of `ALU`, `SVA`, `LINE` or `HERVK`.
- `String prefix`: Prefix for output file names.
- `String palmer_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (5).

Outputs:
- `Array[File] palmer_calls`: Raw PALMER calls, per MEI type.
- `Array[File] palmer_tsd_reads`: PALMER TSD reads, per MEI type.
- `Array[File] palmer_diploid_vcfs`: Per-MEI-type PALMER VCFs.
- `Array[File] palmer_diploid_vcf_idxs`: Indexes for the per-type VCFs.
- `File palmer_combined_vcf`: Final VCF combining all MEI types.
- `File palmer_combined_vcf_idx`: Index for the combined VCF.

### [MergePALMERCallsets](../wdl/tools/MergePALMERCallsets.wdl)
This tool merges multiple PALMER MEI VCFs into a single VCF per contig and concatenates the result across contigs. It outputs the merged PALMER VCF.

Inputs:
- `Array[File] vcfs`: PALMER VCFs to merge.
- `Array[File] vcf_idxs`: Indexes for `vcfs`.
- `Array[String] contigs`: Contigs to process.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File palmer_merged_vcf`: Merged PALMER VCF.
- `File palmer_merged_vcf_idx`: Index for the merged VCF.

### [PAV](../wdl/tools/PAV.wdl)
This tool runs PAV (https://github.com/EichlerLab/pav) in batch mode across multiple samples' phased haplotype assemblies to call variants against the reference. It outputs per-sample VCFs, along with tarballs of the full PAV results and log directories.

Inputs:
- `Array[File] mat_haplotypes`: Maternal haplotype assemblies, one per sample.
- `Array[File] pat_haplotypes`: Paternal haplotype assemblies, one per sample.
- `Array[String] sample_ids`: Sample IDs, aligned by index to `mat_haplotypes`/`pat_haplotypes`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `String prefix`: Prefix for output file names.
- `String pav_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File pav_results_tarball`: Tarball of the full PAV results directory.
- `File pav_log_tarball`: Tarball of the full PAV log directory.
- `Array[File] pav_vcfs`: Per-sample called VCFs.
- `Array[File] pav_vcf_idx`: Indexes for the per-sample VCFs.
- `File? debug_sam`: Optional debug alignment file.
- `Array[File]? debug_temp`: Optional debug intermediate files.

### [PBSV](../wdl/tools/PBSV.wdl)
This tool calls structural variants from a sample's aligned long reads using pbsv (https://github.com/PacificBiosciences/pbsv). Signatures of structural variation are discovered from the alignments and then genotyped into a bgzipped, indexed VCF.

Supplying a tandem repeat BED lets pbsv collapse the alignment noise inside repeats, which reduces false calls at those loci.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for the aligned reads.
- `Boolean is_hifi`: Whether the reads are HiFi, which enables the pbsv optimisations for low-error reads. (default `true`)
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File? tandem_repeat_bed`: Tandem repeat intervals used to suppress alignment noise inside repeats.
- `String prefix`: Prefix for output file names.
- `String pbsv_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File pbsv_vcf`: Structural variant calls for the sample.
- `File pbsv_vcf_idx`: Index for the structural variant calls.
- `File pbsv_svsig`: Structural variant signatures discovered from the alignments.

### [RepeatMasker](../wdl/tools/RepeatMasker.wdl)
This workflow leverages RepeatMasker (https://github.com/Dfam-consortium/RepeatMasker) in order to annotate repeated and mobile-element content in the insertions of an input VCF. It extracts each insertion's inserted sequence to a FASTA, optionally restricted to a minimum length, and runs RepeatMasker over it.

Inputs:
- `File vcf`: VCF whose insertions are masked.
- `File vcf_idx`: Index for VCF.
- `Int? min_length`: Minimum insertion length to extract and mask.
- `String prefix`: Prefix for output file names.
- `String repeatmasker_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File rm_out`: RepeatMasker output table.
- `File rm_fa`: FASTA of the masked insertion sequences.

### [Sawfish](../wdl/tools/Sawfish.wdl)
This tool calls structural variants and copy-number variants from aligned long reads with sawfish. It runs the per-sample discover step in a scatter, then joint-calls across the cohort, using per-sample sex to select the expected copy-number track.

Inputs:
- `Array[File] bams`: Per-sample aligned reads.
- `Array[File] bais`: Indexes for `bams`.
- `Array[String] sexes`: Per-sample sex, aligned to `bams`, selecting the expected copy-number track.
- `Array[String] sample_ids`: Sample IDs, aligned to `bams`.
- `File ref_fa`: Reference FASTA.
- `File ref_fai`: Index for `ref_fa`.
- `File expected_cn_male`: Expected copy-number track for male samples.
- `File expected_cn_female`: Expected copy-number track for female samples.
- `File exclude_bed`: Regions excluded from calling.
- `File exclude_bed_idx`: Index for `exclude_bed`.
- `Int min_sv_size`: Minimum SV length to report. (default `35`)
- `Int min_sv_mapq`: Minimum mapping quality for supporting reads. (default `5`)
- `Boolean fast_cnv_mode`: Use the faster, less sensitive CNV mode. (default `false`)
- `Boolean disable_cnv`: Skip CNV calling entirely. (default `false`)
- `Boolean treat_single_copy_as_haploid`: Emit haploid genotypes on single-copy contigs. (default `false`)
- `Boolean report_supporting_reads`: Also emit the reads supporting each call. (default `false`)
- `String prefix`: Prefix for output file names.
- `String sawfish_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `File sawfish_vcf`: Joint-called SV and CNV VCF.
- `File sawfish_vcf_idx`: Index for `sawfish_vcf`.
- `Array[File] sawfish_bedgraphs`: Per-sample depth bedGraph files.
- `Array[File] sawfish_depth_bws`: Per-sample depth bigWig files.
- `File sawfish_log`: Joint-calling log.
- `File? sawfish_supporting_reads`: Supporting reads per call, emitted only when `report_supporting_reads` is set.

### [Sniffles](../wdl/tools/Sniffles.wdl)
This tool calls structural variants from a sample's aligned long reads using Sniffles2 (https://github.com/fritzsedlazeck/Sniffles). It emits both a bgzipped, indexed single-sample VCF and the sample's SNF file.

The SNF file holds the sample's raw structural variant candidates and is what Sniffles2 population mode re-genotypes across a cohort, so it is retained even though this pipeline merges callsets by other means.

Inputs:
- `File bam`: Aligned reads for the sample.
- `File bai`: Index for the aligned reads.
- `String sample_id`: ID of the sample being called, written to the VCF sample column.
- `Int min_sv_len`: Minimum structural variant length in base pairs to report. (default `50`)
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File? tandem_repeat_bed`: Tandem repeat intervals used to suppress alignment noise inside repeats.
- `String prefix`: Prefix for output file names.
- `String sniffles_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File sniffles_vcf`: Structural variant calls for the sample.
- `File sniffles_vcf_idx`: Index for the structural variant calls.
- `File sniffles_snf`: Structural variant candidates for the sample, for later population-mode calling.

### [TRGT](../wdl/tools/TRGT.wdl)
This workflow leverages TRGT (https://github.com/PacificBiosciences/trgt) in order to genotype short-tandem repeats.

Inputs:
- `File bam`: Aligned reads.
- `File bai`: Index for aligned reads.
- `String sample_id`: ID of the sample being genotyped.
- `String sex`: Sex of sample (one of `M` or `F`).
- `String catalog_name`: Name of the repeat catalog used, included in the output VCF filename.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File repeat_catalog_trgt`: From references.
- `String prefix`: Prefix for output file names.
- `String trgt_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides.

Outputs:
- `File trgt_vcf`: TRGT tandem-repeat genotype VCF.
- `File trgt_vcf_idx`: Index for the TRGT VCF.

### [TRGTLPS](../wdl/tools/TRGTLPS.wdl)
This tool runs the trgt-lps (https://github.com/PacificBiosciences/trgt-lps) tool per contig to compute the longest polymer sequence (LPS) within each TRGT-genotyped tandem-repeat locus for every sample, concatenating the results into a single TSV. Alongside each contig's LPS table it extracts a small TRID-metadata TSV from the same subset VCF, mapping every `(TRID, motif)` to the LocusIds that record covers, which `CreateTRGTHistograms` needs to resolve variation-cluster records whose TRID names several loci. It outputs the LPS TSV and the per-contig TRID-metadata TSVs.

Inputs:
- `File vcf`: TRGT VCF to process.
- `File vcf_idx`: Index for the TRGT VCF.
- `Array[String] contigs`: Contigs to process.
- `Boolean normalize_chry_haploid_genotypes`: Whether to normalize haploid chrY genotypes before computing the longest polymer sequence. Applied to the chrY shard only.
- `Array[Array[String]] filter_trid_motif_pairs`: `(TRID, motif)` pairs to drop from the concatenated LPS table, given as two-element arrays of the LPS table's `trid` and `motif` column values - e.g. `[['X-149631602-149631617-TCC,X-149631685-149631694-GCT,X-149631723-149631735-CGCCGT', 'CGC']]`. Use it for a row trgt-lps emitted from a spurious `INFO/MOTIFS` value, which `CreateTRGTHistograms` cannot resolve because no LocusId in the TRID carries that motif. Every pair must match at least one LPS row or the task fails. Pass an empty array to filter nothing.
- `String prefix`: Prefix for output file names.
- `String trgt_lps_docker`, `String utils_docker`, `String stranalysis_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (7).

Outputs:
- `File trgt_lps_tsv`: TSV of per-locus longest polymer sequences.
- `Array[File] vcf_trid_metadata_tsvs`: Per-contig TRID-metadata TSVs, index-aligned with the `contigs` input array, to be passed straight to `CreateTRGTHistograms`.

### [Vamos](../wdl/tools/Vamos.wdl)
This tool runs Vamos (https://github.com/ChaissonLab/vamos) in order to genotype tandem repeats against a Vamos repeat catalog, in read mode (from an aligned read BAM) and/or assembly mode (from per-haplotype assembly BAMs). It outputs the resulting Vamos VCFs.

Inputs:
- `File? read_bam`: Aligned reads to genotype in read mode.
- `File? read_bai`: Index for `read_bam`.
- `Array[File]? assembly_bams`: Per-haplotype assembly BAMs to genotype in assembly mode.
- `Array[File]? assembly_bais`: Indexes for `assembly_bams`.
- `File repeat_catalog_vamos`: Vamos repeat catalog to genotype against.
- `String sample_id`: ID of the sample being genotyped.
- `String prefix`: Prefix for output file names.
- `String vamos_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `Array[File] vamos_assembly_vcfs`: Per-haplotype assembly-mode Vamos VCFs.
- `Array[File] vamos_assembly_vcf_idxs`: Indexes for the assembly-mode VCFs.
- `File? vamos_reads_vcf`: Read-mode Vamos VCF.
- `File? vamos_reads_vcf_idx`: Index for the read-mode VCF.


## Sub-workflows
These workflows live in `wdl/utils/` and are building blocks rather than entry points. They are never registered in `.dockstore.yml` and are not run directly; a workflow imports one with `import "../utils/<Name>.wdl"` and calls it as `<Name>.<Name>`.


### Depth-based CNV pipeline
`LRCNVs`, `DepthPreprocessing`, `DepthClustering` and `GenotypeDepth` are called in that order by `LongReadCNVs`, which supplies their shared inputs.

### [LRCNVs](../wdl/utils/LRCNVs.wdl)
This component calls copy-number variants across a cohort using GATK germline CNV (gCNV) cohort mode. From per-sample depth profiles over a shared interval list it annotates and filters intervals, determines contig ploidy, fits gCNV across scattered interval shards, post-processes per-sample calls into genotyped interval and segment VCFs, and collects sample- and model-level QC.

Inputs:
- `File intervals`: Interval list over which CNVs are called.
- `Array[String]+ sample_ids`: Sample IDs in the cohort.
- `Array[File]+ depth_profiles`: Per-sample read-depth profiles, aligned to `sample_ids`.
- `Array[String] contigs`: Contigs to call CNVs on. `intervals` is subset to these before annotation and filtering, so a genome-wide interval list may be supplied. GATK copy-number tools reject a non-UNION `--interval-set-rule`, so the subset is taken as a separate step rather than by interval intersection.
- `String cohort_id`: Identifier for the cohort.
- `File contig_ploidy_priors`: Contig ploidy priors used to determine per-sample contig ploidy. May cover more contigs than `contigs`.
- `File ref_fa`: From references.
- `File ref_fai`: From references.
- `File ref_dict`: From references.
- `Int num_intervals_per_scatter`: Number of intervals processed per scatter shard.
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
- `Float? ploidy_mean_bias_standard_deviation`: DetermineGermlineContigPloidy --mean-bias-standard-deviation.
- `Float? ploidy_mapping_error_rate`: DetermineGermlineContigPloidy --mapping-error-rate.
- `Float? ploidy_global_psi_scale`: DetermineGermlineContigPloidy --global-psi-scale.
- `Float? ploidy_sample_psi_scale`: DetermineGermlineContigPloidy --sample-psi-scale.
- `Float? gcnv_p_alt`: GermlineCNVCaller --p-alt.
- `Float? gcnv_p_active`: GermlineCNVCaller --p-active.
- `Float? gcnv_cnv_coherence_length`: GermlineCNVCaller --cnv-coherence-length.
- `Float? gcnv_class_coherence_length`: GermlineCNVCaller --class-coherence-length.
- `Int? gcnv_max_copy_number`: GermlineCNVCaller --max-copy-number.
- `Int? gcnv_max_bias_factors`: GermlineCNVCaller --max-bias-factors.
- `Float? gcnv_mapping_error_rate`: GermlineCNVCaller --mapping-error-rate.
- `Float? gcnv_interval_psi_scale`: GermlineCNVCaller --interval-psi-scale.
- `Float? gcnv_sample_psi_scale`: GermlineCNVCaller --sample-psi-scale.
- `Float? gcnv_depth_correction_tau`: GermlineCNVCaller --depth-correction-tau.
- `Float? gcnv_log_mean_bias_standard_deviation`: GermlineCNVCaller --log-mean-bias-standard-deviation.
- `Float? gcnv_init_ard_rel_unexplained_variance`: GermlineCNVCaller --init-ard-rel-unexplained-variance.
- `Int? gcnv_num_gc_bins`: GermlineCNVCaller --num-gc-bins.
- `Float? gcnv_gc_curve_standard_deviation`: GermlineCNVCaller --gc-curve-standard-deviation.
- `String? gcnv_copy_number_posterior_expectation_mode`: GermlineCNVCaller --copy-number-posterior-expectation-mode.
- `Boolean? gcnv_enable_bias_factors`: GermlineCNVCaller --enable-bias-factors.
- `Int? gcnv_active_class_padding_hybrid_mode`: GermlineCNVCaller --active-class-padding-hybrid-mode.
- `Float? gcnv_learning_rate`: GermlineCNVCaller --learning-rate.
- `Float? gcnv_adamax_beta_1`: GermlineCNVCaller --adamax-beta-1.
- `Float? gcnv_adamax_beta_2`: GermlineCNVCaller --adamax-beta-2.
- `Int? gcnv_log_emission_samples_per_round`: GermlineCNVCaller --log-emission-samples-per-round.
- `Float? gcnv_log_emission_sampling_median_rel_error`: GermlineCNVCaller --log-emission-sampling-median-rel-error.
- `Int? gcnv_log_emission_sampling_rounds`: GermlineCNVCaller --log-emission-sampling-rounds.
- `Int? gcnv_max_advi_iter_first_epoch`: GermlineCNVCaller --max-advi-iter-first-epoch.
- `Int? gcnv_max_advi_iter_subsequent_epochs`: GermlineCNVCaller --max-advi-iter-subsequent-epochs.
- `Int? gcnv_min_training_epochs`: GermlineCNVCaller --min-training-epochs.
- `Int? gcnv_max_training_epochs`: GermlineCNVCaller --max-training-epochs.
- `Float? gcnv_initial_temperature`: GermlineCNVCaller --initial-temperature.
- `Int? gcnv_num_thermal_advi_iters`: GermlineCNVCaller --num-thermal-advi-iters.
- `Int? gcnv_convergence_snr_averaging_window`: GermlineCNVCaller --convergence-snr-averaging-window.
- `Float? gcnv_convergence_snr_trigger_threshold`: GermlineCNVCaller --convergence-snr-trigger-threshold.
- `Int? gcnv_convergence_snr_countdown_window`: GermlineCNVCaller --convergence-snr-countdown-window.
- `Int? gcnv_max_calling_iters`: GermlineCNVCaller --max-calling-iters.
- `Float? gcnv_caller_update_convergence_threshold`: GermlineCNVCaller --caller-update-convergence-threshold.
- `Float? gcnv_caller_internal_admixing_rate`: GermlineCNVCaller --caller-internal-admixing-rate.
- `Float? gcnv_caller_external_admixing_rate`: GermlineCNVCaller --caller-external-admixing-rate.
- `Boolean? gcnv_disable_annealing`: GermlineCNVCaller --disable-annealing.
- `Int ref_copy_number_autosomal_contigs`: Reference copy number for autosomes. (default `2`)
- `Array[String]? allosomal_contigs`: Contigs treated as allosomal.
- `Int maximum_number_events_per_sample`: Maximum number of events permitted per sample. (default `1000`)
- `String prefix`: Prefix for output file names.
- `String gatk_docker`, `String sv_base_mini_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (9).

Outputs:
- `File annotated_intervals`: Intervals annotated with GC content and tracks.
- `File filtered_intervals`: Intervals retained after filtering.
- `File contig_ploidy_model_tar`: Fitted contig-ploidy model.
- `File contig_ploidy_calls_tar`: Per-sample contig-ploidy calls.
- `Array[File] gcnv_model_tars`: Fitted gCNV models, one per scatter shard.
- `Array[Array[File]] gcnv_calls_tars`: Per-shard per-sample gCNV calls.
- `Array[File] gcnv_tracking_tars`: Per-shard model-fitting tracking files.
- `Array[File] genotyped_intervals_vcfs`: Per-sample genotyped interval VCFs.
- `Array[File] genotyped_intervals_vcf_idxs`: Indexes for `genotyped_intervals_vcfs`.
- `Array[File] genotyped_segments_vcfs`: Per-sample genotyped segment VCFs.
- `Array[File] genotyped_segments_vcf_idxs`: Indexes for `genotyped_segments_vcfs`.
- `Array[File] sample_qc_status_files`: Per-sample QC status files.
- `Array[String] sample_qc_status_strings`: Per-sample QC status strings.
- `File model_qc_status_file`: Model-level QC status file.
- `String model_qc_string`: Model-level QC status string.
- `Array[File] denoised_copy_ratios`: Per-sample denoised copy ratios.

### [DepthPreprocessing](../wdl/utils/DepthPreprocessing.wdl)
This sub-workflow converts per-sample gCNV genotyped-segment VCFs into cohort-level deletion and duplication call sets. Each sample's segments are converted to BED, merged per sample and then across the cohort separately for DEL and DUP, and finally rewritten as a single VCF alongside a ploidy table for downstream genotyping.

Inputs:
- `Array[String]+ sample_ids`: Sample IDs in the cohort.
- `Array[File]+ genotyped_segments_vcfs`: Per-sample gCNV genotyped-segment VCFs.
- `Array[File]+ genotyped_segments_vcf_idxs`: Indexes for `genotyped_segments_vcfs`.
- `File contig_ploidy_calls_tar`: Tarred gCNV contig-ploidy calls.
- `Array[String] contigs`: Contigs to process, given in reference dictionary order.
- `File ref_fai`: Reference FASTA index, used for contig ordering.
- `File pedigree`: Pedigree supplying per-sample sex.
- `String batch_id`: Identifier for the batch.
- `String? chr_x`: Allosome contig names, when they differ from the defaults.
- `String? chr_y`: Allosome contig names, when they differ from the defaults.
- `Int gcnv_qs_cutoff`: Minimum gCNV quality score for a segment to be kept.
- `Float? defragment_max_dist`: Maximum gap, as a fraction of call length, across which adjacent calls are defragmented.
- `String prefix`: Prefix for output file names.
- `String sv_base_mini_docker`, `String sv_pipeline_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (6).

Outputs:
- `File del_bed`: Cohort-merged deletion calls.
- `File del_bed_idx`: Index for del_bed.
- `File dup_bed`: Cohort-merged duplication calls.
- `File dup_bed_idx`: Index for dup_bed.
- `File merged_vcf`: Combined depth-based CNV VCF.
- `File merged_vcf_idx`: Index for merged_vcf.
- `File ploidy_table`: Per-sample, per-contig ploidy table consumed by `DepthClustering` and `GenotypeDepth`.

### [DepthClustering](../wdl/utils/DepthClustering.wdl)
This sub-workflow clusters the depth-based CNV calls across samples with GATK `SVCluster`, contig by contig, then optionally drops calls overlapping excluded intervals and converts the GATK representation back to svtk-style VCF before concatenating the per-contig results.

Inputs:
- `File depth_vcf`: Depth CNV VCF from `DepthPreprocessing`.
- `File depth_vcf_idx`: Index for depth_vcf.
- `File ploidy_table`: Ploidy table from `DepthPreprocessing`.
- `String variant_prefix`: Prefix applied to generated variant IDs.
- `Array[String] contigs`: Contigs to cluster over, given in reference dictionary order. These also become the `##contig` lines of the svtk-formatted output.
- `File ref_fa`: Reference FASTA, index and sequence dictionary.
- `File ref_fai`: Reference FASTA, index and sequence dictionary.
- `File ref_dict`: Reference FASTA, index and sequence dictionary.
- `Boolean fast_mode`: Use SVCluster fast mode. (default `true`)
- `String clustering_algorithm`: SVCluster algorithm. (default `SINGLE_LINKAGE`)
- `Boolean? enable_cnv`: SVCluster behavior flags.
- `Boolean? default_no_call`: SVCluster behavior flags.
- `Boolean? omit_members`: SVCluster behavior flags.
- `String? breakpoint_summary_strategy`: SVCluster behavior flags.
- `Float? defrag_padding_fraction`: Defragmentation thresholds.
- `Float? defrag_sample_overlap`: Defragmentation thresholds.
- `Float depth_sample_overlap`: Required sample overlap for depth clustering. (default `0`)
- `Float depth_interval_overlap`: Required reciprocal interval overlap. (default `0.8`)
- `Float? depth_size_similarity`: Required size similarity.
- `Int depth_breakend_window`: Breakend join window in base pairs. (default `10000000`)
- `File? exclude_intervals`: Intervals whose overlapping calls are dropped.
- `Float exclude_overlap_fraction`: Overlap fraction at which a call is excluded. (default `0.5`)
- `File? gatk_to_svtk_script`: Override for the GATK-to-svtk conversion script.
- `Boolean svtk_set_pass`: Set FILTER to PASS during conversion. (default `false`)
- `String prefix`: Prefix for output file names.
- `String gatk_docker`, `String sv_base_mini_docker`, `String sv_pipeline_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (4).

Outputs:
- `File clustered_vcf`: Cohort-clustered depth CNV VCF.
- `File clustered_vcf_idx`: Index for `clustered_vcf`.

### [GenotypeDepth](../wdl/utils/GenotypeDepth.wdl)
This sub-workflow trains a depth genotyping model on a set of training intervals, then genotypes the clustered depth CNV calls per contig with GATK and concatenates the results.

Inputs:
- `File vcf`: Clustered depth CNV VCF from `DepthClustering`.
- `File vcf_idx`: Index for vcf.
- `File training_intervals`: Intervals used to train the genotyping model.
- `File median_coverage`: Per-sample median coverage.
- `File rd_file`: Read-depth evidence matrix.
- `File rd_file_idx`: Index for rd_file.
- `File ref_dict`: Reference sequence dictionary.
- `File ploidy_table`: Ploidy table from `DepthPreprocessing`.
- `Array[String] contigs`: Contigs to genotype over, given in reference dictionary order. Model training is also restricted to these, so `training_intervals` and `rd_file` may be genome-wide.
- `String chr_x`: Allosome contig names (defaults `chrX` and `chrY`). (default `chrX`)
- `String chr_y`: Allosome contig names (defaults `chrX` and `chrY`). (default `chrY`)
- `String prefix`: Prefix for output file names.
- `String gatk_docker`, `String sv_base_mini_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (3).

Outputs:
- `File genotyped_depth_vcf`: Genotyped depth CNV VCF.
- `File genotyped_depth_vcf_idx`: Index for `genotyped_depth_vcf`.
- `File genotyping_rd_table`: Read-depth table produced while training the model.

### Callset matching and sharding
`ExactMatch`, `TruvariMatch` and `BedtoolsClosestSV` are the three comparison rounds driven by `AnnotateCallsetOverlap`; each consumes what the previous round left unmatched. `ScatterVcf` is a general sharding helper.

### [ExactMatch](../wdl/utils/ExactMatch.wdl)
This sub-workflow performs the first callset-comparison round, matching records to a truth callset on exact position and allele. Both callsets are optionally renamed to a common ID scheme, sharded, matched, and the annotations concatenated. Records left unmatched are emitted in the form `TruvariMatch` expects.

Inputs:
- `File vcf`: Callset being compared.
- `File vcf_idx`: Index for vcf.
- `File truth_snv_indel_vcf`: Truth callset.
- `File truth_snv_indel_vcf_idx`: Index for truth_snv_indel_vcf.
- `String contig`: Contig being processed.
- `Int? shard_bin_size_exact_match`: Shard size for the matching step.
- `Int min_sv_length_truvari_vcf`: Minimum lengths applied when emitting the Truvari inputs.
- `Int min_sv_length_truvari_truth_vcf`: Minimum lengths applied when emitting the Truvari inputs.
- `String length_field_vcf`: INFO field holding allele length.
- `String source_tag_truth_snv_indel_vcf`: Tag identifying the truth callset in the annotations.
- `String? rename_id_string_vcf`: ID rename templates.
- `String? rename_id_string_truth_snv_indel_vcf`: ID rename templates.
- `Boolean? rename_id_strip_chr_vcf`: Strip the `chr` prefix while renaming.
- `Boolean? rename_id_strip_chr_truth_snv_indel_vcf`: Strip the `chr` prefix while renaming.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (11).

Outputs:
- `File annotated_tsv`: Exact-match annotations.
- `File truvari_eval_vcf`: Unmatched callset records, passed to `TruvariMatch`.
- `File truvari_eval_vcf_idx`: Index for truvari_eval_vcf.
- `File truvari_truth_vcf`: Unmatched truth records, passed to `TruvariMatch`.
- `File truvari_truth_vcf_idx`: Index for truvari_truth_vcf.

### [TruvariMatch](../wdl/utils/TruvariMatch.wdl)
This sub-workflow performs the second comparison round, matching records left unmatched by `ExactMatch` with Truvari at three decreasing sequence-similarity thresholds (0.9, 0.7, 0.5). Each threshold only sees what the previous one failed to match, so a record is annotated with the strictest threshold that matched it.

Inputs:
- `File vcf`: Unmatched callset records from `ExactMatch`.
- `File vcf_idx`: Index for vcf.
- `File truth_snv_indel_vcf`: Unmatched truth records.
- `File truth_snv_indel_vcf_idx`: Index for truth_snv_indel_vcf.
- `String contig`: Contig being processed.
- `String source_tag`: Tag identifying the truth callset in the annotations. (default `SNV_indel`)
- `Int? shard_bin_size_truvari_match`: Shard size for the matching step.
- `Int min_shard_gap_truvari_match`: Minimum gap between records at which a shard boundary may fall. (default `10000`)
- `File? ref_fa`: Reference FASTA and index, when Truvari is run with reference context.
- `File? ref_fai`: Index for ref_fa.
- `String prefix`: Prefix for output file names.
- `String utils_docker`: Container image.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (9).

Outputs:
- `File annotation_tsv`: Truvari match annotations across all three thresholds.
- `File matched_truth_vcf`: Truth records matched by any threshold.
- `File matched_truth_vcf_idx`: Index for matched_truth_vcf.
- `File unmatched_vcf`: Records still unmatched after 0.5, passed to `BedtoolsClosestSV`.
- `File unmatched_vcf_idx`: Index for unmatched_vcf.

### [BedtoolsClosestSV](../wdl/utils/BedtoolsClosestSV.wdl)
This sub-workflow performs the final comparison round, pairing each still-unmatched record with its nearest truth-callset neighbour using `bedtools closest`. Insertions and CNVs are compared separately, since proximity means different things for each, and the two comparisons are merged into a single annotation table.

Inputs:
- `File vcf`: Records left unmatched by `TruvariMatch`.
- `File vcf_idx`: Index for vcf.
- `File truth_sv_vcf`: Truth SV callset.
- `File truth_sv_vcf_idx`: Index for truth_sv_vcf.
- `Int min_sv_length`: Minimum SV length applied to each callset.
- `Int min_sv_length_truth`: Minimum SV length applied to each callset.
- `String type_field`: INFO field holding variant type.
- `String length_field`: INFO field holding allele length.
- `String source_tag`: Tag identifying the truth callset in the annotations. (default `SV`)
- `String prefix`: Prefix for output file names.
- `String gatk_sv_lr_docker`, `String utils_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (8).

Outputs:
- `File annotation_tsv`: Nearest-neighbour annotations for the remaining records.

### [ScatterVcf](../wdl/utils/ScatterVcf.wdl)
This sub-workflow shards a VCF, either by contig, into a fixed number of record-count shards, or both. It can operate on a localized file or stream a remote one, and is used by `AnnotateVEPHail` to parallelize VEP annotation.

Inputs:
- `File file`: VCF or Hail MatrixTable to shard.
- `Int n_shards`: Target shard count. (default `0`)
- `Int records_per_shard`: Target records per shard. (default `0`)
- `String split_vcf_hail_script`: URL of the Hail sharding script; defaults to this repository's copy on `main`. (default `https://raw.githubusercontent.com/talkowski-lab/lr-pipeline/main/scripts/helper/split_vcf_hail.py`)
- `String genome_build`: Reference genome build. (default `GRCh38`)
- `Boolean localize_vcf`: Localize the input rather than streaming it remotely.
- `Boolean get_chromosome_sizes`: Query contig sizes to size the shards.
- `Boolean split_by_chromosome`: Split by contig.
- `Boolean split_into_shards`: Split into record-count shards.
- `Boolean has_index`: Whether the remote input already has an index.
- `String prefix`: Prefix for output file names.
- `String hail_docker`, `String sv_base_mini_docker`: Container images.
- `RuntimeAttr? runtime_attr_*`: Optional per-task runtime overrides (2).

Outputs:
- `Array[File] vcf_shards`: The resulting shards, or the original file when no splitting was requested.
