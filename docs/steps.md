# Steps
This document is the historical record of the steps actually executed to produce the two released long-read callsets - the combined HPRC/HGSVC set of 292 samples and Phase 1 of All of Us (1027 samples). It is a companion to [Pipeline](pipeline.md), which describes how to run the pipeline on a new cohort; where the two disagree, `pipeline.md` reflects current intent and this document reflects what was run.

Neither callset was produced by one clean linear pass. Several stages were rerun after upstream fixes, and a number of one-off corrections were applied post hoc that are now handled upstream. Those loops are recorded in [Reruns and Refreshes](#reruns-and-refreshes) rather than being smoothed out of the stage sections below. Inputs and outputs are given as Terra data table columns in `code`, and Terra method configuration names are also in `code` where they differ from the workflow name. Dates are the final successful all-contig or all-sample run.


## Release Versions
Each released VCF is a numbered column on the `LR_contig` data table, produced per contig. The HPRC/HGSVC chain is fully recoverable; the All of Us chain is reconstructed from processing notes only and carries no dates.
| Version | Workflow (Terra config) | Input | Date |
| --- | --- | --- | --- |
| `hprc_hgsvc_vcf_V1` | [PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl) (`PostProcess_DropFilters`) | `transformed_vcf` | 2026-08-21 |
| `hprc_hgsvc_vcf_V2` | [PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl) (`PostProcess_FilterAssemblySingletons`) | `hprc_hgsvc_vcf_V1` | 2026-08-27 |
| `hprc_hgsvc_vcf_V3` | [FilterLowCoverageRegions](../wdl/annotation_utils/FilterLowCoverageRegions.wdl) | `hprc_hgsvc_vcf_V2` | 2026-08-31 |
| `hprc_hgsvc_vcf_V4` | [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl) (`AnnotateVcfPostHoc`) | `hprc_hgsvc_vcf_V3` | 2026-09-10 |
| `hprc_hgsvc_vcf_V5` | [PostProcessTRLociHPRCHGSVC](../wdl/annotation_utils/PostProcessTRLociHPRCHGSVC.wdl) (`PostprocessTRLoci`) | `hprc_hgsvc_vcf_V4` | 2026-09-15 |
| `hprc_hgsvc_vcf_V6` | [FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl) | `hprc_hgsvc_vcf_V5`, `glnexus_vcf` | 2026-09-17 |
| `hprc_hgsvc_vcf_V7` | [FilterLowCoverageGenotypes](../wdl/annotation_utils/FilterLowCoverageGenotypes.wdl) | `hprc_hgsvc_vcf_V6` | 2026-09-17 |
| `hprc_hgsvc_vcf_V8` | [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl) (`AnnotateVcfSQGQ`) | `hprc_hgsvc_vcf_V7` | 2026-09-17 |
| `hprc_hgsvc_vcf_V9` | [AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl) (`AnnotateAFPostHoc`) | `hprc_hgsvc_vcf_V8`, `trgt_lps_tsv` | 2026-09-18 |
| `hprc_hgsvc_vcf_V10` | [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl) (`AnnotateVcfSVAnnotate`) | `hprc_hgsvc_vcf_V9` | 2026-09-19 |
| `aou_vcf_V1` | [PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl) (drop filters) | `transformed_vcf` | not recorded |
| `aou_vcf_V2` | [FilterLowCoverageRegions](../wdl/annotation_utils/FilterLowCoverageRegions.wdl) | `aou_vcf_V1` | not recorded |
| `aou_vcf_V3` | [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl) (callset overlap) | `aou_vcf_V2` | not recorded |
| `aou_vcf_V4` | [PostProcessTRLociAoU](../wdl/annotation_utils/PostProcessTRLociAoU.wdl) | `aou_vcf_V3` | not recorded |
| `aou_vcf_V5` | [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl) (SVAnnotate) | `aou_vcf_V4` | not recorded |
| `aou_sites_vcf` | [StripGenotypes](../wdl/annotation_utils/StripGenotypes.wdl) (`DropGenotypes`) | `aou_vcf_V5` | not recorded |


## 1. Cohort Callset Generation
The SNV/indel and SV callsets were produced outside the annotation workspace and entered this pipeline as finished cohort VCFs, which is why the callset diagram starts from `snv_indel_vcf` and `sv_vcf`. Only the TR callset, the depth files and the mobile element calls were generated here.


### Inputs From Upstream Workspaces
- **Cohort SNV/indel VCF** (`snv_indel_vcf`) - [DeepVariant](../wdl/tools/DeepVariant.wdl) per sample then [GLNexus](../wdl/tools/GLNexus.wdl) joint genotyping. No DeepVariant runs appear in this workspace's history; GLNexus was rerun here in September (see [September GLNexus Rerun](#september-glnexus-rerun)).
- **Cohort SV VCF** (`sv_vcf`) - PBSV, Sniffles and [PAV](../wdl/tools/PAV.wdl) per sample, integrated within and then across samples with Truvari. Filtering diverged by cohort: ML-based (XGBoost) filtering for HPRC/HGSVC, and [Kanpig](../wdl/tools/Kanpig.wdl) regenotyping for All of Us. Kanpig was also trialled in this workspace 2026-03-19 to 2026-04-10.

### TRGT Callset
1. **[TRGT](../wdl/tools/TRGT.wdl)** (`TRGT`, `TRGT_Vamos`; `LR_sample`) - genotyped TR loci per sample against both catalogs, giving `trgt_trexplorer_vcf` and `trgt_vamos_vcf`. 2026-01-30 and 2026-02-17.
2. **[CombineTRs](../wdl/annotation_utils/CombineTRs.wdl)** (`LR_sample_set`) - deduplicated loci overlapping between catalogs into `trgt_combined_vcf`. 2026-03-13.
3. **[AnnotateTREndTags](../wdl/annotation_utils/AnnotateTREndTags.wdl)** (`AddEndTRs`) - added `INFO/END`, giving the per-sample `trgt_vcf`. 2026-03-14.

### Callset Integration
**[PreprocessVcfs](../wdl/annotation_utils/PreprocessVcfs.wdl)** (`IntegrateVcfs`; `LR_sample_set`) - normalized multiallelics, flagged short variant calls at or above 50bp, updated INFO and merged `snv_indel_vcf` with `sv_vcf` into `integrated_vcf`. 2026-04-14.

### Depth Files
1. **[MosDepth](../wdl/tools/MosDepth.wdl)** (`LR_sample`) - per-base coverage per sample per contig. 2026-02-12 onward.
2. **[ConcatenateMosDepth](../wdl/annotation_utils/ConcatenateMosDepth.wdl)** - concatenated each sample's shards into `mosdepth_per_base_combined`. 2026-02-25.
3. **[CreateCohortDepthFiles](../wdl/annotation_utils/CreateCohortDepthFiles.wdl)** (`CreateDepthFiles`; `LR_sample_set`) - bincov matrix, median coverage and ploidy estimates. 2026-08-27.
4. **[IdentifyLowCoverageRegions](../wdl/annotation_utils/IdentifyLowCoverageRegions.wdl)** (`LR_sample_set`) - 100bp bins, `sample_proportion_cutoff` 0.9, producing `failed_bins_bed_90` and `sample_cutoffs_tsv`. Run 2026-08-28, then rerun 2026-09-02 at 50%, 70% and 90% cutoffs; the 90% result is the one promoted to `workspace.low_coverage_regions` and consumed downstream.

### Cohort Metadata
**[CreateCohortMetadata](../wdl/annotation_utils/CreateCohortMetadata.wdl)** (`CreateMetadataFile`) - combined the PED and ancestry files into the cohort metadata TSV used by the TR histograms. 2026-03-10.

### Mobile Element Calls
1. **[MinimapAlignment](../wdl/tools/MinimapAlignment.wdl)** (`LR_sample`) - aligned the assembly haplotypes to GRCh38. 2025-12-11 to 2026-01-14.
2. **[PALMERDiploid](../wdl/tools/PALMERDiploid.wdl)** and **[PALMERAssembly](../wdl/tools/PALMERAssembly.wdl)** (`PALMER`, `PALMER_PreCalled`) - per-sample MEI calls. 2025-12-10 to 2026-02-04.
3. **[MergePALMERCallsets](../wdl/tools/MergePALMERCallsets.wdl)** (`PALMERMerge`) - merged into `palmer_merged_vcf`. 2026-02-12.


## 2. Phasing
Physical phasing with HiPhase ran on both cohorts. Backbone phasing ran on HPRC/HGSVC only, because All of Us has no haplotype-resolved base VCFs. The entire phasing chain was executed twice; the dates below are from the second pass (see [HiPhase V2 Rerun](#hiphase-v2-rerun)).


### Physical Phasing
1. **[ExtractSampleVcfs](../wdl/annotation_utils/ExtractSampleVcfs.wdl)** - split `integrated_vcf` into `subset_snv_indel_vcf` and `subset_sv_vcf`. 2026-03-15.
2. **[HiPhase](../wdl/tools/HiPhase.wdl)** (`PhysicalPhasing`, `HiPhase`, `HiPhase_TRGT`; `LR_sample`) - phased each sample's short variants, SVs and TR calls against its reads, giving `hiphase_vcf` and, for the no-TRGT run, `hiphase_notrgt_vcf`. 1426 runs, through 2026-04-22.
3. **[MergeHiPhaseCallsets](../wdl/tools/MergeHiPhaseCallsets.wdl)** (`HiPhaseMerge`) - `bcftools merge` for short variants and SVs, `trgt merge` for TR calls, giving `hiphase_merged_integrated_vcf` and `hiphase_merged_trgt_vcf`. 2026-04-23.
4. **[FillPhasedGenotypes](../wdl/annotation_utils/FillPhasedGenotypes.wdl)** - repopulated `0/0` genotypes from `integrated_vcf`, giving `hiphase_phased_integrated_vcf`. 2026-04-25.
5. **[IntegrateTRs](../wdl/annotation_utils/IntegrateTRs.wdl)** (`AnnotateTRs`) - folded `hiphase_merged_trgt_vcf` back in, flagging TR-overlapping variants, giving `tr_annotated_vcf`. 2026-04-26.

### Backbone Phasing
Applied to HPRC/HGSVC only. Base VCFs were the HGSVC2024v1.0 combined VCF (`truth_hgsvc_vcf`) and the HPRC v2.0 pangenome wave VCF (`truth_hprc_vcf`).
1. **[SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)** - sharded `tr_annotated_vcf` per contig into `full_vcf`. 2026-04-26.
2. **[BackbonePhase](../wdl/tools/BackbonePhase.wdl)** (`BackbonePhase`, `BackbonePhase_NoTRGT`; `LR_contig`) - transferred phase onto each shard, giving `backbone_phased_vcf` and `backbone_phased_notrgt_vcf`. 2026-04-21 to 2026-06-10.
3. **[FillBackbonePhasedGenotypes](../wdl/annotation_utils/FillBackbonePhasedGenotypes.wdl)** (`MergeBackbonePhased`) - filled still-unphased genotypes from the no-TRGT shard, giving `backbone_merged_vcf`. 2026-06-10 to 2026-06-12.
4. **[EvaluateBackbonePhasing](../wdl/annotation_utils/EvaluateBackbonePhasing.wdl)** - phasing QC tables (`outside_tr_table`, `tr_enveloped_table`, `trv_table`).

`backbone_merged_vcf`, not `backbone_phased_vcf`, is what `PostprocessCallset` later consumed as its genotype transfer source.

### TR Histograms
1. **[TRGTLPS](../wdl/tools/TRGTLPS.wdl)** - longest polymer sequence per locus per sample from `hiphase_merged_trgt_vcf`, giving `trgt_lps_tsv` plus per-contig TRID metadata. 2026-09-18.
2. **[CreateTRGTHistograms](../wdl/annotation_utils/CreateTRGTHistograms.wdl)** (`GenerateTRGTJson`) - stratified per-locus histograms for the TR browser, `trgt_histograms_tsv`. 2026-09-18, after several failed attempts from 2026-09-17.

### Methylation
1. **Whatshap** (retired from this repository) - haplotagged each sample's reads against its phased VCF. 309 runs, 2026-03-26 to 2026-07-11. New cohorts should use HiPhase haplotagging instead.
2. **[MethylationProfiling](../wdl/tools/MethylationProfiling.wdl)** - 5mC profiling with pb-CpG-tools, giving combined and per-haplotype CpG BEDs. 366 runs, through 2026-07-12.
3. **[CreateCohortMethylationFile](../wdl/annotation_utils/CreateCohortMethylationFile.wdl)** - per-contig cohort methylation matrices. 2026-07-30 to 2026-07-31.


## 3. Annotation
Annotation ran per contig on the `LR_contig` table. Each characterization workflow writes an `annotations_tsv_*` column, which `AnnotateVcf` later folds into the VCF as INFO fields.


### Variant Characterization
1. **[RepeatMasker](../wdl/tools/RepeatMasker.wdl)** - repeat content of insertion sequences, giving `rm_out` and `rm_fa`. 2026-04-10.
2. Run in parallel on the phased cohort VCF:
   - **[AnnotateL1MEAID](../wdl/annotation/AnnotateL1MEAID.wdl)** - `annotations_tsv_l1meaid`. 2026-04-10.
   - **[AnnotatePALMER](../wdl/annotation/AnnotatePALMER.wdl)** - matched `palmer_merged_vcf` against callset insertions using `rm_out`, giving `annotations_tsv_palmer`. 2026-04-10.
   - **[AnnotateSVAN](../wdl/annotation/AnnotateSVAN.wdl)** - `annotations_tsv_svan` and `annotations_header_svan`. 2026-06-24.
   - **[AnnotateMEDs](../wdl/annotation/AnnotateMEDs.wdl)** - `annotations_tsv_meds`. 2026-04-10.
   - **[AnnotateIndelTRs](../wdl/annotation/AnnotateIndelTRs.wdl)** - `annotations_tsv_trs`. 2026-07-21.
   - **[AnnotateRegion](../wdl/annotation/AnnotateRegion.wdl)** - `annotations_tsv_region`. 2026-05-12.
   - **[AnnotateInSilicoPredictors](../wdl/annotation/AnnotateInSilicoPredictors.wdl)** - `annotations_tsv_insilico`. 2026-04-10.
   - **[AnnotateGnomADSTR](../wdl/annotation/AnnotateGnomADSTR.wdl)** - `annotations_tsv_gnomad_str`. 2026-04-10.
   - **[AnnotateDbSNP](../wdl/annotation/AnnotateDbSNP.wdl)** (`AnnotateDbGaP`) - `annotations_tsv_dbsnp`. 2026-04-10.
   - **[AnnotateVRS](../wdl/annotation/AnnotateVRS.wdl)** - `annotations_tsv_vrs`. 2026-06-11.
   - **[AnnotateTruvariRemap](../wdl/annotation/AnnotateTruvariRemap.wdl)** (`TruvariRemap`) - `annotations_tsv_remap`. 2026-07-17.
   - **[AnnotateAgeMetrics](../wdl/annotation/AnnotateAgeMetrics.wdl)** _(All of Us)_ - `annotations_tsv_age`. Not run for HPRC/HGSVC, which has no age data.

### MEI Consolidation and Allele Typing
1. **[AnnotateMEIs](../wdl/annotation/AnnotateMEIs.wdl)** - reconciled the L1MEAID, PALMER and SVAN evidence into `annotations_tsv_meis`. 2026-04-12.
2. **[SubsetTsvToColumns](../wdl/annotation_utils/SubsetTsvToColumns.wdl)** - pulled the duplications of interest out of the SVAN output into `annotations_tsv_svan_dups`. 2026-04-12.
3. **[AnnotateAlleleType](../wdl/annotation_utils/AnnotateAlleleType.wdl)** - set `allele_type` from the MED, MEI and duplication TSVs, giving `allele_type_annotated_vcf`. 2026-06-21.
4. **[FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl)** - refilled FORMAT fields from `full_vcf`, giving `allele_type_annotated_filled_vcf`. 2026-06-12.
5. **NormalizeDuplicationOrigins** (retired from this repository) - rewrote duplication origin coordinates. This invalidated the existing functional and overlap annotations and forced the refresh described in [Annotation Refresh](#annotation-refresh).

### Functional and Overlap Annotation
Run in parallel on the allele-type-annotated VCF:
- **[AnnotateVEPHail](../wdl/annotation/AnnotateVEPHail.wdl)** - `annotations_tsv_vep`. 2026-04-13.
- **[AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl)** - `annotations_tsv_svannotate`. 2026-09-19.
- **[AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl)** - `annotations_tsv_benchmark` and `annotations_header_benchmark`. 2026-09-02, with OOM retries through 2026-09-06.
- **[AnnotateDbVaR](../wdl/annotation/AnnotateDbVaR.wdl)** - `annotations_tsv_dbvar`. 2026-06-25.

### Integration and Post-Processing
1. **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** (`AnnotateVcf_Functional`, then `AnnotateVcf_Downstream`) - wrote the annotation TSVs into the VCF as INFO fields, giving `annotated_vcf`. 2026-06-14 and 2026-06-27.
2. **[PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)** (`PostProcess`) - transferred phased genotypes from `backbone_merged_vcf`, normalized ploidy, pruned MEIs, flagged homopolymer TRs, giving `post_processed_vcf`. 2026-06-27.
3. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** - allele frequencies from `workspace.ped`, `workspace.par_bed` and `trgt_lps_tsv`, giving `af_annotated_vcf`. 2026-06-27.
4. **FindUntrimmedAlleles** (retired) then **[AnnotateVcfCleared](../wdl/annotation_utils/AnnotateVcfCleared.wdl)** - re-annotated variants whose alleles had been trimmed, giving `untrimmed_vcf` then `cleared_annotated_vcf`. 2026-06-28 and 2026-06-29. See [Untrim Loop](#untrim-loop).
5. **[ResolveHaplotypeOverlaps](../wdl/annotation_utils/ResolveHaplotypeOverlaps.wdl)** _(HPRC/HGSVC)_ - resolved calls overlapping on the same haplotype, giving `overlap_resolved_vcf` and `overlap_tsv`. 2026-06-29 to 2026-06-30.
6. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** (`AnnotateAF_PostHoc`) - recomputed frequencies after the genotype changes above, giving `af_annotated_posthoc_vcf`. 2026-07-01.
7. **[NormalizeAlleleTypes](../wdl/annotation_utils/NormalizeAlleleTypes.wdl)** (`TransformAlleleType`, `TransformDuplications`) - demoted non-tandem duplications, MEIs and NUMTs while recording `allele_subtype`, giving `transformed_vcf`. 2026-08-21.

`transformed_vcf` is the input to `hprc_hgsvc_vcf_V1`; from there the release chain is the [Release Versions](#release-versions) table.


## Reruns and Refreshes
These loops are the reason the stage sections above carry dates spread across several months.


### HiPhase V2 Rerun
The whole phasing and post-processing chain was rerun to exclude TRGT homopolymers from HiPhase. The second pass ran HiPhase through 2026-04-22, `HiPhaseMerge` 2026-04-23, `FillPhasedGenotypes` 2026-04-25, `AnnotateTRs` 2026-04-26, `SplitVcfPerContig` 2026-04-26, `BackbonePhase` through 2026-06-10, `MergeBackbonePhased` 2026-06-10 to 2026-06-12, and `PostProcess` 2026-06-10 to 2026-06-27.

### FillFormatFields Passes
`FillFormatFields` ran four times against different sources, because FORMAT fields kept being lost or superseded as upstream VCFs were regenerated.
1. 2026-04-17 to 2026-04-19 - unfilled `annotated_vcf`, filled from `full_vcf`.
2. 2026-06-09 to 2026-06-12 - unfilled `allele_type_annotated_vcf`, filled from `full_vcf`, giving `allele_type_annotated_filled_vcf`.
3. 2026-09-03 to 2026-09-11 - several attempts, all superseded by the GLNexus rerun.
4. 2026-09-16 to 2026-09-17 - unfilled `hprc_hgsvc_vcf_V5`, filled from `glnexus_vcf`, restricted to `INFO/SOURCE` equal to `DeepVariant`, giving `hprc_hgsvc_vcf_V6`. This is the pass that survives in the release.

### Untrim Loop
The retired `FindUntrimmedAlleles` plus [AnnotateVcfCleared](../wdl/annotation_utils/AnnotateVcfCleared.wdl) pair ran three times - 2026-05-25 on `annotated_vcf`, 2026-06-12 to 2026-06-16 on `post_processed_vcf`, and finally 2026-06-28 to 2026-06-29. Only the last pass feeds the release. Each pass required re-running `AnnotateCallsetOverlap`, `AnnotateDbSNP`, `AnnotateDbVaR`, `AnnotateInSilicoPredictors` and `AnnotateVRS` against the untrimmed variants, recorded in the `untrimmed_annotations_tsv_*` columns.

### Annotation Refresh
After the retired `NormalizeDuplicationOrigins` step rewrote duplication coordinates, the affected annotations were regenerated: `AnnotateSVAnnotate` 2026-05-29 to 2026-06-26, `AnnotateDbVaR` 2026-05-28 to 2026-06-25, `AnnotateVcf_Functional` 2026-06-11 to 2026-06-14, `AnnotateVcf_Downstream` 2026-06-26 to 2026-06-27, and `AnnotateCallsetOverlap` through 2026-07-25.

### September GLNexus Rerun
The SNV/indel joint callset was regenerated late, which invalidated `hprc_hgsvc_vcf_V6` through `V9` and forced them to be rebuilt.
1. **[FilterDuplicateZeroDepthReferenceBlocks](../wdl/annotation_utils/FilterDuplicateZeroDepthReferenceBlocks.wdl)** - cleaned the per-sample gVCFs, 292 samples. 2026-09-11.
2. **[GLNexus](../wdl/tools/GLNexus.wdl)** - rejoint-genotyped all samples. 2026-09-16.
3. **[SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)** - sharded the result into `glnexus_vcf`. 2026-09-16.

Earlier attempts at `FillFormatFields`, `FilterLowCoverageGenotypes`, `AnnotateSQMetrics`, `AnnotateGQMetrics` and `AnnotateAFPostHoc` from 2026-09-03 through 2026-09-14 were all discarded and redone after this rerun.


## Cohort Divergences
| Step | HPRC/HGSVC | All of Us |
| --- | --- | --- |
| SV callset filtering | ML-based (XGBoost) site filtering | [Kanpig](../wdl/tools/Kanpig.wdl) regenotyping, non-ref support required |
| Backbone phasing | Applied, using the HGSVC2024v1.0 and HPRC v2.0 base VCFs | Not applied - no haplotype-resolved base VCFs |
| `AnnotateAgeMetrics` | Not run - no age data | Run |
| TR loci post-processing | [PostProcessTRLociHPRCHGSVC](../wdl/annotation_utils/PostProcessTRLociHPRCHGSVC.wdl), with sequence-agreement phasing | [PostProcessTRLociAoU](../wdl/annotation_utils/PostProcessTRLociAoU.wdl), genotypes emitted unphased |
| `ResolveHaplotypeOverlaps` | Run | Not run |
| `FillFormatFields` | Run | Not run |
| `FilterLowCoverageGenotypes` | Run | Not run |
| Assembly-only singleton filter | Run (`PostProcess_FilterAssemblySingletons`) | Not run |
| Sites-only VCF | Not produced | `aou_sites_vcf` via `DropGenotypes` |
| Release versions | `hprc_hgsvc_vcf_V1` to `V10` | `aou_vcf_V1` to `V5`, plus `aou_sites_vcf` |


## Evidence and Gaps
Reconstructed from four sources: the Terra job history export at `data/archive/migration/LR_GNOMAD-AoU_TALK_Annotation-Pipeline/job_history.tsv` (23800 rows, 2025-12-09 to 2026-07-31; local-only and gitignored), the submission history of the Terra workspace `LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_Annotation-Pipeline` (249 submissions, 2026-08-03 to 2026-09-21), the method configurations in that workspace, and the pipeline diagrams above.

Known gaps and caveats:
- **All of Us lineage is evidence-light.** Both sources cover the HPRC/HGSVC callset; despite its name, the workspace has 292 samples and only `hprc_hgsvc_vcf_*` columns. The All of Us column names, ordering and divergences here come from processing notes rather than a verified data table, and carry no dates.
- **SNV/indel and SV calling are not covered.** DeepVariant and the SV callers ran in workspaces not inspected here, so no dates or configurations are recorded for them.
- **Method configurations were edited in place.** Versions run as high as 26, so a configuration's current input mapping does not necessarily match what a given historical run consumed. Where the two disagree, the version chain in the notes was preferred.
- **History starts 2025-12-09.** Anything earlier is not recorded.
- **Retired workflows.** `FindUntrimmedAlleles`, `NormalizeDuplicationOrigins` and `UpdateGenotypes` were used for these callsets but have since been retired; they are named unlinked above.
- **Terra configuration names drift from workflow names.** `IntegrateVcfs` is `PreprocessVcfs`, `AddEndTRs` is `AnnotateTREndTags`, `AnnotateTRs` is `IntegrateTRs`, `MergeBackbonePhased` is `FillBackbonePhasedGenotypes`, `GenerateTRGTJson` is `CreateTRGTHistograms`, `CreateDepthFiles` is `CreateCohortDepthFiles`, `AnnotateDbGaP` is `AnnotateDbSNP`, and every `PostProcess_*`, `AnnotateVcf_*` and `AnnotateAF*` configuration is the same workflow with different toggles or TSV sets.
- **QC and analysis workflows are excluded.** `CountAnnotations`, `SummarizeAnnotations`, `DiagnoseSingletons`, `ExtractRandomCalls`, `VcfToBedVepParsed`, `ExtractIGVVariants`, `RunVisualizePlots`, `BenchmarkIndividualVcf`, `PlotPhasingResults`, `EvaluateOverlappingTRLoci`, `Automop`, the depth-based CNV workflows and unrelated work in the same workspace (Paraphase, LPA and SMN1 assembly, Himito, Immuannot) are not part of the callset path.
