# Steps
This document is the ordered record of the steps actually run to produce the two released long-read callsets - the combined HPRC/HGSVC set of 292 samples and Phase 1 of All of Us (1027 samples). It is a companion to [Pipeline](pipeline.md), which describes how to run the pipeline on a new cohort; where the two disagree, `pipeline.md` reflects current intent and this document reflects what was run.

Steps are listed in the order they were run. Inputs and outputs are Terra data table columns, in `code`; Terra method configuration names are given in brackets where they differ from the workflow name. Dates are the final successful run of each step. Several stages were executed more than once as upstream inputs were regenerated - only the run that survives in the release is listed, noted where relevant. Cohort-specific steps are tagged; [Cohort Divergences](#cohort-divergences) summarizes them.


## 1. Cohort Callsets
The SNV/indel and SV callsets were built in other workspaces and entered this pipeline as finished cohort VCFs, which is why steps 1 to 3 have no dates recorded here.
1. **[DeepVariant](../wdl/tools/DeepVariant.wdl)** - short variants per sample.
2. **[GLNexus](../wdl/tools/GLNexus.wdl)** - joint genotyping into `snv_indel_vcf`.
3. **SV calling** - PBSV, Sniffles and [PAV](../wdl/tools/PAV.wdl) per sample, integrated within then across samples with Truvari into `sv_vcf`. Filtered by XGBoost _(HPRC/HGSVC)_ or [Kanpig](../wdl/tools/Kanpig.wdl) _(All of Us)_.
4. **[PALMERDiploid](../wdl/tools/PALMERDiploid.wdl)** and **[PALMERAssembly](../wdl/tools/PALMERAssembly.wdl)** - per-sample MEI calls, the latter from assembly BAMs produced by **[MinimapAlignment](../wdl/tools/MinimapAlignment.wdl)**. 2026-02-04.
5. **[MergePALMERCallsets](../wdl/tools/MergePALMERCallsets.wdl)** - `palmer_merged_vcf`. 2026-02-12.
6. **[MosDepth](../wdl/tools/MosDepth.wdl)** - per-base coverage per sample, per contig. 2026-02-12.
7. **[ConcatenateMosDepth](../wdl/annotation_utils/ConcatenateMosDepth.wdl)** - `mosdepth_per_base_combined`. 2026-02-25.
8. **[CreateCohortMetadata](../wdl/annotation_utils/CreateCohortMetadata.wdl)** [`CreateMetadataFile`] - PED and ancestry combined into the cohort metadata TSV. 2026-03-10.
9. **[TRGT](../wdl/tools/TRGT.wdl)** [`TRGT`, `TRGT_Vamos`] - TR genotyping per sample against both catalogs, giving `trgt_trexplorer_vcf` and `trgt_vamos_vcf`. 2026-02-17.
10. **[CombineTRs](../wdl/annotation_utils/CombineTRs.wdl)** - deduplicated overlapping loci into `trgt_combined_vcf`. 2026-03-13.
11. **[AnnotateTREndTags](../wdl/annotation_utils/AnnotateTREndTags.wdl)** [`AddEndTRs`] - added `INFO/END`, giving `trgt_vcf`. 2026-03-14.
12. **[PreprocessVcfs](../wdl/annotation_utils/PreprocessVcfs.wdl)** [`IntegrateVcfs`] - normalized multiallelics, flagged short variant calls at or above 50bp and merged `snv_indel_vcf` with `sv_vcf` into `integrated_vcf`. 2026-04-14.


## 2. Phasing
The whole chain below was run twice; the dates are from the second pass, which reran HiPhase to exclude TRGT homopolymers.
1. **[ExtractSampleVcfs](../wdl/annotation_utils/ExtractSampleVcfs.wdl)** - `integrated_vcf` into `subset_snv_indel_vcf` and `subset_sv_vcf`. 2026-03-15.
2. **[HiPhase](../wdl/tools/HiPhase.wdl)** [`PhysicalPhasing`, `HiPhase_TRGT`] - phased each sample's short variants, SVs and TR calls against its reads. Run with and without the TR VCF, giving `hiphase_vcf` and `hiphase_notrgt_vcf`. 2026-04-22.
3. **[MergeHiPhaseCallsets](../wdl/tools/MergeHiPhaseCallsets.wdl)** [`HiPhaseMerge`] - `bcftools merge` for short variants and SVs, `trgt merge` for TR calls, giving `hiphase_merged_integrated_vcf` and `hiphase_merged_trgt_vcf`. 2026-04-23.
4. **[FillPhasedGenotypes](../wdl/annotation_utils/FillPhasedGenotypes.wdl)** - repopulated `0/0` genotypes from `integrated_vcf`, giving `hiphase_phased_integrated_vcf`. 2026-04-25.
5. **[IntegrateTRs](../wdl/annotation_utils/IntegrateTRs.wdl)** [`AnnotateTRs`] - folded the TR calls back in, giving `tr_annotated_vcf`. 2026-04-26.
6. **[SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)** - sharded per contig into `full_vcf`; all downstream steps run per contig. 2026-04-26.
7. **[BackbonePhase](../wdl/tools/BackbonePhase.wdl)** _(HPRC/HGSVC)_ - transferred phase from `truth_hgsvc_vcf` and `truth_hprc_vcf`, giving `backbone_phased_vcf` and `backbone_phased_notrgt_vcf`. 2026-06-10.
8. **[FillBackbonePhasedGenotypes](../wdl/annotation_utils/FillBackbonePhasedGenotypes.wdl)** [`MergeBackbonePhased`] _(HPRC/HGSVC)_ - filled still-unphased genotypes from the no-TRGT shard, giving `backbone_merged_vcf`. 2026-06-12.
9. **[MethylationProfiling](../wdl/tools/MethylationProfiling.wdl)** - 5mC profiling with pb-CpG-tools, from reads haplotagged by the now-retired Whatshap workflow. 2026-07-12.
10. **[CreateCohortMethylationFile](../wdl/annotation_utils/CreateCohortMethylationFile.wdl)** - per-contig cohort methylation matrices. 2026-07-31.


## 3. Annotation
Each characterization workflow writes an `annotations_tsv_*` column that [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl) later folds into the VCF as INFO fields.
1. **[RepeatMasker](../wdl/tools/RepeatMasker.wdl)** - repeat content of insertion sequences, giving `rm_out` and `rm_fa`. 2026-04-10.
2. **Variant characterization**, run in parallel: [AnnotateL1MEAID](../wdl/annotation/AnnotateL1MEAID.wdl), [AnnotatePALMER](../wdl/annotation/AnnotatePALMER.wdl), [AnnotateSVAN](../wdl/annotation/AnnotateSVAN.wdl), [AnnotateMEDs](../wdl/annotation/AnnotateMEDs.wdl), [AnnotateIndelTRs](../wdl/annotation/AnnotateIndelTRs.wdl), [AnnotateRegion](../wdl/annotation/AnnotateRegion.wdl), [AnnotateInSilicoPredictors](../wdl/annotation/AnnotateInSilicoPredictors.wdl), [AnnotateGnomADSTR](../wdl/annotation/AnnotateGnomADSTR.wdl), [AnnotateDbSNP](../wdl/annotation/AnnotateDbSNP.wdl) [`AnnotateDbGaP`], [AnnotateVRS](../wdl/annotation/AnnotateVRS.wdl), [AnnotateTruvariRemap](../wdl/annotation/AnnotateTruvariRemap.wdl) [`TruvariRemap`], and [AnnotateAgeMetrics](../wdl/annotation/AnnotateAgeMetrics.wdl) _(All of Us)_.
3. **[AnnotateMEIs](../wdl/annotation/AnnotateMEIs.wdl)** - reconciled the L1MEAID, PALMER and SVAN evidence into `annotations_tsv_meis`. 2026-04-12.
4. **[SubsetTsvToColumns](../wdl/annotation_utils/SubsetTsvToColumns.wdl)** - duplications of interest into `annotations_tsv_svan_dups`. 2026-04-12.
5. **[AnnotateAlleleType](../wdl/annotation_utils/AnnotateAlleleType.wdl)** - set `allele_type`, giving `allele_type_annotated_vcf`. 2026-06-21.
6. **[FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl)** - refilled FORMAT fields from `full_vcf`, giving `allele_type_annotated_filled_vcf`. 2026-06-12.
7. **NormalizeDuplicationOrigins** (retired) - rewrote duplication origin coordinates, which required the functional and overlap annotations below to be regenerated.
8. **Functional and overlap annotation**, run in parallel: [AnnotateVEPHail](../wdl/annotation/AnnotateVEPHail.wdl), [AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl), [AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl) and [AnnotateDbVaR](../wdl/annotation/AnnotateDbVaR.wdl). 2026-06-26.
9. **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** [`AnnotateVcf_Functional`, then `AnnotateVcf_Downstream`] - wrote the annotation TSVs into the VCF, giving `annotated_vcf`. 2026-06-27.
10. **[PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)** [`PostProcess`] - transferred phased genotypes from `backbone_merged_vcf`, normalized ploidy, pruned MEIs and flagged homopolymer TRs, giving `post_processed_vcf`. 2026-06-27.
11. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** - allele frequencies from `workspace.ped`, `workspace.par_bed` and `trgt_lps_tsv`, giving `af_annotated_vcf`. 2026-06-27.
12. **FindUntrimmedAlleles** (retired) then **[AnnotateVcfCleared](../wdl/annotation_utils/AnnotateVcfCleared.wdl)** - re-annotated variants whose alleles had been trimmed, giving `untrimmed_vcf` then `cleared_annotated_vcf`. 2026-06-29.
13. **[ResolveHaplotypeOverlaps](../wdl/annotation_utils/ResolveHaplotypeOverlaps.wdl)** _(HPRC/HGSVC)_ - resolved calls overlapping on the same haplotype, giving `overlap_resolved_vcf`. 2026-06-30.
14. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** [`AnnotateAF_PostHoc`] - recomputed frequencies after the genotype changes above, giving `af_annotated_posthoc_vcf`. 2026-07-01.
15. **[NormalizeAlleleTypes](../wdl/annotation_utils/NormalizeAlleleTypes.wdl)** [`TransformAlleleType`] - demoted non-tandem duplications, MEIs and NUMTs while recording `allele_subtype`, giving `transformed_vcf`. 2026-08-21.


## 4. Release
Each step below produces the next numbered release column. The All of Us chain ran the same way but produced `aou_vcf_V1` to `V5`, skipping the steps tagged _(HPRC/HGSVC)_; its dates are not recorded.
1. **[PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)** [`PostProcess_DropFilters`] - dropped the working FILTER values. `transformed_vcf` to `hprc_hgsvc_vcf_V1`. 2026-08-21.
2. **[PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)** [`PostProcess_FilterAssemblySingletons`] _(HPRC/HGSVC)_ - filtered assembly-only singletons. `hprc_hgsvc_vcf_V2`. 2026-08-27.
3. **[CreateCohortDepthFiles](../wdl/annotation_utils/CreateCohortDepthFiles.wdl)** [`CreateDepthFiles`] - bincov matrix, median coverage and ploidy estimates. 2026-08-27.
4. **[IdentifyLowCoverageRegions](../wdl/annotation_utils/IdentifyLowCoverageRegions.wdl)** - 100bp bins at a 90% sample cutoff, giving `failed_bins_bed_90` and `sample_cutoffs_tsv`. 2026-08-28.
5. **[FilterLowCoverageRegions](../wdl/annotation_utils/FilterLowCoverageRegions.wdl)** - filtered variants in the recurrently low-coverage bins. `hprc_hgsvc_vcf_V3`. 2026-08-31.
6. **[AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl)** then **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** [`AnnotateVcfPostHoc`] - refreshed the gnomAD overlap annotations. `hprc_hgsvc_vcf_V4`. 2026-09-10.
7. **[PostProcessTRLociHPRCHGSVC](../wdl/annotation_utils/PostProcessTRLociHPRCHGSVC.wdl)** [`PostprocessTRLoci`] - recovered disease-associated TR loci against the TRExplorer catalog. `hprc_hgsvc_vcf_V5`. 2026-09-15. All of Us used [PostProcessTRLociAoU](../wdl/annotation_utils/PostProcessTRLociAoU.wdl) instead.
8. **[FilterDuplicateZeroDepthReferenceBlocks](../wdl/annotation_utils/FilterDuplicateZeroDepthReferenceBlocks.wdl)**, then **[GLNexus](../wdl/tools/GLNexus.wdl)** and **[SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)** - regenerated the SNV/indel joint callset into `glnexus_vcf`. 2026-09-16.
9. **[FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl)** _(HPRC/HGSVC)_ - refilled `AD`, `PL`, `GQ`, `DP` and `RNC` for DeepVariant records from the regenerated `glnexus_vcf`. `hprc_hgsvc_vcf_V6`. 2026-09-17.
10. **[FilterLowCoverageGenotypes](../wdl/annotation_utils/FilterLowCoverageGenotypes.wdl)** _(HPRC/HGSVC)_ - no-called genotypes below each sample's cutoff in `sample_cutoffs_tsv`. `hprc_hgsvc_vcf_V7`. 2026-09-17.
11. **[AnnotateSQMetrics](../wdl/annotation/AnnotateSQMetrics.wdl)** and **[AnnotateGQMetrics](../wdl/annotation/AnnotateGQMetrics.wdl)**, then **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** [`AnnotateVcfSQGQ`] - recomputed site and genotype quality metrics after the genotype changes above. `hprc_hgsvc_vcf_V8`. 2026-09-17.
12. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** [`AnnotateAFPostHoc`] - recomputed allele frequencies. `hprc_hgsvc_vcf_V9`. 2026-09-18.
13. **[TRGTLPS](../wdl/tools/TRGTLPS.wdl)** then **[CreateTRGTHistograms](../wdl/annotation_utils/CreateTRGTHistograms.wdl)** - `trgt_lps_tsv` and the `trgt_histograms_tsv` browser histograms. 2026-09-18.
14. **[AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl)** then **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** [`AnnotateVcfSVAnnotate`] - refreshed SV functional consequences. `hprc_hgsvc_vcf_V10`, the final HPRC/HGSVC release. 2026-09-19.
15. **[StripGenotypes](../wdl/annotation_utils/StripGenotypes.wdl)** [`DropGenotypes`] _(All of Us)_ - dropped genotypes to produce `aou_sites_vcf`.


## Cohort Divergences
| Step | HPRC/HGSVC | All of Us |
| --- | --- | --- |
| SV callset filtering | XGBoost site filtering | [Kanpig](../wdl/tools/Kanpig.wdl) regenotyping |
| Backbone phasing | Applied | Not applied - no haplotype-resolved base VCFs |
| `AnnotateAgeMetrics` | Not run - no age data | Run |
| TR loci post-processing | [PostProcessTRLociHPRCHGSVC](../wdl/annotation_utils/PostProcessTRLociHPRCHGSVC.wdl), with sequence-agreement phasing | [PostProcessTRLociAoU](../wdl/annotation_utils/PostProcessTRLociAoU.wdl), genotypes unphased |
| `ResolveHaplotypeOverlaps` | Run | Not run |
| `FillFormatFields` | Run | Not run |
| `FilterLowCoverageGenotypes` | Run | Not run |
| Assembly-only singleton filter | Run | Not run |
| Sites-only VCF | Not produced | `aou_sites_vcf` |
| Release versions | `hprc_hgsvc_vcf_V1` to `V10` | `aou_vcf_V1` to `V5`, plus `aou_sites_vcf` |


## Evidence and Gaps
Reconstructed from the Terra job history export at `data/archive/migration/LR_GNOMAD-AoU_TALK_Annotation-Pipeline/job_history.tsv` (23800 rows, 2025-12-09 to 2026-07-31; local-only and gitignored), the submission history and method configurations of the Terra workspace `LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_Annotation-Pipeline` (249 submissions, 2026-08-03 to 2026-09-21), and the pipeline diagrams.

- **All of Us lineage is evidence-light.** Both sources cover the HPRC/HGSVC callset; despite its name, that workspace holds 292 samples and only `hprc_hgsvc_vcf_*` columns. The All of Us column names and divergences come from processing notes rather than a verified data table, and carry no dates.
- **SNV/indel and SV calling are not covered.** They ran in workspaces not inspected here, so steps 1 to 3 of [Cohort Callsets](#1-cohort-callsets) have no dates or configurations.
- **Method configurations were edited in place**, with versions as high as 26, so a configuration's current input mapping does not necessarily match what a given historical run consumed.
- **History starts 2025-12-09.** Anything earlier is not recorded.
- **Retired workflows.** `FindUntrimmedAlleles`, `NormalizeDuplicationOrigins`, `UpdateGenotypes` and `Whatshap` were used for these callsets but have since been retired, so they are named unlinked above.
- **QC and analysis workflows are excluded**, along with unrelated work sharing the same workspace (Paraphase, LPA and SMN1 assembly, Himito, Immuannot) and the depth-based CNV workflows.
