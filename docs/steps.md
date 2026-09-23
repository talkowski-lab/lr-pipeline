# Steps
This document is the ordered record of the steps actually run to produce the two released long-read callsets - the combined HPRC/HGSVC set of 292 samples and Phase 1 of All of Us (1027 samples). It is a companion to [Pipeline](pipeline.md), which describes how to run the pipeline on a new cohort; where the two disagree, `pipeline.md` reflects current intent and this document reflects what was run.

Steps are listed in the order they were run. Inputs and outputs are Terra data table columns, in `code`. Workflows are named as they appear in this repository and linked to their WDL; the Terra method configuration that ran a step is often named differently and is not recorded here. Dates are the final successful run of each step. Several stages were executed more than once as upstream inputs were regenerated - only the run that survives in the release is listed, noted where relevant.

The two callsets share sections 1 to 3 apart from the steps tagged _(HPRC/HGSVC)_ or _(All of Us)_, which ran for that cohort only. They diverge substantially in section 4, so [Release](#4-release) gives each step's output per cohort side by side. [Cohort Divergences](#cohort-divergences) summarizes every divergence in one place.


## 1. Cohort Callsets
Both callsets entered this pipeline as finished cohort-level VCFs - `snv_indel_vcf` and `sv_vcf` - produced in other workspaces; for HPRC/HGSVC these were a joint-called short variant VCF and `GRCh38_INSDEL_1218.vcf.gz`, both read straight out of the upstream workspace bucket. Everything below prepares the remaining per-sample evidence and integrates it into the single cohort VCF that phasing consumes.
1. **[PALMERDiploid](../wdl/tools/PALMERDiploid.wdl)** and **[PALMERAssembly](../wdl/tools/PALMERAssembly.wdl)** - per-sample MEI calls, the latter from assembly BAMs produced by **[MinimapAlignment](../wdl/tools/MinimapAlignment.wdl)**. 2026-02-04.
2. **[MergePALMERCallsets](../wdl/tools/MergePALMERCallsets.wdl)** - `palmer_merged_vcf`. 2026-02-12.
3. **[MosDepth](../wdl/tools/MosDepth.wdl)** - per-base coverage per sample, per contig. 2026-02-12.
4. **[ConcatenateMosDepth](../wdl/annotation_utils/ConcatenateMosDepth.wdl)** - `mosdepth_per_base_combined`. 2026-02-25.
5. **[CreateCohortMetadata](../wdl/annotation_utils/CreateCohortMetadata.wdl)** - PED and ancestry combined into the cohort metadata TSV. 2026-03-10.
6. **[TRGT](../wdl/tools/TRGT.wdl)** - TR genotyping per sample against both catalogs, giving `trgt_trexplorer_vcf` and `trgt_vamos_vcf`. 2026-02-17.
7. **[CombineTRs](../wdl/annotation_utils/CombineTRs.wdl)** - deduplicated overlapping loci into `trgt_combined_vcf`. 2026-03-13.
8. **[AnnotateTREndTags](../wdl/annotation_utils/AnnotateTREndTags.wdl)** - added `INFO/END`, giving `trgt_vcf`. 2026-03-14.
9. **[Kanpig](../wdl/tools/Kanpig.wdl)** _(HPRC/HGSVC)_ - regenotyped every SV site against each sample's reads, giving `kanpig_vcf`. The supplied `sv_vcf` carried no `AD`, `DP` or `PL` on its reference and empty calls, so these calls had to be regenotyped locally to recover those FORMAT fields, which reached the callset through [AnnotateSvCallerSupport](../wdl/annotation_utils/AnnotateSvCallerSupport.wdl) in [Annotation](#3-annotation) step 6. 2026-04-10.
10. **[PreprocessVcfs](../wdl/annotation_utils/PreprocessVcfs.wdl)** - normalized multiallelics, flagged short variant calls at or above 50bp and merged `snv_indel_vcf` with `sv_vcf` into `integrated_vcf`. 2026-04-14.


## 2. Phasing
The whole chain below was run twice, the second pass rerunning HiPhase to exclude TRGT homopolymers. The dates are the runs whose output survives in the release: tracing the released VCF back shows its sites came through the first pass, and the second pass reaches the release only through the phased genotypes that [PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl) transfers from `backbone_merged_vcf` in [Annotation](#3-annotation) step 11.
1. **[ExtractSampleVcfs](../wdl/annotation_utils/ExtractSampleVcfs.wdl)** - `integrated_vcf` into `subset_snv_indel_vcf` and `subset_sv_vcf`. 2026-03-15.
2. **[HiPhase](../wdl/tools/HiPhase.wdl)** - phased each sample's short variants, SVs and TR calls against its reads. Run with and without the TR VCF, giving `hiphase_vcf` and `hiphase_notrgt_vcf`. 2026-03-16.
3. **[MergeHiPhaseCallsets](../wdl/tools/MergeHiPhaseCallsets.wdl)** - `bcftools merge` for short variants and SVs, `trgt merge` for TR calls, giving `hiphase_merged_integrated_vcf` and `hiphase_merged_trgt_vcf`. 2026-03-18.
4. **[FillPhasedGenotypes](../wdl/annotation_utils/FillPhasedGenotypes.wdl)** - repopulated `0/0` genotypes from `integrated_vcf`, giving `hiphase_phased_integrated_vcf`. 2026-03-19.
5. **[IntegrateTRs](../wdl/annotation_utils/IntegrateTRs.wdl)** - folded the TR calls back in, giving `tr_annotated_vcf`. 2026-04-02.
6. **[SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)** - sharded per contig into `full_vcf`; all downstream steps run per contig. 2026-04-06.
7. **[BackbonePhase](../wdl/tools/BackbonePhase.wdl)** _(HPRC/HGSVC)_ - transferred phase from `truth_hgsvc_vcf` and `truth_hprc_vcf`, giving `backbone_phased_vcf` and `backbone_phased_notrgt_vcf`. 2026-06-10.
8. **[FillBackbonePhasedGenotypes](../wdl/annotation_utils/FillBackbonePhasedGenotypes.wdl)** _(HPRC/HGSVC)_ - filled still-unphased genotypes from the no-TRGT shard, giving `backbone_merged_vcf`. 2026-06-12.
9. **[MethylationProfiling](../wdl/tools/MethylationProfiling.wdl)** - 5mC profiling with pb-CpG-tools, from reads haplotagged by **[Whatshap](../archive/wdl/tools/Whatshap.wdl)**, since archived. 2026-07-12.
10. **[CreateCohortMethylationFile](../wdl/annotation_utils/CreateCohortMethylationFile.wdl)** - per-contig cohort methylation matrices. 2026-07-31.


## 3. Annotation
Each characterization workflow writes an `annotations_tsv_*` column that [AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl) later folds into the VCF as INFO fields.
1. **[RepeatMasker](../wdl/tools/RepeatMasker.wdl)** - repeat content of insertion sequences, giving `rm_out` and `rm_fa`. 2026-04-10.
2. **Variant characterization**, run in parallel: [AnnotateL1MEAID](../wdl/annotation/AnnotateL1MEAID.wdl), [AnnotatePALMER](../wdl/annotation/AnnotatePALMER.wdl), [AnnotateSVAN](../wdl/annotation/AnnotateSVAN.wdl), [AnnotateMEDs](../wdl/annotation/AnnotateMEDs.wdl), [AnnotateIndelTRs](../wdl/annotation/AnnotateIndelTRs.wdl), [AnnotateRegion](../wdl/annotation/AnnotateRegion.wdl), [AnnotateInSilicoPredictors](../wdl/annotation/AnnotateInSilicoPredictors.wdl), [AnnotateGnomADSTR](../wdl/annotation/AnnotateGnomADSTR.wdl), [AnnotateDbSNP](../wdl/annotation/AnnotateDbSNP.wdl), [AnnotateVRS](../wdl/annotation/AnnotateVRS.wdl), [AnnotateTruvariRemap](../wdl/annotation/AnnotateTruvariRemap.wdl), and [AnnotateAgeMetrics](../wdl/annotation/AnnotateAgeMetrics.wdl) _(All of Us)_.
3. **[AnnotateMEIs](../wdl/annotation/AnnotateMEIs.wdl)** - reconciled the L1MEAID, PALMER and SVAN evidence into `annotations_tsv_meis`. 2026-04-12.
4. **[SubsetTsvToColumns](../wdl/annotation_utils/SubsetTsvToColumns.wdl)** - duplications of interest into `annotations_tsv_svan_dups`. 2026-04-12.
5. **[AnnotateAlleleType](../wdl/annotation_utils/AnnotateAlleleType.wdl)** - set `allele_type`, giving `allele_type_annotated_vcf`. 2026-06-21.
6. **[AnnotateSvCallerSupport](../wdl/annotation_utils/AnnotateSvCallerSupport.wdl)** _(HPRC/HGSVC)_ - merged the per-sample `kanpig_vcf` calls from [Cohort Callsets](#1-cohort-callsets) step 9 back onto the cohort `sv_vcf`, then **[SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)** re-sharded the result. This Kanpig-bearing `full_vcf` is the donor the next step copies from, and is how the missing `AD`, `DP` and `PL` reached the callset. 2026-06-10.
7. **[FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl)** - copied the FORMAT fields across from the `full_vcf` built in step 6, giving `allele_type_annotated_filled_vcf`. 2026-06-12.
8. **[NormalizeDuplicationOrigins](../wdl/annotation_utils/NormalizeDuplicationOrigins.wdl)** - rewrote duplication origin coordinates, which required the functional and overlap annotations below to be regenerated.
9. **Functional and overlap annotation**, run in parallel: [AnnotateVEPHail](../wdl/annotation/AnnotateVEPHail.wdl), [AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl), [AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl) and [AnnotateDbVaR](../wdl/annotation/AnnotateDbVaR.wdl). 2026-06-26.
10. **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** - run twice, writing first the functional and then the remaining annotation TSVs into the VCF, giving `annotated_vcf`. 2026-06-27.
11. **[PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)** - transferred phased genotypes from `backbone_merged_vcf`, normalized ploidy, pruned MEIs and flagged homopolymer TRs, giving `post_processed_vcf`. 2026-06-27.
12. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** - allele frequencies from `workspace.ped`, `workspace.par_bed` and `trgt_lps_tsv`, giving `af_annotated_vcf`. 2026-06-27.
13. **[FindUntrimmedAlleles](../wdl/annotation_utils/FindUntrimmedAlleles.wdl)** then **[AnnotateVcfCleared](../wdl/annotation_utils/AnnotateVcfCleared.wdl)** - re-annotated variants whose alleles had been trimmed, giving `untrimmed_vcf` then `cleared_annotated_vcf`. 2026-06-29.
14. **[ResolveHaplotypeOverlaps](../wdl/annotation_utils/ResolveHaplotypeOverlaps.wdl)** _(HPRC/HGSVC)_ - resolved calls overlapping on the same haplotype, giving `overlap_resolved_vcf`. 2026-06-30.
15. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** _(HPRC/HGSVC)_ - recomputed frequencies after the genotype changes above, giving `af_annotated_posthoc_vcf`. 2026-07-01.
16. **[NormalizeAlleleTypes](../wdl/annotation_utils/NormalizeAlleleTypes.wdl)** - demoted non-tandem duplications, MEIs and NUMTs while recording `allele_subtype`, giving `transformed_vcf`. 2026-08-21.


## 4. Release
This is where the two callsets diverge most, and where their version numbering stops lining up - HPRC/HGSVC ran ten release steps, All of Us five. The table gives, for each step, the release column that step wrote for each cohort. `–` means the step was not run for that cohort at all. `run, no column` means the step produced supporting files rather than a release VCF. Dates are the HPRC/HGSVC runs; the All of Us runs are undated.

| # | Step | HPRC/HGSVC | All of Us | Date |
| --- | --- | --- | --- | --- |
| 1 | **[PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)** - dropped the working FILTER values, from `transformed_vcf` | `hprc_hgsvc_vcf_V1` | `aou_vcf_V1` | 2026-08-21 |
| 2 | **[PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl)** - filtered assembly-only singletons | `hprc_hgsvc_vcf_V2` | – | 2026-08-27 |
| 3 | **[CreateCohortDepthFiles](../wdl/annotation_utils/CreateCohortDepthFiles.wdl)** - bincov matrix, median coverage and ploidy estimates | run, no column | run, no column | 2026-08-27 |
| 4 | **[IdentifyLowCoverageRegions](../wdl/annotation_utils/IdentifyLowCoverageRegions.wdl)** - 100bp bins at a 90% sample cutoff, giving `failed_bins_bed_90` and `sample_cutoffs_tsv` | run, no column | run, no column | 2026-08-28 |
| 5 | **[FilterLowCoverageRegions](../wdl/annotation_utils/FilterLowCoverageRegions.wdl)** - filtered variants in the recurrently low-coverage bins | `hprc_hgsvc_vcf_V3` | `aou_vcf_V2` | 2026-08-31 |
| 6 | **[AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl)** then **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** - refreshed the gnomAD overlap annotations | `hprc_hgsvc_vcf_V4` | `aou_vcf_V3` | 2026-09-10 |
| 7 | TR loci post-processing - recovered disease-associated TR loci. HPRC/HGSVC ran **[PostProcessTRLociHPRCHGSVC](../wdl/annotation_utils/PostProcessTRLociHPRCHGSVC.wdl)** against the TRExplorer catalog; All of Us ran **[PostProcessTRLociAoU](../wdl/annotation_utils/PostProcessTRLociAoU.wdl)** | `hprc_hgsvc_vcf_V5` | `aou_vcf_V4` | 2026-09-15 |
| 8 | **[FilterDuplicateZeroDepthReferenceBlocks](../wdl/annotation_utils/FilterDuplicateZeroDepthReferenceBlocks.wdl)**, then **[GLNexus](../wdl/tools/GLNexus.wdl)** and **[SplitVcfPerContig](../wdl/annotation_utils/SplitVcfPerContig.wdl)** - regenerated the SNV/indel joint callset into `glnexus_vcf`, solely to feed step 9 | run, no column | – | 2026-09-15 |
| 9 | **[FillFormatFields](../wdl/annotation_utils/FillFormatFields.wdl)** - copied `AD`, `PL`, `GQ`, `DP` and `RNC` for the short variant records across from the regenerated `glnexus_vcf` | `hprc_hgsvc_vcf_V6` | – | 2026-09-16 |
| 10 | **[FilterLowCoverageGenotypes](../wdl/annotation_utils/FilterLowCoverageGenotypes.wdl)** - no-called genotypes below each sample's cutoff in `sample_cutoffs_tsv` | `hprc_hgsvc_vcf_V7` | – | 2026-09-17 |
| 11 | **[AnnotateSQMetrics](../wdl/annotation/AnnotateSQMetrics.wdl)** and **[AnnotateGQMetrics](../wdl/annotation/AnnotateGQMetrics.wdl)**, then **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** - recomputed site and genotype quality metrics after the genotype changes in steps 9 and 10 | `hprc_hgsvc_vcf_V8` | – | 2026-09-17 |
| 12 | **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** - recomputed allele frequencies | `hprc_hgsvc_vcf_V9` | – | 2026-09-18 |
| 13 | **[TRGTLPS](../wdl/tools/TRGTLPS.wdl)** then **[CreateTRGTHistograms](../wdl/annotation_utils/CreateTRGTHistograms.wdl)** - `trgt_lps_tsv` and the `trgt_histograms_tsv` browser histograms | run, no column | not recorded | 2026-09-18 |
| 14 | **[AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl)** then **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** - refreshed SV functional consequences. Final release for both cohorts | `hprc_hgsvc_vcf_V10` | `aou_vcf_V5` | 2026-09-19 |
| 15 | **[DropGenotypes](../wdl/annotation_utils/DropGenotypes.wdl)** - dropped genotypes to produce a sites-only VCF | – | `aou_sites_vcf` | – |

Steps 2 and 8 to 12 are HPRC/HGSVC only, which is why its chain reaches `V10` while All of Us stops at `V5`. Steps 1, 5, 6, 7 and 14 advanced both.


## Cohort Divergences
| Step | HPRC/HGSVC | All of Us |
| --- | --- | --- |
| SV FORMAT fields | Absent on reference and empty calls in the supplied `sv_vcf`, so [Kanpig](../wdl/tools/Kanpig.wdl) was run locally and its calls merged back by [AnnotateSvCallerSupport](../wdl/annotation_utils/AnnotateSvCallerSupport.wdl) | Already present |
| Backbone phasing | Applied | Not applied - no haplotype-resolved base VCFs |
| `AnnotateAgeMetrics` | Not run - no age data | Run |
| TR loci post-processing | [PostProcessTRLociHPRCHGSVC](../wdl/annotation_utils/PostProcessTRLociHPRCHGSVC.wdl), with sequence-agreement phasing | [PostProcessTRLociAoU](../wdl/annotation_utils/PostProcessTRLociAoU.wdl), genotypes unphased |
| `ResolveHaplotypeOverlaps` | Run | Not run |
| Post-hoc allele frequencies | Recomputed twice, after the overlap resolution and again after the genotype filters | Not recomputed |
| SNV/indel joint callset regeneration | Rerun through `GLNexus` to feed `FillFormatFields` | Not rerun |
| `FillFormatFields` | Run | Not run |
| `FilterLowCoverageGenotypes` | Run | Not run |
| Assembly-only singleton filter | Run | Not run |
| Sites-only VCF | Not produced | `aou_sites_vcf` via [DropGenotypes](../wdl/annotation_utils/DropGenotypes.wdl) |
| Release versions | `hprc_hgsvc_vcf_V1` to `V10` | `aou_vcf_V1` to `V5`, plus `aou_sites_vcf` |


## Evidence and Gaps
Reconstructed from the Terra job history export at `data/archive/migration/LR_GNOMAD-AoU_TALK_Annotation-Pipeline/job_history.tsv` (23800 rows, 2025-12-09 to 2026-07-31; local-only and gitignored), the submission history and method configurations of the Terra workspace `LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_Annotation-Pipeline` (249 submissions, 2026-08-03 to 2026-09-21), and the pipeline diagrams.

- **All of Us lineage is evidence-light.** The Terra sources cover the HPRC/HGSVC callset only; despite its name, that workspace holds 292 samples and only `hprc_hgsvc_vcf_*` columns. Every All of Us column name and per-step divergence above comes from [`scratch.md`](scratch.md) rather than a verified data table, and carries no dates. Sections 1 to 3 are assumed shared except where tagged, since those notes only begin at the release chain.
- **All of Us supporting steps are inferred.** `CreateCohortDepthFiles` and `IdentifyLowCoverageRegions` are marked as run for All of Us because `FilterLowCoverageRegions`, which did run, cannot proceed without them; the notes do not name them. The GLNexus regeneration is marked as not run because its only consumer, `FillFormatFields`, is HPRC/HGSVC only.
- **SNV/indel and SV calling are not covered.** Both callsets were built in workspaces not inspected here and arrive as `snv_indel_vcf` and `sv_vcf`, so nothing upstream of [Cohort Callsets](#1-cohort-callsets) is recorded - including how the SV callset was filtered.
- **The HPRC/HGSVC chain is verified by lineage trace.** Every step above was confirmed by taking a released `hprc_hgsvc_vcf_V10` shard - `chr20.chr20.annotated.vcf.gz` from the final `AnnotateVcf` submission - and walking its input VCF back 27 hops to the upstream cohort VCFs, using the Terra submission API for the August and September runs and the `inputs` column of the job history export for everything earlier. Paths under the current bucket and the older `fc-fd42e80c` bucket both resolve, since the copy preserved the `submissions/<id>/` layout and therefore the submission IDs.
- **Two Kanpig transfer attempts did not survive.** `ReplaceKanpigGT`, still in [`archive/`](../archive/wdl/annotation_utils/ReplaceKanpigGT.wdl), ran three times on chr22 alone on 2026-06-01 and 2026-06-02 and nothing ever consumed its output; a `FillSVFormatFields` configuration consumed Kanpig output once on 2026-04-13 and likewise does not appear in the traced lineage. Neither is listed as a step.
- **Method configurations were edited in place**, with versions as high as 26, so a configuration's current input mapping does not necessarily match what a given historical run consumed.
- **History starts 2025-12-09.** Anything earlier is not recorded.
- **Workflow renames.** Several steps ran under Terra method configuration names that no longer match any workflow, and some workflows have been renamed since. Matches above were inferred from inputs, tasks and purpose. In particular `UpdateGenotypes`, used for these callsets and now archived, is superseded by [PostprocessCallset](../wdl/annotation_utils/PostprocessCallset.wdl), which performs the same operations behind its `run_*` flags, so `PostprocessCallset` is named for those steps.
- **`Whatshap` is still archived.** It haplotagged the reads for methylation profiling, but [HiPhase](../wdl/tools/HiPhase.wdl) now emits haplotagged BAMs, so it was deliberately retired rather than restored; section 2 links it in place.
- **QC and analysis workflows are excluded**, along with unrelated work sharing the same workspace (Paraphase, LPA and SMN1 assembly, Himito, Immuannot) and the depth-based CNV workflows.
