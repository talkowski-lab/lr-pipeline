# Pipeline 
This document describes how the gnomAD long-read callsets were generated end-to-end, covering three major phases: callset generation, callset processing phasing and callset annotation.


---


## 1. Callset Generation
Raw reads were aligned to GRCh38 at 65x coverage using Minimap2, then downsampled to 30x with samtools. Three variant callsets were produced in parallel from the downsampled reads.

### SNV/Indel Callset
1. **DeepVariant**: Call short variants per sample from the downsampled reads.
2. **GLNexus** : Joint call and integrate per-sample DeepVariant gVCFs into a cohort short-variant VCF.

### SV Callset
1. **SV Callers**: Multiple callers (cuteSV, Sniffles, Delly, pbsv, Sawfish, dipcall, hapdiff) were run per sample to generate structural variant calls.
2. **Integration**: Per-sample SV calls were merged using Truvari and then filtered via differing mthods depending on the cohort:
   - HPRC/HGSVC: ML-based filtering.
   - All of Us: [Kanpig](../wdl/tools/Kanpig.wdl) ([docs](workflows.md#kanpig)) was used to regenotype across all sites, with non-ref support needed to retain a call.

### TRGT Callset
**[TRGT](../wdl/tools/TRGT.wdl)** ([docs](workflows.md#trgt)): Genotyped tandem repeat loci per sample from the diploid reads, running separately against two catalogs:
   - TRExplorer v1.0.1.
   - Vamos v2.1.


---


## 2. Callset Processing
The three callsets were integrated and phased. Steps were run as described below.

### VCF Integration
1. **[IntegrateVcfs](../wdl/annotation_utils/IntegrateVcfs.wdl)** ([docs](workflows.md#integratevcfs)): Normalized multiallelics, flagged short variants ≥ 50bp, added core INFO fields (`allele_type`, `allele_length`), renamed variant IDs and merged the SNV/indel and SV VCFs into a single cohort integrated VCF.
2. **[ExtractSampleVcfs](../wdl/annotation_utils/ExtractSampleVcfs.wdl)** ([docs](workflows.md#extractsamplevcfs)): Split the cohort VCF into per-sample SNV/indel and SV VCFs for phasing.

### TR Callset Integration
Concurrently with VCF integration, the per-sample TRGT callsets were processed:
- **[CombineTRs](../wdl/annotation_utils/CombineTRs.wdl)** ([docs](workflows.md#combinetrs)) - deduplicated overlapping loci across the TRExplorer and Vamos catalogs per sample.
- **[AddEndTRs](../wdl/annotation_utils/AddEndTRs.wdl)** ([docs](workflows.md#addendtrs)) - added the `END` INFO field required by `trgt merge` downstream.

### Phasing
HiPhase was run twice per sample - once with the TRGT VCF and once without - then both runs were merged:
1. **[HiPhase](../wdl/tools/HiPhase.wdl)** ([docs](workflows.md#hiphase)) - jointly phased per-sample SNV/indel, SV and (for the TRGT run) TRGT VCFs against the aligned reads.
2. **[HiPhaseMerge](../wdl/tools/HiPhaseMerge.wdl)** ([docs](workflows.md#hiphasemerge)) - merged per-sample phased VCFs into a cohort phased VCF; when `merge_trgt` is enabled, also merged TRGT calls separately.
3. **[FillPhasedGenotypes](../wdl/annotation_utils/FillPhasedGenotypes.wdl)** ([docs](workflows.md#fillphasedgenotypes)) - populated the cohort integrated VCF with phased genotypes from the HiPhase output.
4. **[IntegrateTRs](../wdl/annotation_utils/IntegrateTRs.wdl)** ([docs](workflows.md#integratetrs)) - combined the TRGT VCF with the base phased VCF, flagging overlapping variants.

### Backbone Phasing
Backbone phasing was used to extend HiPhase phase blocks by linking them using a phased truth set - the HGSVC2024v1.0 combined truth VCF and the HPRC v2.0 pangenome ("wave") VCF:
1. **[BackbonePhase](../wdl/tools/BackbonePhase.wdl)** ([docs](workflows.md#backbonephase)) - transferred phasing from a set of truth base VCFs onto the target cohort VCF to connect phase blocks.
2. **[MergeBackbonePhased](../wdl/annotation_utils/MergeBackbonePhased.wdl)** ([docs](workflows.md#mergebackbonephased)) - merged the backbone-phased VCFs into the final cohort phased VCF.

### TRGT Downstream Processing
Concurrently with phasing, the merged TRGT VCF was processed for downstream use:
- **[TRGTLPS](../wdl/tools/TRGTLPS.wdl)** ([docs](workflows.md#trgtlps)) - computed the longest polymer sequence (LPS) per TRGT locus per sample; the resulting TSV is used downstream in allele frequency annotation.
- **[GenerateTRGTJson](../wdl/annotation_utils/GenerateTRGTJson.wdl)** ([docs](workflows.md#generatetrgtjson)) - generated per-locus allele-frequency histograms for the TR browser.


---


## 3. Callset Annotation
The phased cohort VCF was annotated in two broad stages: a set of parallel variant characterization workflows, followed by downstream integration and post-processing.

### Variant Characterization (parallel)
The following workflows were run concurrently on the phased cohort VCF:
- **[RepeatMasker](../wdl/tools/RepeatMasker.wdl)** ([docs](workflows.md#repeatmasker)) - annotated mobile-element content in insertions; used as input to AnnotatePALMER.
- **[AnnotateL1MEAID](../wdl/annotation/AnnotateL1MEAID.wdl)** ([docs](workflows.md#annotatel1meaid)) - identified MEI calls using L1ME-AID and INTACT_MEI.
- **[AnnotateSVAN](../wdl/annotation/AnnotateSVAN.wdl)** ([docs](workflows.md#annotatesvan)) - annotated MEIs, MEDs, tandem duplications, dispersed duplications and NUMTs using SVAN.
- **[AnnotateMEDs](../wdl/annotation/AnnotateMEDs.wdl)** ([docs](workflows.md#annotatemeds)) - identified mobile element deletions by intersecting deletions against the MEI catalog.
- **[AnnotateIndelTRs](../wdl/annotation/AnnotateIndelTRs.wdl)** ([docs](workflows.md#annotateindeltrs)) - flagged short indels that represent tandem repeats.
- **[AnnotateSQMetrics](../wdl/annotation/AnnotateSQMetrics.wdl)** ([docs](workflows.md#annotatesqmetrics)) - recomputed site-level quality metrics (HWE, inbreeding coefficient, AS fields).
- **[AnnotateGQMetrics](../wdl/annotation/AnnotateGQMetrics.wdl)** ([docs](workflows.md#annotategqmetrics)) - computed binned distributions of genotype quality and allele balance.
- **[AnnotateAgeMetrics](../wdl/annotation/AnnotateAgeMetrics.wdl)** ([docs](workflows.md#annotateagemetrics)) _(All of Us only)_ - computed age-bin carrier distributions per variant.
- **[AnnotateRegion](../wdl/annotation/AnnotateRegion.wdl)** ([docs](workflows.md#annotateregion)) - assigned each variant a genomic region class (SR, SD, RM or US).
- **[AnnotateInSilicoPredictors](../wdl/annotation/AnnotateInSilicoPredictors.wdl)** ([docs](workflows.md#annotateinsilicopredictors)) - looked up CADD, Pangolin, PhyloP, REVEL and SpliceAI scores for SNVs and indels.
- **[AnnotateGnomADSTR](../wdl/annotation/AnnotateGnomADSTR.wdl)** ([docs](workflows.md#annotategnomadstr)) - matched tandem repeat calls against the gnomAD V4 TR catalog.
- **[AnnotateDbSNP](../wdl/annotation/AnnotateDbSNP.wdl)** ([docs](workflows.md#annotatedbsnp)) - assigned dbSNP rsIDs to matched variants.
- **[AnnotateVRS](../wdl/annotation/AnnotateVRS.wdl)** ([docs](workflows.md#annotatevrs)) - annotated each variant with GA4GH VRS identifiers.

Additionally, the Cohort PALMER VCF was processed in parallel:

- **[AnnotatePALMER](../wdl/annotation/AnnotatePALMER.wdl)** ([docs](workflows.md#annotatepalmer)) - matched MEI calls from PALMER against the cohort VCF insertions using the RepeatMasker output.

### MEI Consolidation and Allele Typing
1. **[AnnotateMEIs](../wdl/annotation/AnnotateMEIs.wdl)** ([docs](workflows.md#annotatemeis)) - consolidated the L1MEAID, PALMER and SVAN MEI annotations into a single harmonized MEI TSV.
2. **[SubsetTsvToColumns](../wdl/annotation_utils/SubsetTsvToColumns.wdl)** ([docs](workflows.md#subsettsvtocolumns)) - extracted the duplication-specific columns from the SVAN output.
3. **[AnnotateAlleleType](../wdl/annotation_utils/AnnotateAlleleType.wdl)** ([docs](workflows.md#annotatealleletype)) - updated the `allele_type` INFO field using the MED, MEI and duplication annotation TSVs.

### Functional and Overlap Annotation (parallel)
The following workflows were run concurrently on the allele-type-annotated VCF:
- **[AnnotateVEPHail](../wdl/annotation/AnnotateVEPHail.wdl)** ([docs](workflows.md#annotatevephail)) - predicted functional effects using VEP via Hail.
- **[AnnotateSVAnnotate](../wdl/annotation/AnnotateSVAnnotate.wdl)** ([docs](workflows.md#annotatesvannotate)) - predicted functional effects for SVs using GATK SVAnnotate.
- **[AnnotateCallsetOverlap](../wdl/annotation/AnnotateCallsetOverlap.wdl)** ([docs](workflows.md#annotatecallsetoverlap)) - identified variants overlapping gnomAD V4 using exact, Truvari and bedtools closest matching.
- **[AnnotateDbVaR](../wdl/annotation/AnnotateDbVaR.wdl)** ([docs](workflows.md#annotatedbvar)) - matched SVs against dbVar records.

### Downstream Integration and Post-Processing
1. **[AnnotateVcf](../wdl/annotation_utils/AnnotateVcf.wdl)** ([docs](workflows.md#annotatevcf)) _(AnnotateVcf_Downstream)_ - integrated all annotation TSVs into the VCF as INFO fields.
2. **[PostProcess](../wdl/annotation_utils/PostProcess.wdl)** ([docs](workflows.md#postprocess)) - applied final genotype updates using phased genotypes from the backbone-phased VCF, normalized ploidy, pruned MEIs, flagged homopolymer TRs and filtered singletons.
3. **[AnnotateAF](https://github.com/broadinstitute/gatk-sv/blob/kj_project_gnomad_lr/wdl/AnnotateAF.wdl)** ([docs](workflows.md#annotateaf)) - annotated allele frequencies using sample ancestries and the TRGT LPS TSV.
4. **[TransformDuplications](../wdl/annotation_utils/TransformDuplications.wdl)** ([docs](workflows.md#transformduplications)) - consolidated `allele_type` values for duplication variants and classified their `dup_type` as tandem, interspersed, complex or inv.
