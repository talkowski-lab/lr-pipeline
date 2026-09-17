# Processing Notes
This document tracks ad hoc processing notes across generated callsets.


## Phased `_V2` callset lineage
- Purpose: rerun _HiPhase_ without TRGT homopolymers (remove one-base repeat motifs before phasing).
- _HiPhase_.
- _HiPhaseMerge_ --> hiphase_merged_integrated_vcf_V2, hiphase_merged_trgt_vcf_V2.
- _FillPhasedGenotypes_ --> hiphase_phased_integrated_vcf_V2.
- _AnnotateTRs_ --> tr_annotated_vcf_V2.
- _SplitVcfPerContig_.
- _BackbonePhase_.
- _MergeBackbonePhased_.
- _PostprocessCallset_.


## V1
- _FillFormatFields_ on _allele_type_annotated_vcf_.
- _NormalizeDuplicationOrigins_.
- Annotation: _AnnotateCallsetOverlap_, _AnnotateDbVaR_, _AnnotateGQMetrics_, _AnnotateSQMetrics_, _AnnotateSVAnnotate_, (AoU Only)_AnnotateAgeMetrics_.
- _AnnotateVcfFunctional_.
- _AnnotateAF_.
- _AnnotateVcfDownstream_.
- _PostprocessCallset_.
- Untrim variants:
	- _FindUntrimmedAlleles_.
	- Annotation: _AnnotateCallsetOverlap_, _AnnotateDbSNP_, _AnnotateDbVaR_, _AnnotateInSilicoPredictors_, _AnnotateVRS_.
	- _AnnotateVcfCleared_.


## V2
- _NormalizeDuplicationOrigins_ on _allele_type_annotated_filled_vcf_.
- Annotation: _AnnotateCallsetOverlap_, _AnnotateSVAnnotate_
- _AnnotateVcfDownstream_.
- _PostprocessCallset_.
- _AnnotateAF_.
- Untrim variants:
	- _FindUntrimmedAlleles_.
	- _AnnotateCallsetOverlap_.
	- _AnnotateVcfCleared_.
- (HPRC/HGSVC Only) _ResolveHaplotypeOverlaps_.
- (HPRC/HGSVC Only) _AnnotateAFPostHoc_.
- _NormalizeAlleleTypes_.
- _PostprocessCallset-DropFilters_ --> hprc_hgsvc_vcf_V1, aou_vcf_V1.
- (HPRC/HGSVC Only) _PostprocessCallset-FilterAssemblySingletons_ --> hprc_hgsvc_vcf_V2.
- _FilterLowCoverageRegions_ --> hprc_hgsvc_vcf_V3, aou_vcf_V2.
- _AnnotateCallsetOverlap_.
- _AnnotateVcfCallsetOverlap_ --> hprc_hgsvc_vcf_V4, aou_vcf_V3.
- (HPRC/HGSVC Only) _PostProcessTRLociHPRCHGSVC_ --> hprc_hgsvc_vcf_V5.
- (AoU Only) _PostProcessTRLociAoU_ --> aou_vcf_V4.
- (AoU Only) _StripGenotypes_ --> aou_sites_vcf.
- (HPRC/HGSVC Only) _FillFormatFields_ --> hprc_hgsvc_vcf_V6.
- (HPRC/HGSVC Only) _FilterLowCoverageGenotypes_ --> hprc_hgsvc_vcf_V7.
- (HPRC/HGSVC Only) _AnnotatSQMetrics_, _AnnotatGQMetrics_.
- (HPRC/HGSVC Only) _AnnotateVcfGQSQ_ --> hprc_hgsvc_vcf_V8.
- (HPRC/HGSVC Only) _AnnotateAFPostHoc_ --> hprc_hgsvc_vcf_V9.
