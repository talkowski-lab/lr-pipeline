# Processing Notes

## Callset Generation - V1
- _FillFormatFields_ on _allele_type_annotated_vcf_.
- _NormalizeDuplicationOrigins_.
- Annotation: _AnnotateCallsetOverlap_, _AnnotateDbVaR_, _AnnotateGQMetrics_, _AnnotateSQMetrics_, _AnnotateSVAnnotate_, (AoU Only)_AnnotateAgeMetrics_.
- _AnnotateVcf_Functional_.
- _AnnotateAF_.
- _AnnotateVcf_Downstream_.
- _PostprocessCallset_.
- Untrim variants:
	- _FindUntrimmedAlleles_.
	- Annotation: _AnnotateCallsetOverlap_, _AnnotateDbSNP_, _AnnotateDbVaR_, _AnnotateInSilicoPredictors_, _AnnotateVRS_.
	- _AnnotateVcfCleared_.

## Callset Generation - V2
- _NormalizeDuplicationOrigins_ on _allele_type_annotated_filled_vcf_.
- Annotation: _AnnotateCallsetOverlap_, _AnnotateSVAnnotate_
- _AnnotateVcf_Downstream_.
- _PostprocessCallset_.
- _AnnotateAF_.
- Untrim variants:
	- _FindUntrimmedAlleles_.
	- _AnnotateCallsetOverlap_.
	- _AnnotateVcfCleared_.
- (HPRC/HGSVC Only) _ResolveHaplotypeOverlaps_.
- (HPRC/HGSVC Only) _AnnotateAF_.
- _NormalizeAlleleTypes_.
- _PostprocessCallset-DropFilters_ --> hprc_hgsvc_vcf, aou_vcf.
- (HPRC/HGSVC Only) _PostprocessCallset-FilterAssemblySingletons_ --> hprc_hgsvc_vcf_V2.
- _FilterLowCoverageRegions_ --> hprc_hgsvc_vcf_V3, aou_vcf_V2.
- (AoU Only) _StripGenotypes_.
- (HPRC/HGSVC Only)_FillFormatFields_ --> hprc_hgsvc_vcf_V4.
- (HPRC/HGSVC Only) Annotation: _AnnotatSQMetrics_, _AnnotatGQMetrics_, _AnnotateCallsetOverlap_.
- (HPRC/HGSVC Only) _AnnotateVcfCleared_ --> hprc_hgsvc_vcf_V5.
