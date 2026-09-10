# Processing Notes

## Callset Generation - V1
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

## Callset Generation - V2
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
- _AnnotateVcfPostHoc_ --> hprc_hgsvc_vcf_V4, aou_vcf_V3.
- (AoU Only) _StripGenotypes_ --> aou_sites_vcf.
- (HPRC/HGSVC Only)_FillFormatFields_ --> hprc_hgsvc_vcf_V5.
- (HPRC/HGSVC Only) _FilterLowCoverageGenotypes_ --> hprc_hgsvc_vcf_V6.
- (HPRC/HGSVC Only) _AnnotatSQMetrics_, _AnnotatGQMetrics_.
- (HPRC/HGSVC Only) _AnnotateVcfGQSQ_ --> hprc_hgsvc_vcf_V7.
- (HPRC/HGSVC Only) _AnnotateAFPostHoc_ --> hprc_hgsvc_vcf_V8.
