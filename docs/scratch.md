# Processing Notes

## Callset Regeneration - V1
- _FillFormatFields_ on _allele_type_annotated_vcf_.
- _ParseAbsoluteOrigin_.
- Annotation: (AoU Only) _AnnotateAgeMetrics_, _AnnotateCallsetOverlap_, _AnnotateDbVaR_, _AnnotateGQMetrics_, _AnnotateSQMetrics_, _AnnotateSVAnnotate_
- _AnnotateVcf_Functional_.
- _AnnotateAF_.
- _AnnotateVcf_Downstream_.
- _PostProcess_.
- Untrim variants:
	- _FindUntrimmedAlleles_.
	- Annotation: _AnnotateCallsetOverlap_, _AnnotateDbSNP_, _AnnotateDbVaR_, _AnnotateInSilicoPredictors_, _AnnotateVRS_.
	- _AnnotateVcfCleared_.

## Callset Regeneration - V2
- _ParseAbsoluteOrigin_ on _allele_type_annotated_filled_vcf_.
- Annotation: _AnnotateCallsetOverlap_, _AnnotateSVAnnotate_
- _AnnotateVcf_Downstream_.
- _PostProcess_.
- _AnnotateAF_.
- Untrim variants:
	- _FindUntrimmedAlleles_.
	- _AnnotateCallsetOverlap_.
	- _AnnotateVcfCleared_.
- (HPRC/HGSVC Only) _ResolveHaplotypeOverlaps_.
- (HPRC/HGSVC Only) _AnnotateAF_.
- _TransformAlleleType_.
- (AoU Only) _DropGenotypes_.
