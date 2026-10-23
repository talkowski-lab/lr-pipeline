# Cohort
This document describes the sample cohorts used to build and validate the released callsets - HPRC and HGSVC, combined into one callset, and All of Us Phase 1 - and where each cohort's sample list differs from its public metadata. How each callset was produced is recorded in [Trace](trace.md).


## HPRC/HGSVC
- HPRC.
  - 232 total samples.
  - [Metadata File](https://github.com/human-pangenomics/hprc_intermediate_assembly/blob/main/data_tables/sample/hprc_release2_sample_metadata.csv).
    - 234 total samples.
    - Additionally includes the GRCh38 and CHM13 references.
    - The SV integration sample table and the DeepVariant callset also include HG03492.
- HGSVC.
  - 65 total samples.
  - [Metadata File](https://www.internationalgenome.org/data-portal/data-collection/hgsvc3).
    - 67 total samples.
    - Renames NA21487 to GM21487.
    - Duplicates NA19129 (also includes GM19129) and NA20355 (also includes GM20355).
    - Additionally includes GM19320.
    - Misses NA24385 (HG002 in this callset).
- Combined.
  - 292 total samples.
  - Overlapping samples: HG002, HG00733, HG02818, NA19036, NA19240.


## All of Us Phase 1
- 1027 total samples.
- Metadata file not publicly accessible.
- All samples are unrelated and of African ancestry.
