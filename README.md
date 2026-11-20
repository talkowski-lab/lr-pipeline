# LR Annotation
This repository contains all code for processing long-read sequencing cohorts end-to-end - from variant calling and callset integration through phasing and annotation. The pipeline was run on two cohorts: a combined HPRC/HGSVC set (292 samples) and All of Us Phase 1 (1027 samples).

The pipeline covers three broad stages. 
1. Callset Generation: Per-sample variant calls from DeepVariant (for SNVs/indels), a series of SV callers and TRGT (for tandem repeats) are integrated into cohort-level VCFs and filtered.
2. Callset Processing: A series of preprocessing steps convert these VCFs into the desired format. Then, variants are physically phased per sample using HiPhase, and these phase blocks are then linked via backbone phasing. 
3. Callset Annotation: A suite of annotation workflows characterizes each variant - calling mobile element insertions (MEIs) and deletions (MEDs) using PALMER, SVAN and L1ME-AID, identifying duplications and NUMTs, annotating functional effects via VEP and SVAnnotate, computing in-silico predictor scores, assigning allele frequencies, linking variants to external databases (dbSNP, dbVaR, gnomAD) and more.

**Stack:** All workflows are written in WDL and executed via Cromwell on Terra, which is built on top of GCP. The code logic is a mixture of Python, R, Bash and Hail. 


## Repository Layout
```
wdl/
  00_annotate_SV/                 # Annotations added to long-read variants for the gnomAD-LR resource (prefix: Annotate*)
  01_insertion_source_analyses/   # Origin of inserted sequences (AnnotateInsertionOrigins, AnnotateSmallInsertionOrigins)
  02_variant_TR_constraint/       # Per-sample burden, per-gene variant / TR constraint, CADD scores
  03_expression_alterations/      # Expression alterations (PhasedHetPerGene: phased hets per gene for Iso-Seq phasing)
  04_methylation_analyses/        # Methylation analyses (placeholder)
  05_CNV_depth_calling/           # Depth CNV calling, genomic-disorder breakpoints, read recurrence, local assembly
  06_compound_het/                # Compound-het / double-hit analyses (placeholder)
  07_ancestry_matching/           # PCA ancestry projection and ancestry-matched control selection
  annotation_utils/  # VCF manipulation and utility workflows
  tools/             # Individual bioinformatics tool wrappers
  utils/             # Shared structs (Structs.wdl) and helper tasks (Helpers.wdl)
scripts/
  vep/               # VEP and Hail-based annotation scripts
  mei/               # MEI analysis scripts
  benchmark/         # Benchmarking scripts
  miscellaneous/     # Other utility scripts
dockerfiles/         # Dockerfile.<ToolName> for each container
data/
  references/        # Reference genomes, catalogs, BED files
  base_vcfs/         # Test/base VCF files
  metadata/          # Sample metadata, pedigrees, ancestry
  qc_annotations/    # QC annotation outputs
docs/                # Extended documentation
lrGS_TR_analyses/    # Tandem-repeat variant (TRV) analyses and workflows (see its README)
archive/             # Deprecated workflows, scripts, and dockerfiles
```


## Documentation
- [Annotations](docs/annotations.md) - VCF INFO fields, FORMAT fields and filter definitions.
- [Cohort](docs/cohort.md) - sample cohorts, sizes and metadata sources.
- [Conventions](docs/conventions.md) - WDL, Python and codebase conventions.
- [Pipeline](docs/pipeline.md) - end-to-end pipeline description covering callset generation, preprocessing and annotation.
- [References](docs/references.md) - all reference files and their GCS locations.
- [Workflows](docs/workflows.md) - annotation workflows, annotation utilities and tool wrappers with inputs and outputs.
