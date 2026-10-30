# Getting Started
This document covers what you need to run the pipeline on a new cohort, from raw unaligned lrGS reads to an annotated VCF for each contig. The quickest route is to clone the template Terra workspace, which already holds every method config in run order, along with the reference and Docker image workspace data they read, and then replace its example cohort with your own. [Pipeline](pipeline.md) describes what each step does.


## Template Workspace
Clone [LR_GNOMAD-AoU_TALK_Pipeline-Template](https://app.terra.bio/#workspaces/LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_Pipeline-Template) into your own billing project. It runs the whole pipeline on three CEPH 1463 Platinum Pedigree samples - NA12877 and NA12878 and their child NA12882 - so every data table column is filled in for reference. The clone holds:
- **Method configs** - one per workflow run, numbered by the step it runs in, as listed in [Run Order](pipeline.md#run-order). Configs sharing a number can be launched together once the previous step finishes. `SplitVcfPerContig` runs twice, as `13-SplitVcfPerContig_Phased` and `13-SplitVcfPerContig_PALMER`, `AnnotateVcf` twice, as `20-AnnotateVcf_Base` and `25-AnnotateVcf_Quality`, and `TRGT` once per repeat catalog.
- **Workspace data** - every reference file in [References](references.md) apart from the per-contig ones on the `contig` table, every `*_docker` image in [Docker images](dockers.md), and the cohort-wide inputs below.
- **Data tables**:
  - `sample` - one row per sample, keyed by `sample_id`, with its `unaligned_bam`. Per-sample workflows read from and write to these rows.
  - `sample_set` - the `cohort` set, listing the samples to run. Per-cohort workflows run on this row.
  - `contig` - one row per contig, keyed by `contig_id`, with a `sample_set` column linking to the cohort row, which lets contig-row configs read cohort-wide outputs as `this.sample_set.<column>`, and the contig's per-contig references - `gnomad_snv_indel_vcf`, `dbsnp_vcf` and `dbvar_vcf`, each with its `_idx`. Every workflow from step 13 runs on these rows, apart from the methylation branch, so the pipeline ends with one VCF per contig.
  - `contig_set` - `all_contigs`, every contig row, which `15-ConcatenateVcfsAcrossContigs` runs on to rebuild the cohort VCF for methylation.

Outside the lab, most reference files and the repository images are not publicly readable. Copy each reference from the Source listed in [References](references.md) into your own bucket, build the repository images with `dockerfiles/build_docker.sh` as described in [Docker images](dockers.md), and point the workspace data at your copies.


## Inputs
Everything below is cohort-specific and must be supplied per run. Reference files, repeat catalogs and Docker images are shared across cohorts and supplied separately - see [References](references.md) and [Docker images](dockers.md).
- **Unaligned lrGS reads** - one unaligned BAM per sample, plus the sample IDs. On Terra these are the `sample` table's `unaligned_bam` and `sample_id` columns.
- **PED** - six-column pedigree covering every sample, providing each sample's sex (`1` male, `2` female). Its family and parent columns are read only by the trio-based workflows - `EvaluateNCRCutoffs` in step 27 and the inheritance summaries in [Optional Analyses](pipeline.md#optional-analyses) - so a cohort of unrelated samples can leave them `0`. Required by the coverage, post-processing and allele frequency steps. On Terra this is the `ped` workspace data.
- **Ancestry** - per-sample genetic similarity group assignments, used for allele frequency annotation. On Terra this is the `ancestry` workspace data.
- **Cohort metadata TSV** - PED and ancestry combined into the single table the TR histograms consume. [CreateCohortMetadata](../wdl/annotation_utils/CreateCohortMetadata.wdl) builds it in step 01 if you only have the two files separately.
- **Phased truth base VCFs** _(optional)_ - haplotype-resolved VCFs covering the cohort's samples. Supplying these enables backbone phasing in step 14; the template uses `truth_platinum_pedigree_vcf`.
- **Sample age data** _(optional)_ - per-sample ages. Supplying this enables age-bin carrier annotation; the template has none, so it does not run `AnnotateAgeMetrics`.


## Running a New Cohort
1. Clone the template workspace and, if you work outside the lab, repoint its workspace data and the `contig` table's reference columns at your own reference and image copies.
2. Replace the `sample` rows with your samples and their unaligned BAMs, and the `cohort` set's members with their IDs.
3. Upload your PED and ancestry files and set the `ped` and `ancestry` workspace data to them.
4. Point `14-BackbonePhase` at your phased truth VCFs, or skip it and point every config that reads `backbone_phased_vcf` at `13-SplitVcfPerContig_Phased`'s `phased_vcf` column instead.
5. Launch the configs in step order, waiting for every config of a step to finish before starting the next. Step 13 onward is launched on the `contig` rows, all 24 at once, except the methylation branch: `15-ConcatenateVcfsAcrossContigs` runs on the `all_contigs` set, `16-ExtractSampleVcfs_Backbone` and `17-MethylationProfiling` on the samples, and `18-CreateCohortMethylationFile` on the cohort.

To run outside Terra, each workflow has an example input JSON with placeholder values in [`inputs/`](../inputs/), and [Workflows](workflows.md) describes every input and output. The smoke tests in [tests/README.md](../tests/README.md) check a local setup.
