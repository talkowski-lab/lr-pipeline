# Getting Started
This document covers what you need to run the pipeline on a new cohort, from raw unaligned lrGS reads to an annotated VCF for each contig. The quickest route is to clone the template Terra workspace, which already holds every method config in run order and the workspace data they read, and replace its example cohort with your own. [Pipeline](pipeline.md) describes what each step does.


## Template Workspace
Clone [LR_GNOMAD-AoU_TALK_Pipeline-Template](https://app.terra.bio/#workspaces/LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_Pipeline-Template) into your own billing project. It runs the whole pipeline on a small example cohort, so every data table column is filled in for reference. The clone holds:
- **Method configs** - one per workflow run, numbered by its step in [Run Order](pipeline.md#run-order); configs sharing a number can be launched together. Unnumbered configs run the [Optional Analyses](pipeline.md#optional-analyses).
- **Workspace data** - the reference files in [References](references.md), the `*_docker` images in [Docker images](dockers.md), the cohort-wide inputs below, and the few cohort-wide outputs that later configs read back.
- **Data tables** - `sample`, one row per sample with its `unaligned_bam`; `sample_set`, the cohort; `contig`, one row per contig with its per-contig references, which every workflow from step 13 runs on; and `contig_set`, with `all_contigs`, every contig, for the configs that read all contigs at once, and `truth_contigs`, every contig but chrY, for `EvaluateCallsetAccuracy`.

Outside the lab, most reference files and the repository images are not publicly readable. Copy each reference from the Source listed in [References](references.md) into your own bucket, build the repository images with `dockerfiles/build_docker.sh` as described in [Docker images](dockers.md), and point the workspace data at your copies.


## Inputs
Everything below is cohort-specific. Reference files, repeat catalogs and Docker images are shared across cohorts - see [References](references.md) and [Docker images](dockers.md).
- **Unaligned lrGS reads** - one unaligned BAM per sample, plus the sample IDs: the `sample` table's `unaligned_bam` and `sample_id`.
- **PED** - six-column pedigree covering every sample, giving each sample's sex. Its family and parent columns are read only by the trio-based workflows, so a cohort of unrelated samples can leave them `0`. The `ped` workspace data.
- **Ancestry** - per-sample genetic similarity group assignments, for allele frequency annotation. The `ancestry` workspace data.
- **Cohort metadata TSV** - PED and ancestry combined, built by `CreateCohortMetadata` in step 01.
- **Phased truth base VCFs** _(optional)_ - haplotype-resolved VCFs covering the cohort's samples, which enable backbone phasing in step 14.
- **Sample age data** _(optional)_ - per-sample ages, which enable age-bin carrier annotation.


## Running a New Cohort
1. Clone the template workspace and, outside the lab, repoint its workspace data and the `contig` table's reference columns at your own copies.
2. Replace the `sample` rows with your samples and their unaligned BAMs, and the cohort set's members with their IDs.
3. Upload your PED and ancestry files and set the `ped` and `ancestry` workspace data to them.
4. Point `14-BackbonePhase` at your phased truth VCFs, or skip it and point every config that reads `backbone_phased_vcf` at `13-SplitVcfPerContig_Phased`'s `phased_vcf` column instead.
5. Launch the configs in step order, waiting for every config of a step to finish before starting the next. From step 13, per-contig configs are launched on all `contig` rows at once.

To run outside Terra, each workflow has an example input JSON in [`inputs/`](../inputs/), [Workflows](workflows.md) describes every input and output, and the smoke tests in [tests/README.md](../tests/README.md) check a local setup.


## CNV Pipeline
Depth-based CNV calling is a separate pipeline, starting from aligned BAMs, with its own template workspace, [LR_GNOMAD-AoU_TALK_CNV-Pipeline-Template](https://app.terra.bio/#workspaces/LR_GNOMAD_1_CO-AoU_TALK/LR_GNOMAD-AoU_TALK_CNV-Pipeline-Template). [CNV Pipeline](pipeline_cnvs.md) describes its steps and the choices to make before running it.
