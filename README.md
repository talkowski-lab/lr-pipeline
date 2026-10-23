# LR Pipeline
[![Validate Workflows](https://github.com/talkowski-lab/lr-pipeline/actions/workflows/wdl-validation.yml/badge.svg)](https://github.com/talkowski-lab/lr-pipeline/actions/workflows/wdl-validation.yml)
[![Tests](https://github.com/talkowski-lab/lr-pipeline/actions/workflows/tests.yml/badge.svg)](https://github.com/talkowski-lab/lr-pipeline/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Dockstore](https://img.shields.io/badge/Dockstore-lr--pipeline-blue)](https://dockstore.org/search?entryType=workflows&search=lr-pipeline)

This repository contains code for processing cohorts sequenced with long-read genomic sequencing methods - from variant calling and callset integration, all the way through to phasing and annotation. The pipeline was run on two cohorts as part of its initial release: a combined set of 292 samples from the HPRC and HGSVC consortia, and 1027 samples from Phase 1 of the All of Us Research Program - see [Cohort](docs/cohort.md).

The pipeline covers three broad stages.
1. Callset Generation: Per-sample variant calls from DeepVariant (for SNVs/indels), a series of SV callers and TRGT (for tandem repeats) are integrated into cohort-level VCFs and filtered.
2. Callset Processing: A series of preprocessing steps convert these VCFs into the desired format. Then, variants are physically phased per sample using HiPhase, and these phase blocks are then linked via backbone phasing.
3. Callset Annotation: A suite of annotation workflows characterizes each variant - calling mobile element insertions (MEIs) and deletions (MEDs) using PALMER, SVAN and L1ME-AID, identifying duplications and NUMTs, annotating functional effects via VEP and SVAnnotate, computing in-silico predictor scores, assigning allele frequencies, linking variants to external databases (dbSNP, dbVar, gnomAD) and more.

**Stack:** All workflows are written in WDL 1.0 and executed via Cromwell on Terra, which is built on top of GCP. The code logic is a mixture of Python, R, Bash and Hail.


## Scope
The pipeline targets the GRCh38 reference and PacBio HiFi long reads. Cromwell on Terra, on GCP, is the supported backend; other backends and references are untested.


## Quickstart
1. Import the workflows you need into Terra from [Dockstore](https://dockstore.org/search?entryType=workflows&search=lr-pipeline), where each is registered as `github.com/talkowski-lab/lr-pipeline/<WorkflowName>`; use the `main` version.
2. Set the workspace data for every `*_docker` input as described in [Docker images](docs/dockers.md), and for the reference files as described in [References](docs/references.md), fetching each from its public source where the listed copy is not publicly readable. The images in the maintainers' registry are private, so build your own with `dockerfiles/build_docker.sh`.
3. Follow [Pipeline](docs/pipeline.md) to run the stages end to end on a new cohort. Each workflow runs on its own, with its outputs passed as inputs to the next, and [Workflows](docs/workflows.md) describes every input and output.
4. Example input JSONs with standard placeholder values for every workflow are in [`inputs/`](inputs/), one per workflow.
5. To check a local setup, run the smoke tests described in [tests/README.md](tests/README.md).


## Documentation
Running the pipeline:
- [Pipeline](docs/pipeline.md) - how to run the pipeline end-to-end on a new cohort, covering callset generation, phasing and annotation.
- [CNV Pipeline](docs/pipeline_cnvs.md) - how to call depth-based CNVs on a new cohort, from per-sample depth to a genotyped cohort CNV VCF.
- [References](docs/references.md) - all reference files, their sources and their GCS locations.
- [Docker images](docs/dockers.md) - every image the pipeline uses, and how repository images are built, tagged and pushed.

Reference:
- [Workflows](docs/workflows.md) - every workflow and shared sub-workflow, with its inputs and outputs.
- [Annotations](docs/annotations.md) - VCF FILTER values, INFO fields and FORMAT fields.

Released callsets:
- [Cohort](docs/cohort.md) - sample cohorts, sizes and metadata sources.
- [Trace](docs/trace.md) - the steps actually run to produce the HPRC/HGSVC and All of Us callsets, with data table columns and cohort divergences.

Development:
- [Repository Structure](docs/repository-structure.md) - directory layout, the scripts inventory, and where each kind of code lives.
- [Conventions](docs/conventions.md) - WDL, Python and Markdown style conventions, and the branch workflow.
- [CI/CD](docs/ci-cd.md) - GitHub Actions checks, the matching local pre-flight commands, and Dockstore registration.


## Contributing
Contributions are welcome through a fork and a pull request; see [CONTRIBUTING.md](CONTRIBUTING.md) for the local checks and the review requirements, and the [Code of Conduct](CODE_OF_CONDUCT.md).
- Style rules for WDL, Python and Markdown are in [Conventions](docs/conventions.md), and the checklist for adding a workflow or tool is in [AGENTS.md](AGENTS.md).
- Workflows are imported into Terra from Dockstore via [`.dockstore.yml`](.dockstore.yml); every directly-run workflow needs an entry there, and listing a branch in an entry is how a workflow is tested from that branch - see [Dockstore Registration](docs/ci-cd.md#dockstore-registration).


## Citation
If you use this pipeline, please cite it using [CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button reads the same file.


## License
Released under the [MIT License](LICENSE). Report security problems as described in [SECURITY.md](SECURITY.md).
