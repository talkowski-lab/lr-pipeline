# Repository Structure
This document is a map of the repository - how the WDL workflows are organized, where Dockerfiles and scripts live, and what each top-level directory holds. Image builds are covered in [Docker images](dockers.md) and the checks that gate `main` in [CI/CD](ci-cd.md).


## File Structure
```
wdl/
  annotation/        # Main annotation workflows (prefix: Annotate*)
  annotation_utils/  # VCF manipulation and utility workflows
  tools/             # Individual bioinformatics tool wrappers
  utils/             # Shared structs (Structs.wdl) and helper tasks (Helpers.wdl)
scripts/
  annotation/        # Standalone annotation scripts (e.g. genomic context)
  helper/            # Shared standalone utility scripts, including Hail helpers
  mei/               # MEI analysis scripts
  benchmark/         # Benchmarking scripts
dockerfiles/         # Dockerfile.<image-name> (lowercase) for each container, plus build_docker.sh and versions.env
data/                # Local-only analysis inputs, outputs and scratch work; gitignored and never referenced by WDL
docs/                # Extended documentation (this file included)
.github/
  workflows/         # Active GitHub Actions CI
  scripts/           # CI check scripts, also run locally, plus merge_branch.sh for the branch merge sequence
archive/             # Retired workflows, scripts, Dockerfiles, notebooks and their reference docs
  docs/workflows.md  # Generated from the retired workflows, the counterpart to docs/workflows.md
  notebooks/         # Retired ad hoc Jupyter notebooks (Terra cost analysis via the API and BigQuery)
.dockstore.yml       # Dockstore registration for every directly-run workflow
AGENTS.md            # Canonical agent instructions; .claude/CLAUDE.md and .github/copilot-instructions.md are symlinks to it
```


## WDL Workflows
Workflows are split by role:
- **`wdl/annotation/`** - top-level annotation workflows, always prefixed `Annotate*`. Each characterizes one aspect of the callset (MEIs, functional consequence, external-database overlap, etc.) and typically outputs a TSV rather than a VCF.
- **`wdl/annotation_utils/`** - VCF manipulation utilities used to glue the annotation workflows together (splitting, merging, applying TSV annotations back onto a VCF, post-processing).
- **`wdl/tools/`** - thin wrappers around individual bioinformatics tools (PALMER, TRGT, HiPhase, mosdepth, etc.) that aren't annotation-specific.
- **`wdl/utils/`** - not directly run. `Structs.wdl` defines the shared `RuntimeAttr` struct; `Helpers.wdl` is a task library holding reusable tasks (subsetting, concatenation, sharding, etc.) imported by the other three directories, with its tasks declared alphabetically. It also holds importable sub-workflows, which are never Dockstore-registered: `BedtoolsClosestSV.wdl`, `ExactMatch.wdl`, `ScatterVcf.wdl` and `TruvariMatch.wdl` are callset-matching and sharding building blocks, while `LRCNVs.wdl`, `DepthPreprocessing.wdl`, `DepthClustering.wdl` and `GenotypeDepth.wdl` form the depth-based CNV pipeline driven by `tools/LongReadCNVs.wdl`. All eight are described in [Sub-workflows](workflows.md#sub-workflows).

Every workflow directly run in the pipeline (i.e. everything in `annotation/`, `annotation_utils/` and `tools/`) must have a matching entry in [`.dockstore.yml`](../.dockstore.yml), under its corresponding `# Annotation Workflows` / `# Annotation Utilities` / `# Tools` section. This is enforced by CI (see [CI/CD](ci-cd.md)), and that file's ordering is also the section order of [Workflows](workflows.md). For the full WDL/task/input style convention, see [Conventions](conventions.md).

Every workflow also documents itself: a `meta` block holds its description and a `parameter_meta` block describes every input and output, and [Workflows](workflows.md) is generated from those blocks rather than written by hand. See [Workflow documentation](conventions.md#workflow-documentation).


## Dockerfiles
Every container image is built from a `dockerfiles/Dockerfile.<image-name>` file, where `<image-name>` is **all-lowercase** and is, with one documented exception, exactly the name the image is pushed under - e.g. `Dockerfile.utils` builds `utils`. Pinned tool and library versions live in [`dockerfiles/versions.env`](../dockerfiles/versions.env).

The full image inventory, the registry path, and the build, tag and push process via `build_docker.sh` are documented in [Docker images](dockers.md).


## Scripts
`scripts/` holds standalone Python and R scripts that run inside a Docker container as CLI tools. They are not imported as a library across containers - each script is self-contained within the container that runs it.

Scripts reach a running task by one of two routes. Most are baked into the `utils` image, which copies the whole tree to `/opt/scripts` (see `COPY ./scripts /opt/scripts` in `Dockerfile.utils`), so a change to one of them only takes effect once `utils` is rebuilt and every image inheriting from it is refreshed. The three Hail scripts are instead fetched at run time over HTTPS from the `main` branch of this repository, so changes to those go live on push without any rebuild.

| Script | Called by | Purpose |
|---|---|---|
| `annotation/annotate_genomic_context.R` | AnnotateRegion | Assign genomic context labels to variants from breakpoint and CNV body coverage tables. |
| `annotation/annotate_insilico_predictors.py` | AnnotateInSilicoPredictors, PostProcessTRLociAoU, PostProcessTRLociHPRCHGSVC (fetched by URL) | Annotate a VCF with gnomAD in-silico predictor scores (CADD, REVEL, phyloP, SpliceAI, Pangolin) from Hail tables. |
| `benchmark/R1.bedtools_closest_CNV.R` | BedtoolsClosestSV | Compare CNVs against their closest reference callset records. |
| `benchmark/R2.bedtools_closest_INS.R` | BedtoolsClosestSV | Compare insertions against their closest reference callset records. |
| `helper/estimatePloidy.R` | CreateCohortDepthFiles | Estimate per-sample ploidy from binned coverage. Vendored unmodified from GATK-SV `src/WGD/bin/`; preserve its license header. |
| `helper/estimated_CN_denoising.py` | CreateCohortDepthFiles | Denoise estimated copy number across samples. Vendored from GATK-SV; preserve its license header. |
| `helper/medianCoverage.R` | CreateCohortDepthFiles | Calculate median bin coverage per sample from a bincov matrix. Loads the matrix into memory, so split large inputs by contig. |
| `helper/split_vcf_hail.py` | AnnotateVEPHail, ScatterVcf (fetched by URL) | Split a VCF into a fixed number of shards using Hail. |
| `helper/vep_annotate_hail.py` | AnnotateVEPHail (fetched by URL) | Run VEP over a VCF via Hail and emit the annotated result. |
| `mei/PALMER_to_vcf.py` | Helpers.wdl | Convert raw PALMER mobile-element calls into a VCF. |
| `mei/PALMER_transfer_annotations.py` | AnnotatePALMER | Transfer PALMER annotations onto matching callset records, using edit distance to resolve near-matches. |

Some workflows invoke scripts under `/opt/sv-pipeline/` or `/app/SVAN/`. Those belong to the GATK-SV and SVAN images respectively and are not part of this repository.


## CI-CD
Seven GitHub Actions workflows gate `main`, covering WDL syntax, WDL style, Markdown style, Python linting, Dockstore registration, agent-instruction sync and workflow-documentation generation. See [CI/CD](ci-cd.md) for the full table, the matching commands to run locally before pushing, and how Dockstore registration and image releases work.
