# Repository Structure

This document is a map of the repository: how the WDL workflows are organized, how Dockerfiles are built and consumed, what lives under `scripts/`, and how CI/CD is wired up. See [Conventions](conventions.md) for the detailed WDL/Python style rules, and [Workflows](workflows.md) for per-workflow input/output documentation.

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
notebooks/           # Ad hoc Jupyter notebooks (e.g. Terra cost analysis)
docs/                # Extended documentation (this file included)
.github/
  workflows/         # Active GitHub Actions CI
  scripts/           # Helper scripts invoked by CI
archive/             # Historical workflows, scripts, Dockerfiles, and their reference docs
```


## WDL Workflows
Workflows are split by role:
- **`wdl/annotation/`** - top-level annotation workflows, always prefixed `Annotate*`. Each characterizes one aspect of the callset (MEIs, functional consequence, external-database overlap, etc.) and typically outputs a TSV rather than a VCF.
- **`wdl/annotation_utils/`** - VCF manipulation utilities used to glue the annotation workflows together (splitting, merging, applying TSV annotations back onto a VCF, post-processing).
- **`wdl/tools/`** - thin wrappers around individual bioinformatics tools (PALMER, TRGT, HiPhase, mosdepth, etc.) that aren't annotation-specific.
- **`wdl/utils/`** - not directly run. `Structs.wdl` defines the shared `RuntimeAttr` struct; `Helpers.wdl` holds reusable tasks (subsetting, concatenation, sharding, etc.) imported by the other three directories. It also holds importable (but not directly dockstore-registered) sub-workflows used as building blocks by other workflows: `BedtoolsClosestSV.wdl`, `ExactMatch.wdl`, `ScatterVcf.wdl`, `TruvariMatch.wdl`.

Every workflow directly run in the pipeline (i.e. everything in `annotation/`, `annotation_utils/` and `tools/`) must have a matching entry in [`.dockstore.yml`](../.dockstore.yml), under its corresponding `# Annotation Workflows` / `# Annotation Utilities` / `# Tools` section. This is enforced by CI (see below). For the full WDL/task/input style convention, see [Conventions](conventions.md).


## Dockerfiles
Every container image is built from a `dockerfiles/Dockerfile.<image-name>` file, where `<image-name>` is **all-lowercase** and is exactly the name the image is pushed under - e.g. `Dockerfile.utils` builds `utils`, `Dockerfile.stranalysis` builds `stranalysis`. There is no separate name-mapping file; the Dockerfile suffix mechanically is the image name.

Images are pushed to Artifact Registry at:
```
us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/<image-name>
```

**Building and pushing** an image is done via the helper script:
```bash
dockerfiles/build_docker.sh <image-name>
# e.g. dockerfiles/build_docker.sh utils
```
This script:
- Resolves any `ARG` declared in the Dockerfile (with no inline default) from [`dockerfiles/versions.env`](../dockerfiles/versions.env), keyed as `<image-name>__<ARG_NAME>` - this file is the single source of truth for pinned tool/library versions.
- Auto-increments a `kj_V<N>` tag from the current highest version already pushed for that image, builds with `podman`, and pushes the image under both the new `kj_V<N>` tag and `:latest`.

**Which tag to use:** WDL tasks never hardcode a docker image URI - the `String docker` task input is always supplied by the caller (via Terra workspace data, per [Conventions](conventions.md)). By convention, workspace data should point at the `:latest` tag for each image, since every `build_docker.sh` run retags `:latest` to the newest build. The `kj_V<N>` tags exist purely as an immutable version history if you need to pin or roll back to a specific prior build.


## Scripts
`scripts/` holds standalone Python/R scripts that run inside a Docker container as CLI tools (not imported as a library across containers - each script is self-contained within the container that runs it):
- **`helper/`** - shared standalone utility scripts, including Hail helpers.
- **`mei/`** - PALMER call post-processing (raw-call-to-VCF conversion, annotation transfer).
- **`benchmark/`** - `bedtools closest` benchmarking scripts.
- **`annotation/`** - standalone annotation scripts not tied to a Hail/VEP or MEI-specific workflow (e.g. genomic context annotation).


## GitHub Actions / CI-CD
Three CI workflows currently run on push/PR to `main`, each gated on the paths they check:

| Workflow | Trigger paths | What it does |
|---|---|---|
| [`wdl-validation.yml`](../.github/workflows/wdl-validation.yml) | `wdl/**` | Runs `womtool validate` over every `.wdl` file. |
| [`python-linting.yaml`](../.github/workflows/python-linting.yaml) | `scripts/**` | Runs `flake8` over `scripts/`. |
| [`dockstore-sync.yml`](../.github/workflows/dockstore-sync.yml) | `wdl/**`, `.dockstore.yml` | Runs [`check_dockstore_sync.py`](../.github/scripts/check_dockstore_sync.py), which fails if any active workflow in `wdl/annotation`, `wdl/annotation_utils` or `wdl/tools` is missing a `.dockstore.yml` entry (or vice versa), aside from the allowlisted external-repo workflows `AnnotateAF` and `QcAnnotations`. |

A [`pyproject.toml`](../pyproject.toml) `[tool.black]` section pins `line-length = 88` for local formatting with `black`, but this is not run in CI - only `flake8` is enforced.

A Docker build/push CI pipeline (triggering on `dockerfiles/**` changes, authenticating to GCP via `google-github-actions/auth`) was also designed, but is intentionally **not active** - it would require a `GCP_SA_KEY` service-account secret to be configured in the repo, which hasn't been set up. Until then, images are built and pushed locally via `dockerfiles/build_docker.sh`.
