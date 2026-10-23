# Docker images
This document lists the Docker images used by the pipeline - repo-built, collaborator-built, and published third-party - and describes how the repository's own images are built, tagged and pushed.


## Building and Pushing
Every container image is built from a `dockerfiles/Dockerfile.<image-name>` file, where `<image-name>` is all-lowercase and is exactly the name the image is pushed under - e.g. `Dockerfile.utils` builds `utils`, `Dockerfile.stranalysis` builds `stranalysis`. There is no separate name-mapping file; the Dockerfile suffix mechanically is the image name. The one exception is `Dockerfile.trgtlps`, which is pushed as `trgt-lps` with a hyphen, so do not assume the two always match.

Images are pushed to Artifact Registry at:
```
us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/<image-name>
```

This registry is private. To run the pipeline elsewhere, build the images you need from [`dockerfiles/`](../dockerfiles/), push them to a registry your Terra workspace can pull from, and point the workspace data at them.

Build and push an image with the helper script, which takes the image name as its only argument:
```bash
dockerfiles/build_docker.sh [--no-push] <image-name>
# e.g. dockerfiles/build_docker.sh utils
```

The script does the following:
- Resolves the build arguments through [`dockerfiles/build_args.sh`](../dockerfiles/build_args.sh) `<image-name>`, which prints a `NAME=value` line for every `ARG` declared in the Dockerfile without an inline default, taken from [`dockerfiles/versions.env`](../dockerfiles/versions.env) where it is keyed as `<image-name>__<ARG_NAME>`. CI uses the same script, so local and CI builds resolve identically. A missing key fails the build.
- Builds with `podman build --platform linux/amd64` from the repository root - meaning the worktree it is run from, so each branch builds its own tree.
- Tags the build by the branch it was built on, since branches are developed in parallel worktrees - see [Branches](conventions.md#branches).
- With `--no-push`, builds and tags locally only, as `localhost/<image-name>:<tag>`, and pushes nothing. This is how contributors without registry access build images, and how the utils image for the [smoke tests](ci-cd.md#tests) is made.

| Branch | Tags pushed | How the tag is chosen |
|---|---|---|
| `main` | `<YYYYMMDD>-<short-sha>` and `:latest` | The HEAD commit's date and 7-character SHA, e.g. `20261022-72ee249`, so every tag maps to exactly one commit |
| any other branch | `:<branch>` only | The branch name itself, so parallel branches never overwrite each other's images or `:latest` |

A `main` build refuses to push from a dirty tree, because its tag names the HEAD commit. A branch build also passes `--from <registry>/<base>:<branch>` when the Dockerfile's first `FROM` is a `:latest` image from this registry and that base has a build on the same branch, so a branch that changes what the base image bakes in is testable end to end. `merge_branch.sh` deletes the branch's image tags once the branch merges, and the merged change reaches `:latest` only when the image is rebuilt from the main checkout. The earlier `kj_V<N>` auto-increment tags are no longer produced; the existing ones remain in Artifact Registry as history.

Every upstream tool is pinned in `versions.env`: released tools by version, and tools built from a git clone by commit, as `<image>__<NAME>_COMMIT` keys - for example `svan__SVAN_COMMIT` and `stranalysis__STR_ANALYSIS_COMMIT`. The commits are the ones baked into the current `:latest` images, so a rebuild reproduces them rather than picking up the upstream default branch. `Dockerfile.palmer` is intentionally left unpinned. [`tests/unit/test_repo_invariants.py`](../tests/unit/test_repo_invariants.py) fails on an `ARG` with no key, a key no Dockerfile declares, or any other unpinned `git clone`. [`.hadolint.yaml`](../.hadolint.yaml) configures `hadolint` for the Dockerfiles, which CI runs - see [CI/CD](ci-cd.md#github-actions).

Two consequences worth knowing:
- Because the build context is the repository root, `Dockerfile.utils`, `Dockerfile.hail` and `Dockerfile.vephail` bake the live `scripts/` tree into the image at `/opt/scripts`, so a change under `scripts/` reaches a running task only once those images and the images inheriting from `utils` are rebuilt - see [Scripts](repository-structure.md#scripts) for how each script reaches a task.
- The script needs an authenticated `gcloud` for the branch base-image lookup and a `podman` logged in to the registry for the push; neither is needed with `--no-push` on a branch.

**Which tag to use:** WDL tasks never hardcode a docker image URI - the `String docker` task input is always supplied by the caller via Terra workspace data, per [Conventions](conventions.md). Workspace data should point at the `:latest` tag for each image, since every `main` build retags `:latest` to the newest build. The dated tags exist purely as an immutable version history, for pinning or rolling back to a specific prior build. While testing a branch, point that branch's Terra method config at the `:<branch>` tag instead, and drop back to `:latest` after the merge and rebuild. Workspace data `hail_docker` should point at `kj-dockers/hail:latest` once the first `main` build of the `hail` image has been pushed.


## Repository
Built from [dockerfiles/](../dockerfiles/) in this repository.

| Argument Name | Image | Dockerfile |
|---|---|---|
| `intact_mei_docker` | kj-dockers/intactmei:latest | `Dockerfile.intactmei` |
| `palmer_docker` | kj-dockers/palmer:latest | `Dockerfile.palmer` |
| `cpg_docker` | kj-dockers/cpgtools:latest | `Dockerfile.cpgtools` |
| `vrs_docker` | kj-dockers/vrs:latest | `Dockerfile.vrs` |
| `l1meaid_docker` | kj-dockers/l1meaid:latest | `Dockerfile.l1meaid` |
| `vep_hail_docker` | kj-dockers/vephail:latest | `Dockerfile.vephail` |
| `hail_docker` | kj-dockers/hail:latest | `Dockerfile.hail` |
| `stranalysis_docker` | kj-dockers/stranalysis:latest | `Dockerfile.stranalysis` |
| `utils_docker` | kj-dockers/utils:latest | `Dockerfile.utils` |
| `pav_docker` | kj-dockers/pav:latest | `Dockerfile.pav` |
| `vamos_docker` | kj-dockers/vamos:latest | `Dockerfile.vamos` |
| `trgt_lps_docker` | kj-dockers/trgt-lps:latest | `Dockerfile.trgtlps` |
| `trgt_docker` | kj-dockers/trgt:latest | `Dockerfile.trgt` |
| `mosdepth_docker` | kj-dockers/mosdepth:latest | `Dockerfile.mosdepth` |
| `mosdepthstream_docker` | kj-dockers/mosdepthstream:latest | `Dockerfile.mosdepthstream` |
| `kanpig_docker` | kj-dockers/kanpig:latest | `Dockerfile.kanpig` |
| `svan_docker` | kj-dockers/svan:latest | `Dockerfile.svan` |
| `hificnv_docker` | kj-dockers/hificnv:latest | `Dockerfile.hificnv` |
| `sawfish_docker` | kj-dockers/sawfish:latest | `Dockerfile.sawfish` |
| `minimap2_docker` | kj-dockers/minimap2:latest | `Dockerfile.minimap2` |
| `hifiasm_docker` | kj-dockers/hifiasm:latest | `Dockerfile.hifiasm` |
| `pbsv_docker` | kj-dockers/pbsv:latest | `Dockerfile.pbsv` |
| `sniffles_docker` | kj-dockers/sniffles:latest | `Dockerfile.sniffles` |

`Dockerfile.utils` is the base image for most other repo Dockerfiles. `Dockerfile.hail` is `hailgenetics/hail:0.2.122` plus the `scripts/` tree, and `Dockerfile.vephail` builds on the VEP image with Hail pinned through `HAIL_VERSION`; neither inherits from `utils`. Every Dockerfile in `dockerfiles/` maps to exactly one argument above; retired ones live in [`archive/dockerfiles/`](../archive/dockerfiles/) and are not built.


## Collaborators
Built/maintained by collaborators.

| Argument Name | Image | Source |
|---|---|---|
| `sv_base_mini_docker` | gatk-sv/sv-base-mini:2024-10-25-... | GATK-SV |
| `sv_pipeline_docker` | gatk-sv/sv-pipeline:2025-10-02-... | GATK-SV |
| `gatk_docker` | gatk-sv/gatk:mw-gatk-sv-672d85 | GATK-SV |
| `gatk_sv_lr_docker` | kj-dockers/sv-pipeline:kj_V65 | `kj_project_gnomad_lr` branch of GATK-SV |
| `hiphase_docker` | broad-dsp-lrma/hangsuunc/hiphase:v1.5.0 | Hang Su |
| `hiphase_preprocess_docker` | hangsuunc/cleanvcf:v1 | Hang Su |
| `remap_docker` | quay.io/ymostovoy/lr-remap | Yulia Mostovoy |
| `sv_integration_docker` | broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages | Fabio Cunial |
| `xgb_scoring_docker` | broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb | GATK snapshot |
| `minimap_docker` | eichlerlab/assembly_eval:0.2 | Eichler Lab |
| `automop_docker` | broad-dsde-methods/automop:0.1 | Broad DSP |


## Published
Public images.

| Argument Name | Image | Source |
|---|---|---|
| `repeatmasker_docker` | dfam/tetools:1.8 | Dfam consortium |
| `glnexus_docker` | ghcr.io/dnanexus-rnd/glnexus:v1.4.1 | DNAnexus |
