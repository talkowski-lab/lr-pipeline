# Docker images
This document lists the Docker images used by the pipeline - repo-built, collaborator-built, and published third-party - and describes how the repository's own images are built, tagged and pushed.


## Building and Pushing
Every container image is built from a `dockerfiles/Dockerfile.<image-name>` file, where `<image-name>` is all-lowercase and is exactly the name the image is pushed under - e.g. `Dockerfile.utils` builds `utils`, `Dockerfile.stranalysis` builds `stranalysis`. There is no separate name-mapping file; the Dockerfile suffix mechanically is the image name. The one exception is `Dockerfile.trgtlps`, which is pushed as `trgt-lps` with a hyphen, so do not assume the two always match.

Images are pushed to Artifact Registry at:
```
us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/<image-name>
```

This registry is private. To run the pipeline elsewhere, build the images you need from [`dockerfiles/`](../dockerfiles/), push them to a registry your Terra workspace can pull from, and point the workspace data at them. Set `DOCKER_REGISTRY` to that registry - e.g. `DOCKER_REGISTRY=us-central1-docker.pkg.dev/<project>/<repository> dockerfiles/build_docker.sh utils` - and the script pushes there instead. Build `utils` first, then `trgt`, since the other repository images start `FROM` one of them; their `FROM` lines name this registry, so replace it with yours before building.

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
Built from [dockerfiles/](../dockerfiles/) in this repository and pushed to the registry above, so each image below is `us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/<image>:latest`.

| Argument Name | Image | Dockerfile |
|---|---|---|
| `cpg_docker` | `cpgtools` | [`Dockerfile.cpgtools`](../dockerfiles/Dockerfile.cpgtools) |
| `hail_docker` | `hail` | [`Dockerfile.hail`](../dockerfiles/Dockerfile.hail) |
| `hifiasm_docker` | `hifiasm` | [`Dockerfile.hifiasm`](../dockerfiles/Dockerfile.hifiasm) |
| `hificnv_docker` | `hificnv` | [`Dockerfile.hificnv`](../dockerfiles/Dockerfile.hificnv) |
| `intact_mei_docker` | `intactmei` | [`Dockerfile.intactmei`](../dockerfiles/Dockerfile.intactmei) |
| `kanpig_docker` | `kanpig` | [`Dockerfile.kanpig`](../dockerfiles/Dockerfile.kanpig) |
| `l1meaid_docker` | `l1meaid` | [`Dockerfile.l1meaid`](../dockerfiles/Dockerfile.l1meaid) |
| `minimap2_docker` | `minimap2` | [`Dockerfile.minimap2`](../dockerfiles/Dockerfile.minimap2) |
| `mosdepth_docker` | `mosdepth` | [`Dockerfile.mosdepth`](../dockerfiles/Dockerfile.mosdepth) |
| `mosdepthstream_docker` | `mosdepthstream` | [`Dockerfile.mosdepthstream`](../dockerfiles/Dockerfile.mosdepthstream) |
| `palmer_docker` | `palmer` | [`Dockerfile.palmer`](../dockerfiles/Dockerfile.palmer) |
| `pav_docker` | `pav` | [`Dockerfile.pav`](../dockerfiles/Dockerfile.pav) |
| `pbsv_docker` | `pbsv` | [`Dockerfile.pbsv`](../dockerfiles/Dockerfile.pbsv) |
| `sawfish_docker` | `sawfish` | [`Dockerfile.sawfish`](../dockerfiles/Dockerfile.sawfish) |
| `sniffles_docker` | `sniffles` | [`Dockerfile.sniffles`](../dockerfiles/Dockerfile.sniffles) |
| `stranalysis_docker` | `stranalysis` | [`Dockerfile.stranalysis`](../dockerfiles/Dockerfile.stranalysis) |
| `svan_docker` | `svan` | [`Dockerfile.svan`](../dockerfiles/Dockerfile.svan) |
| `trgt_docker` | `trgt` | [`Dockerfile.trgt`](../dockerfiles/Dockerfile.trgt) |
| `trgt_lps_docker` | `trgt-lps` | [`Dockerfile.trgtlps`](../dockerfiles/Dockerfile.trgtlps) |
| `utils_docker` | `utils` | [`Dockerfile.utils`](../dockerfiles/Dockerfile.utils) |
| `vamos_docker` | `vamos` | [`Dockerfile.vamos`](../dockerfiles/Dockerfile.vamos) |
| `vep_hail_docker` | `vephail` | [`Dockerfile.vephail`](../dockerfiles/Dockerfile.vephail) |
| `vrs_docker` | `vrs` | [`Dockerfile.vrs`](../dockerfiles/Dockerfile.vrs) |

`Dockerfile.utils` is the base image for most other repo Dockerfiles. `Dockerfile.hail` is `hailgenetics/hail:0.2.122` plus the `scripts/` tree, and `Dockerfile.vephail` builds on the VEP image with Hail pinned through `HAIL_VERSION`; neither inherits from `utils`. Every Dockerfile in `dockerfiles/` maps to exactly one argument above; retired ones live in [`archive/dockerfiles/`](../archive/dockerfiles/) and are not built.


## Collaborators
Built and maintained outside this repository.

| Argument Name | Image | Source |
|---|---|---|
| `automop_docker` | `us.gcr.io/broad-dsde-methods/automop:0.1` | Broad Data Sciences Platform |
| `deepvariant_docker` | `us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.6.0` | [long-read-pipelines](https://github.com/broadinstitute/long-read-pipelines), the source of the `DeepVariant` workflow |
| `deepvariant_gpu_docker` | `us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.6.0-gpu` | [long-read-pipelines](https://github.com/broadinstitute/long-read-pipelines) |
| `gatk_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-672d85` | [GATK-SV](https://github.com/broadinstitute/gatk-sv) |
| `gatk_sv_lr_docker` | `us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/sv-pipeline:latest` | `sv-pipeline` image built from the `kj_project_gnomad_lr` branch of [GATK-SV](https://github.com/broadinstitute/gatk-sv/tree/kj_project_gnomad_lr) |
| `hiphase_docker` | `us.gcr.io/broad-dsp-lrma/hangsuunc/hiphase:v1.5.0` | Hang Su |
| `hiphase_preprocess_docker` | `hangsuunc/cleanvcf:v1` | Hang Su |
| `minimap_docker` | `eichlerlab/assembly_eval:0.2` | Eichler Lab |
| `remap_docker` | `quay.io/ymostovoy/lr-remap:latest` | Yulia Mostovoy |
| `resource_visualization_docker` | `us.gcr.io/broad-dsp-lrma/lr-resource-visual:0.1.1` | [long-read-pipelines](https://github.com/broadinstitute/long-read-pipelines) |
| `sv_base_mini_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52` | [GATK-SV](https://github.com/broadinstitute/gatk-sv) |
| `sv_integration_docker` | `broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages` | Fabio Cunial |
| `sv_pipeline_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-10-02-v1.1-483973d6` | [GATK-SV](https://github.com/broadinstitute/gatk-sv) |
| `xgb_scoring_docker` | `broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb` | GATK snapshot |


## Published
Public images.

| Argument Name | Image | Source |
|---|---|---|
| `glnexus_docker` | `ghcr.io/dnanexus-rnd/glnexus:v1.4.1` | DNAnexus |
| `repeatmasker_docker` | `dfam/tetools:1.8` | Dfam consortium |
