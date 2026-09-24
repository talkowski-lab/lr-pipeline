# Docker images
This document lists the Docker images used by the pipeline - repo-built, collaborator-built, and published third-party - and describes how the repository's own images are built, tagged and pushed.


## Building and Pushing
Every container image is built from a `dockerfiles/Dockerfile.<image-name>` file, where `<image-name>` is all-lowercase and is exactly the name the image is pushed under - e.g. `Dockerfile.utils` builds `utils`, `Dockerfile.stranalysis` builds `stranalysis`. There is no separate name-mapping file; the Dockerfile suffix mechanically is the image name. The one exception is `Dockerfile.trgtlps`, which is pushed as `trgt-lps` with a hyphen, so do not assume the two always match.

Images are pushed to Artifact Registry at:
```
us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/<image-name>
```

Build and push an image with the helper script, which takes the image name as its only argument:
```bash
dockerfiles/build_docker.sh <image-name>
# e.g. dockerfiles/build_docker.sh utils
```

The script does the following:
- Resolves any `ARG` declared in the Dockerfile without an inline default from [`dockerfiles/versions.env`](../dockerfiles/versions.env), keyed as `<image-name>__<ARG_NAME>`, and passes each as a `--build-arg`. That file is the single source of truth for pinned tool and library versions.
- Builds with `podman build --platform linux/amd64` from the repository root - meaning the worktree it is run from, so each branch builds its own tree.
- Tags the build by the branch it was built on, since branches are developed in parallel worktrees - see [Branches](conventions.md#branches).

| Branch | Tags pushed | How the tag is chosen |
|---|---|---|
| `main` | `kj_V<N>` and `:latest` | `gcloud artifacts docker tags list` gives the highest `kj_V<N>` already pushed, which the script increments |
| `kj-<topic>` | `:<branch>` only | The branch name itself, so parallel branches never overwrite each other's images or `:latest` |

A branch build also passes `--from <registry>/<base>:<branch>` when the Dockerfile's first `FROM` is a `:latest` image from this registry and that base has a build on the same branch, so a branch that changes what the base image bakes in is testable end to end. `merge_branch.sh` deletes the branch's image tags once the branch merges, and the merged change reaches `:latest` only when the image is rebuilt from the main checkout. Two `main` builds of the same image started at once would still race for the same `kj_V<N>`; build one at a time.

Two consequences worth knowing:
- Because the build context is the repository root, `Dockerfile.utils` bakes the live `scripts/` tree into the image via `COPY ./scripts /opt/scripts`. Any change under `scripts/` therefore requires rebuilding `utils`, and then any image that inherits from it, before the change reaches a running task. The exception is the Hail scripts, which workflows fetch by URL at run time - see [Scripts](repository-structure.md#scripts).
- The script needs an authenticated `gcloud` for tag discovery and a `podman` logged in to the registry for the push.

**Which tag to use:** WDL tasks never hardcode a docker image URI - the `String docker` task input is always supplied by the caller via Terra workspace data, per [Conventions](conventions.md). Workspace data should point at the `:latest` tag for each image, since every `main` build retags `:latest` to the newest build. The `kj_V<N>` tags exist purely as an immutable version history, for pinning or rolling back to a specific prior build. While testing a branch, point that branch's Terra method config at the `:<branch>` tag instead, and drop back to `:latest` after the merge and rebuild.


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

`Dockerfile.utils` is the base image for most other repo Dockerfiles. Every Dockerfile in `dockerfiles/` maps to exactly one argument above; retired ones live in [`archive/dockerfiles/`](../archive/dockerfiles/) and are not built.


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
| `minimap_docker` | eichlerlab/assembly_eval:0.2 | Eichler Lab |
| `automop_docker` | broad-dsde-methods/automop:0.1 | Broad DSP |


## Published
Public images.

| Argument Name | Image | Source |
|---|---|---|
| `repeatmasker_docker` | dfam/tetools:1.8 | Dfam consortium |
| `hail_docker` | hailgenetics/hail:0.2.122 | Hail team |
| `glnexus_docker` | ghcr.io/dnanexus-rnd/glnexus:v1.4.1 | DNAnexus |
