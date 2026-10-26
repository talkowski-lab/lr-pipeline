# Docker images
This document lists the Docker images used by the pipeline - repo-built, collaborator-built, and published third-party - and describes how the repository's own images are built, tagged and pushed.


## Building and Pushing
Every container image is built from a `dockerfiles/Dockerfile.<image-name>` file, where `<image-name>` is all-lowercase and is exactly the name the image is pushed under - e.g. `Dockerfile.utils` builds `utils`, `Dockerfile.stranalysis` builds `stranalysis`. There is no separate name-mapping file; the Dockerfile suffix mechanically is the image name, so `Dockerfile.trgtlps` is pushed as `trgtlps`. An older `trgt-lps` image in the registry predates the build script and is not rebuilt.

Images are pushed to Artifact Registry at:
```
us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/<image-name>
```

This registry is private. To run the pipeline elsewhere, build the images you need from [`dockerfiles/`](../dockerfiles/), push them to a registry your Terra workspace can pull from, and point the workspace data at them. Set `DOCKER_REGISTRY` to that registry - e.g. `DOCKER_REGISTRY=us-central1-docker.pkg.dev/<project>/<repository> dockerfiles/build_docker.sh utils` - and the script pushes there instead. Build `utils` first, then `trgt`, since the other repository images start from one of them. Each of those Dockerfiles declares `ARG BASE_IMAGE=<registry>/<base>:latest` before `FROM ${BASE_IMAGE}`, so point it at your registry with `--build-arg BASE_IMAGE=<your-registry>/utils:latest` instead of editing the Dockerfile.

Build and push an image with the helper script, which takes the image name as its only argument, plus an optional `--no-push` in any position:
```bash
dockerfiles/build_docker.sh [--no-push] <image-name>
# e.g. dockerfiles/build_docker.sh utils
```

The script does the following:
- Resolves the build arguments through [`dockerfiles/build_args.sh`](../dockerfiles/build_args.sh) `<image-name>`, which prints a `NAME=value` line for every `ARG` declared in the Dockerfile without an inline default, taken from [`dockerfiles/versions.env`](../dockerfiles/versions.env) where it is keyed as `<image-name>__<ARG_NAME>`. CI uses the same script, so local and CI builds resolve identically. A missing key fails the build.
- Builds with `podman build --platform linux/amd64` from the repository root - meaning the worktree it is run from, so each branch builds its own tree.
- Labels the build with the OCI keys `org.opencontainers.image.source`, `.revision` (the full HEAD commit), `.version` (the tag it is built as) and `.created`, so `podman inspect` on any pulled image names the commit behind it.
- Tags the build by the branch it was built on, since branches are developed in parallel worktrees - see [Branches](conventions.md#branches).
- With `--no-push`, builds and tags locally only, as `localhost/<image-name>:<tag>`, and pushes nothing. A push ends by printing the pushed manifest digest, `<registry>/<image-name>@sha256:...`. This is how contributors without registry access build images, and how the utils image for the [smoke tests](ci-cd.md#tests) is made.

| Branch | Tags pushed | How the tag is chosen |
|---|---|---|
| `main` | `<YYYYMMDD>-<short-sha>` and `:latest` | The HEAD commit's date and 7-character SHA, e.g. `20261022-72ee249`, so every tag maps to exactly one commit |
| any other branch | `:<branch>` only | The branch name itself, so parallel branches never overwrite each other's images or `:latest` |

A `main` build refuses to push from a dirty tree, because its tag names the HEAD commit. A branch build also passes `--build-arg BASE_IMAGE=<registry>/<base>:<branch>` when the Dockerfile's `ARG BASE_IMAGE` default is a `:latest` image from this registry and that base has a build on the same branch, so a branch that changes what the base image bakes in is testable end to end; otherwise the `:latest` default is used. `merge_branch.sh` deletes the branch's image tags once the branch merges, and the merged change reaches `:latest` only when the image is rebuilt from the main checkout.

Every upstream tool is pinned in `versions.env`: released tools by version, and tools built from a git clone by commit, as `<image>__<NAME>_COMMIT` keys - for example `svan__SVAN_COMMIT` and `stranalysis__STR_ANALYSIS_COMMIT`. The commits are the ones baked into the current `:latest` images, so a rebuild reproduces them rather than picking up the upstream default branch. Language packages, toolchains and SDKs are pinned there too, and R packages come from a dated CRAN snapshot, so a rebuild reproduces the same package versions. The only pins written inline are the apt package versions, since they follow the Ubuntu release the base image uses. [`tests/unit/test_repo_invariants.py`](../tests/unit/test_repo_invariants.py) fails on an `ARG` with no key, a key no Dockerfile declares, any other unpinned `git clone`, an unpinned `pip install`, `apt-get upgrade` or `|| true`, or a Dockerfile that names the private registry in `FROM` instead of `${BASE_IMAGE}`. [`.hadolint.yaml`](../.hadolint.yaml) configures `hadolint` for the Dockerfiles, which CI runs - see [CI/CD](ci-cd.md#github-actions).

[`.github/workflows/docker-build.yml`](../.github/workflows/docker-build.yml) builds the changed Dockerfiles in CI, without pushing, on top of a `utils:ci` image it builds itself with `BASE_IMAGE=utils:ci` - see [CI/CD](ci-cd.md#docker-images).

Two consequences worth knowing:
- Because the build context is the repository root, `Dockerfile.utils`, `Dockerfile.hail` and `Dockerfile.vephail` bake the live `scripts/` tree into the image at `/opt/scripts`, so a change under `scripts/` reaches a running task only once those images and the images inheriting from `utils` are rebuilt - see [Scripts](repository-structure.md#scripts) for how each script reaches a task.
- The script needs an authenticated `gcloud` for the branch base-image lookup and a `podman` logged in to the registry for the push; neither is needed with `--no-push` on a branch.

**Which tag to use:** WDL tasks never hardcode a docker image URI - the `String docker` task input is always supplied by the caller via Terra workspace data, per [Conventions](conventions.md). For development, workspace data points at the `:latest` tag of each image, since every `main` build retags `:latest` to the newest build. A production or publication run should instead pin each image by its dated tag, or by the digest the build printed, so the exact images behind a callset stay recoverable after `:latest` moves; the dated tags are also how to roll back to a prior build. [`dockerfiles/audit_images.sh`](../dockerfiles/audit_images.sh) lists the digest and dated tag behind every image's `:latest`, and exits non-zero naming any image whose `:latest` has no dated tag. While testing a branch, point that branch's Terra method config at the `:<branch>` tag instead, and drop back to `:latest` after the merge and rebuild.


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
| `trgt_lps_docker` | `trgtlps` | [`Dockerfile.trgtlps`](../dockerfiles/Dockerfile.trgtlps) |
| `utils_docker` | `utils` | [`Dockerfile.utils`](../dockerfiles/Dockerfile.utils) |
| `vamos_docker` | `vamos` | [`Dockerfile.vamos`](../dockerfiles/Dockerfile.vamos) |
| `vep_hail_docker` | `vephail` | [`Dockerfile.vephail`](../dockerfiles/Dockerfile.vephail) |
| `vrs_docker` | `vrs` | [`Dockerfile.vrs`](../dockerfiles/Dockerfile.vrs) |

Every Dockerfile in `dockerfiles/` maps to exactly one argument above; retired ones live in [`archive/dockerfiles/`](../archive/dockerfiles/) and are not built.


## Collaborators
Built and maintained outside this repository. The Digest column records the manifest digest of every image that is untagged or whose tag is moved, so the exact image behind a callset can be recovered; workspace data should pin such images by digest (`<image>@sha256:...`) rather than by tag. Digests were read from the registries on 2026-10-23.

| Argument Name | Image | Source | Digest |
|---|---|---|---|
| `automop_docker` | `us.gcr.io/broad-dsde-methods/automop:0.1` | Broad Data Sciences Platform | |
| `gatk_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/gatk:mw-gatk-sv-672d85` | [GATK-SV](https://github.com/broadinstitute/gatk-sv) development build, kept over a GATK release because it carries `TrainSVGenotyping` and `GenotypeSVs`, which no release includes | |
| `hiphase_docker` | `us.gcr.io/broad-dsp-lrma/hangsuunc/hiphase:v1.5.0` | Hang Su | `sha256:a7cca48696f6134ea8ba32ba5ddb0bbc0e55013357071f915ca937993c23a23c` |
| `hiphase_preprocess_docker` | `hangsuunc/cleanvcf:v1` | Hang Su | `sha256:73fd5d2f982f95e97e67286783646f3c5c9d5c85cf2b43e8d2cd717bbf8f98ae` |
| `remap_docker` | `quay.io/ymostovoy/lr-remap:latest` | Yulia Mostovoy; a modified Truvari whose `anno remap` takes `--mm2-threshold` and `--cov-threshold`, which no Truvari release has | `sha256:3410df5da3631e9eab42e6e7cb750a0a91d1e9b94a523bd632586a1c565a6e6d` |
| `sv_integration_docker` | `us.gcr.io/broad-dsp-lrma/fcunial/callset_integration_phase2_workpackages` | Fabio Cunial; bundles the compiled Java and scripts of [callset_integration_phase2](https://github.com/fabio-cunial/callset_integration_phase2) | untagged, `:latest` moves; `sha256:ad7202bbeb8eff86df592e334026ff262ec66b0057c3f7e5fa256e7761f1e9e9` |
| `xgb_scoring_docker` | `us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk:sl_aou_lr_intrasample_filtering_xgb` | GATK snapshot from the `sl_aou_lr_intrasample_filtering_xgb` branch, running the XGBoost scoring of `MergeSampleSVCallsets` | `sha256:a35a80cf67bd366fdc6a34a5b99e2632a607cebd7f007bb607de2ce3f1de8e8b` |


## Published
Public images.

| Argument Name | Image | Source |
|---|---|---|
| `deepvariant_docker` | `google/deepvariant:1.6.0` | [DeepVariant](https://github.com/google/deepvariant) |
| `glnexus_docker` | `ghcr.io/dnanexus-rnd/glnexus:v1.4.1` | DNAnexus |
| `repeatmasker_docker` | `dfam/tetools:1.8` | Dfam consortium |
| `sv_pipeline_docker` | `us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-10-02-v1.1-483973d6` | [GATK-SV](https://github.com/broadinstitute/gatk-sv) |
