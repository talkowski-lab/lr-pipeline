# Docker images

Docker images used across lr-pipeline workflows, as configured in the Terra
workspace (`attributes.tsv`). Grouped by who builds/maintains each image.

## Repo dockers

Built from `Dockerfiles` in this repo ([dockerfiles/](../dockerfiles/)), pushed to
`us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/`.

| Attribute | Image | Dockerfile |
|---|---|---|
| `intact_mei_docker` | kj-dockers/intactmei:latest | `Dockerfile.intactmei` |
| `palmer_docker` | kj-dockers/palmer:latest | `Dockerfile.palmer` |
| `cpg_docker` | kj-dockers/cpgtools:latest | `Dockerfile.cpgtools` |
| `vrs_docker` | kj-dockers/vrs:latest | `Dockerfile.vrs` |
| `l1meaid_docker` | kj-dockers/l1meaid:latest ⚠️ | `Dockerfile.l1meaid` |
| `vep_hail_docker` | kj-dockers/vephail:latest | `Dockerfile.vephail` |
| `stranalysis_docker` | kj-dockers/stranalysis:latest | `Dockerfile.stranalysis` |
| `utils_docker` | kj-dockers/utils:latest | `Dockerfile.utils` |
| `pav_docker` | kj-dockers/pav:latest | `Dockerfile.pav` |
| `vamos_docker` | kj-dockers/vamos:latest | `Dockerfile.vamos` |
| `trgt_lps_docker` | kj-dockers/trgt-lps:latest | `Dockerfile.trgtlps` |
| `trgt_docker` | kj-dockers/trgt:latest | `Dockerfile.trgt` |
| `mosdepth_docker` | kj-dockers/mosdepth:latest | `Dockerfile.mosdepth` |
| `kanpig_docker` | kj-dockers/kanpig:latest | `Dockerfile.kanpig` |
| `svan_docker` | kj-dockers/svan:latest | `Dockerfile.svan` |
| `whatshap_docker` | kj-dockers/whatshap:latest | `Dockerfile.whatshap` |

`Dockerfile.utils` is the base image for most other repo Dockerfiles.

## Collaborator dockers

Built/maintained by collaborators, outside this repo.

| Attribute | Image | Source |
|---|---|---|
| `sv_base_mini_docker` | gatk-sv/sv-base-mini:2024-10-25-... | gatk-sv, see [dockers.json](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers.json) |
| `sv_pipeline_docker` | gatk-sv/sv-pipeline:2025-10-02-... | gatk-sv, dockers.json |
| `gatk_docker` | gatk-sv/gatk:mw-gatk-sv-672d85 | gatk-sv fork build (tag `mw-...`) |
| `gatk_sv_lr_docker` | kj-dockers/sv-pipeline:kj_V65 | built off `kj_project_gnomad_lr` branch of gatk-sv (not upstream dockers.json — see [[project_cleanvcf_migration]]) |
| `hiphase_docker` | broad-dsp-lrma/hangsuunc/hiphase:v1.5.0 | Hangsu |
| `hiphase_preprocess_docker` | hangsuunc/cleanvcf:v1 | Hangsu (CleanVcf) |
| `remap_docker` | quay.io/ymostovoy/lr-remap | Yulia Mostovoy |
| `minimap_docker` | eichlerlab/assembly_eval:0.2 | Eichler lab |
| `automop_docker` | broad-dsde-methods/automop:0.1 | Broad DSP team (same group as `hiphase_docker`/`hiphase_preprocess_docker`) |

## Published

Public images, not custom-built for this project.

| Attribute | Image | Source |
|---|---|---|
| `repeatmasker_docker` | dfam/tetools:1.8 | Dfam consortium |
| `hail_docker` | hailgenetics/hail:0.2.105 | Hail team |
| `glnexus_docker` | ghcr.io/dnanexus-rnd/glnexus:v1.4.1 | DNAnexus |

`Dockerfile.sawfish`, `Dockerfile.hificnv`, and `Dockerfile.mosdepthstream` also
exist in [dockerfiles/](../dockerfiles/) with no matching workspace attribute
currently wired up.

⚠️ `l1meaid_docker` points at `kj-dockers/l1meaid:latest`, but that tag has
never been pushed — the registry only has `l1meaid:kj_V1`. Any task using this
attribute will fail to pull. Fix by running `dockerfiles/build_docker.sh l1meaid`
(pushes a new version and retags `:latest`), or point the attribute at
`l1meaid:kj_V1` directly until then.
