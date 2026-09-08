# Docker images

Every `String *_docker` input declared across `wdl/` (source of truth — checked
against every workflow, not just the Terra workspace snapshot), mapped to
where it's built. Grouped by who builds/maintains each image.

## Repo dockers

Built from `Dockerfiles` in this repo ([dockerfiles/](../dockerfiles/)), pushed to
`us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers/`.

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
| `whatshap_docker` | kj-dockers/whatshap:latest | `Dockerfile.whatshap` |
| `hificnv_docker` | kj-dockers/hificnv:latest | `Dockerfile.hificnv` |
| `sawfish_docker` | kj-dockers/sawfish:latest | `Dockerfile.sawfish` |

`Dockerfile.utils` is the base image for most other repo Dockerfiles. Every
Dockerfile currently in `dockerfiles/` maps to exactly one argument above —
none are unused.

Note: `Dockerfile.trgtlps` is pushed under the image name `trgt-lps` (hyphenated),
not `trgtlps` as `build_docker.sh`'s naming convention would derive from the
filename — a pre-existing drift from [repository-structure.md](repository-structure.md)'s
"Dockerfile suffix mechanically is the image name" rule. Works today since the
attribute matches the actual pushed name; just don't assume the two always match.

## Collaborator dockers

Built/maintained by collaborators, outside this repo.

| Argument Name | Image | Source |
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
| `vcfdist_docker` | timd1/vcfdist:v2.6.4 | [TimD1/vcfdist](https://github.com/TimD1/vcfdist) author's own published image. Only consumer, `VcfDist`/`VcfDistCohort`, is archived — see [archive/docs/workflows.md](../archive/docs/workflows.md#tools) |

## Published

Public images, not custom-built for this project.

| Argument Name | Image | Source |
|---|---|---|
| `repeatmasker_docker` | dfam/tetools:1.8 | Dfam consortium |
| `hail_docker` | hailgenetics/hail:0.2.105 | Hail team |
| `glnexus_docker` | ghcr.io/dnanexus-rnd/glnexus:v1.4.1 | DNAnexus |
