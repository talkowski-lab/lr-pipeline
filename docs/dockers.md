# Docker images
This document lists the Docker images used by the pipeline - repo-built, collaborator-built, and published third-party.


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
