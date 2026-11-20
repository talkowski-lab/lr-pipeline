# Instructions

## Project Overview
This repository contains all scripts, workflows, and processes for annotating long-read variant callsets. It processes structural variants (SVs), mobile element insertions (MEIs), tandem repeats (TRs), and other complex variants from cohorts including HPRC (232 samples), HGSVC (65 samples), and All of Us Phase 1 (1027 samples).

**Stack:** WDL 1.0 (executed via Cromwell on Terra), Python 3.8+, Bash, Hail. Cloud platform is GCP (GCS URIs `gs://`). Reference genome is GRCh38.


## Repository Layout
```
wdl/
  00_annotate_SV/                 # Annotations added to long-read variants for the gnomAD-LR resource (prefix: Annotate*)
  01_insertion_source_analyses/   # Origin of inserted sequences (AnnotateInsertionOrigins, AnnotateSmallInsertionOrigins)
  02_variant_TR_constraint/       # Per-sample burden, per-gene variant / TR constraint, CADD scores
  03_expression_alterations/      # Expression alterations (Step01_CountPhasedHetPerGene: phased hets per gene for Iso-Seq phasing)
  04_methylation_analyses/        # Methylation analyses (placeholder)
  05_CNV_depth_calling/           # Depth CNV calling, genomic-disorder breakpoints, read recurrence, local assembly
  06_compound_het/                # Compound-het / double-hit analyses (placeholder)
  07_ancestry_matching/           # PCA ancestry projection and ancestry-matched control selection
  annotation_utils/  # VCF manipulation and utility workflows
  tools/             # Individual bioinformatics tool wrappers
  utils/             # Shared structs (Structs.wdl) and helper tasks (Helpers.wdl)
scripts/
  vep/               # VEP and Hail-based annotation scripts
  mei/               # MEI analysis scripts
  benchmark/         # Benchmarking scripts
  miscellaneous/     # Other utility scripts
dockerfiles/         # Dockerfile.<ToolName> for each container
data/
  references/        # Reference genomes, catalogs, BED files
  base_vcfs/         # Test/base VCF files
  metadata/          # Sample metadata, pedigrees, ancestry
  qc_annotations/    # QC annotation outputs
archive/             # Deprecated workflows, scripts, and dockerfiles (do not modify)
```

## Build, Validate & Lint
### WDL Validation
CI validates all WDL files on push/PR to `main` when `wdl/**` changes:
```bash
wget -O womtool.jar https://github.com/broadinstitute/cromwell/releases/download/87/womtool-87.jar
find wdl -type f -name "*.wdl" -exec java -jar womtool.jar validate {} \;
```
Always run `womtool validate` on any new or modified WDL before committing. You can find a local copy of the `womtool-87.jar` at `/Users/kjaising/Desktop/Work/Miscellaneous/Software/womtool-87.jar` - whenever possible, opt to use this rather than downloading the remote version.

### Python Linting
CI lints `scripts/**` on push/PR to `main`:
```bash
flake8 scripts/ --max-line-length=130
```
The `.flake8` config extends ignore for `E203` and `W503`, and excludes `data/*` and `archive/*`.

### Dockstore
All workflows that are directly run in the pipeline must have entries in `.dockstore.yml`. Annotation workflows go under the `# Annotation Workflows` section, utilities under `# Annotation Utilities`, and tools under `# Tool Workflows`. Each entry follows this format:
```yaml
- subclass: WDL
  name: WorkflowName
  primaryDescriptorPath: /wdl/<subdirectory>/WorkflowName.wdl
  filters:
    branches:
      - main
    tags:
      - /.*/
```


## WDL Conventions
### Version & Imports
- Always use `version 1.0`.
- Import `Structs.wdl` and `Helpers.wdl` from `../utils/` using relative paths.
- Never rename imports with the `as` operator.

### Naming
- **Workflows & tasks**: PascalCase. Proper noun acronyms (e.g., PALMER, TRGT, L1MEAID, SVAN) retain their capitalization.
- **Variables, inputs & outputs**: snake_case. Capitalized nouns are exceptions (e.g., `PALMER_vcf`).
- **Workflow file names** must match the workflow name defined within them.
- Annotation workflows in `wdl/00_annotate_SV/` must be prefixed with `Annotate`.

### File-Type Suffixes
- Use `fa` not `fasta` — e.g., `ref_fa`.
- Use `fai` not `fasta_index`, `fasta_fai`, or `fa_fai` — e.g., `ref_fai`.
- Use `vcf_idx` not `vcf_index` or `vcf_tbi` — e.g., `input_vcf_idx`.
- VCF variables use the `_vcf` suffix and are always paired with a `_vcf_idx` index.

### Workflow Structure
Workflows are ordered as follows, with each section separated by a single blank line:
1. Imports.
2. Inputs.
3. Dynamically generated variables.
4. Task calls.
5. Outputs.

### Workflow Inputs
Inputs are grouped in this order, with blank lines between groups:
1. **Core files** — VCFs, BAMs, indexes, contigs, prefix.
2. **Parameters** — prefixes, modes, tool arguments, PED files, metadata.
3. **Reference files** — reference FASTA, indexes, annotation catalogs.
4. **Runtime non-RuntimeAttr** — docker image paths, core counts, sharding info.
5. **RuntimeAttr? overrides** — one per task, named to reflect the task's function.
Every workflow takes a `prefix` input, passed to every task for naming output files.

### Task Structure
Tasks are ordered as follows, with each section separated by a single blank line:
1. Inputs.
2. Dynamically generated variables.
3. Command block.
4. Outputs.
5. Runtime — default RuntimeAttr, `select_first` with override, then the `runtime {}` block.

### Task Inputs
- No blank lines between inputs within a task's input section.
- Every task must have `String docker` and `RuntimeAttr? runtime_attr_override` inputs.
- Every task must have a `String prefix` input, used with a descriptive suffix for output file names.
- When calling a task, docker and runtime_attr_override should be explicitly named — e.g., `docker = utils_docker`, `runtime_attr_override = runtime_attr_annotate_svan`.
- Inputs passed to a task call have spaces on either side of `=`.
- The order of inputs at a task call site should match the order defined in the task's input block.

### Command Blocks
- Always start with `set -euo pipefail` followed by a blank line.
- Use `~{var}` for WDL string interpolation (not `${var}`).

### RuntimeAttr Pattern
Every task defines a default `RuntimeAttr` and merges it with the override:
```wdl
RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4,
    disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
}
RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
}
```
Defaults to follow:
- `disk_gb`: Calculate dynamically from input file sizes. Define inline unless the formula is complex.
- `mem_gb`, `boot_disk_gb`, `cpu_cores`: Set based on task compute intensity, not input file size.
- `boot_disk_gb`: Always `10`.
- `preemptible_tries`: Always `2`.
- `max_retries`: Always `0`.
- Use `default_attr` and `runtime_attr` as variable names (not `runtime_default`/`runtime_override`).
- Memory unit is `GiB` (not `GB`).

### Formatting
- No blank comments (e.g., `#########################`).
- No consecutive blank lines — maximum one blank line at a time.
- No extra indentation for alignment purposes. Indentation is only at line start.
- No blank lines between inputs passed to a task call.

### Shared Tasks
Reusable tasks live in `wdl/utils/Helpers.wdl`. Tasks that can be generalized across workflows should be added there rather than duplicated in individual workflow files. Key existing helper tasks: `AddFilter`, `AddInfo`, `BedtoolsClosest`, `ConcatVcfs`, `SubsetVcfToContig`, `SubsetTsvToContig`.

### Annotation Workflow Outputs
Annotation workflows should output a TSV file rather than a VCF, unless underlying tool is specifically designed to annotate VCFs inline.


## Python Conventions
- Lint with `flake8`. CI uses `--max-line-length=130`. The `.flake8` config at repo root ignores `E203` and `W503`.
- Scripts live in `scripts/` and run as standalone CLI tools inside Docker containers. Use `argparse` for argument parsing.
- Use `pysam` for VCF manipulation. For Hail-based workflows, use `hl.import_vcf()` with `force_bgz=True`.
- Do not import across script directories — each script is self-contained within its Docker container.


## Dockerfiles
- Named `Dockerfile.<ToolName>` (e.g., `Dockerfile.Utils`, `Dockerfile.PALMER`).
- Common base images: `ubuntu:22.04` for general-purpose containers, tool-specific bases for specialized containers.
- Docker image strings are never hardcoded in WDL — they are passed as `String` inputs (e.g., `String utils_docker`).
- The `env/` subdirectory in `dockerfiles/` contains conda environment YAML files.


## Change Propagation
When implementing a new annotation or tool:
1. **Scripts** (`scripts/`): Implement core logic if needed - prefer to use in-line Python within workflows though, unless I explicitly ask you to create a separate script.
2. **WDL** (`wdl/`): Create or update workflow and tasks.
3. **Docker** (`dockerfiles/`): Add or update Dockerfile if new dependencies are needed.
4. **Dockstore** (`.dockstore.yml`): Add entry for any new directly-run workflow.
5. **README** (`README.md`): Document new references, output schema fields, or pipeline steps.

## Key Files
- `wdl/utils/Structs.wdl` — Defines the `RuntimeAttr` struct used by every task.
- `wdl/utils/Helpers.wdl` — Shared reusable tasks (filtering, annotation, VCF concat, subsetting).
- `.dockstore.yml` — Registry of all pipeline workflows.
- `.flake8` / `pyproject.toml` — Python formatting and linting config.
- `.github/workflows/wdl-validation.yml` — CI for WDL syntax validation.
- `.github/workflows/python-linting.yaml` — CI for Python linting.


## Things to Avoid
- Do not modify anything in `archive/` unless explicitly asked to — it is retained for historical reference only.
- Do not manually hardcode Docker image URIs in WDL tasks.
- Do not use `${var}` interpolation in WDL command blocks — use `~{var}`.
- Do not add unnecessary edge-case handling or fallback logic.
- Do not create backup or versioned copies of files.
- Do not add blank comment lines or consecutive blank lines in WDL files.
