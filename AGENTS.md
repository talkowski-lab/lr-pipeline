# Instructions
<!-- Canonical project instructions for Claude Code, Codex, and GitHub Copilot. Edit this file only. -->

- Repo annotates long-read variant callsets (SVs, MEIs, TRs, other complex variants) for HPRC/HGSVC/AoU — cohort/pipeline details in [docs/cohort.md](docs/cohort.md) and [docs/pipeline.md](docs/pipeline.md).
- Stack: WDL 1.0 on Cromwell/Terra, Python 3.8+, Bash, Hail. GCP (`gs://` URIs). Reference genome GRCh38.
- Repo layout, Dockerfile/image conventions, and CI wiring: [docs/repository-structure.md](docs/repository-structure.md).
- WDL/Python style conventions: [docs/conventions.md](docs/conventions.md) — always follow when writing or editing WDL/Python. Match the dominant style of the surrounding file (4-space indent, no tabs, task section order `input` → `command` → `output` → `runtime`); never introduce a second style into a file.
- After changing any file under `wdl/`, run both of these and fix what they report before reporting the change complete:
  - Syntax: `find wdl -type f -name "*.wdl" -exec java -jar womtool.jar validate {} \;` (local jar at `/Users/kjaising/Desktop/Work/Miscellaneous/Software/womtool-87.jar`).
  - Style: `python .github/scripts/check_wdl_style.py` (checks `wdl/`; the same script runs in CI).
- Every workflow must document **every input and every output** in a `parameter_meta` block, in declaration order, plus a `meta { description: [...] }` block — `docs/workflows.md` is generated from these and must never be hand-edited. Rules: [docs/conventions.md](docs/conventions.md#workflow-documentation). Exceptions: never document `RuntimeAttr?`, `*_docker` or `prefix` inputs.
- After changing any file under `scripts/` or `.github/scripts/`, run `flake8 scripts/ .github/scripts/` (config in `.flake8`, max-line-length 130) and fix what it reports.
- After changing any `.md` file, run `python .github/scripts/check_markdown_style.py` (checks `README.md`, `AGENTS.md` and `docs/`; the same script runs in CI) and fix what it reports.
- New directly-run workflows need a `.dockstore.yml` entry — match the format of existing entries.
- New annotation/tool checklist: implement in `wdl/`, adding to `scripts/` only if inline Python in the workflow isn't enough → add/update Dockerfile if new deps needed → register in `.dockstore.yml` → write `meta` and `parameter_meta` covering every input and output → update any other affected `docs/`.
- `archive/` contains retired workflows, scripts, Dockerfiles, and reference documentation. It is not active pipeline code and is excluded from active validation and Dockstore registration; do not use it for new work. Its workflows still carry `meta` and `parameter_meta` blocks, which generate `archive/docs/workflows.md` via `python .github/scripts/generate_workflows_doc.py --site archive`; never hand-edit that document either.
- Don't hardcode Docker image URIs in WDL — always pass as a `String` input.
