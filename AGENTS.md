# Instructions

- Repo annotates long-read variant callsets (SVs, MEIs, TRs, other complex variants) for HPRC/HGSVC/AoU — cohort/pipeline details in [docs/cohort.md](docs/cohort.md) and [docs/pipeline.md](docs/pipeline.md).
- Stack: WDL 1.0 on Cromwell/Terra, Python 3.8+, Bash, Hail. GCP (`gs://` URIs). Reference genome GRCh38.
- Repo layout, Dockerfile/image conventions, and CI wiring: [docs/repository-structure.md](docs/repository-structure.md).
- WDL/Python style conventions: [docs/conventions.md](docs/conventions.md) — always follow when writing or editing WDL/Python.
- Validate WDL before committing: `find wdl -type f -name "*.wdl" -exec java -jar womtool.jar validate {} \;` (local jar at `/Users/kjaising/Desktop/Work/Miscellaneous/Software/womtool-87.jar`).
- Lint Python before committing: `flake8 scripts/` (config in `.flake8`, max-line-length 130).
- New directly-run workflows need a `.dockstore.yml` entry — match the format of existing entries.
- New annotation/tool checklist: implement in `wdl/` (add to `scripts/` only if inline Python in the workflow isn't enough) → add/update Dockerfile if new deps → register in `.dockstore.yml` → document in `docs/`.
- Don't modify `archive/` or `data/` (gitignored, out of scope) unless explicitly asked.
- Don't hardcode Docker image URIs in WDL — always pass as a `String` input.
