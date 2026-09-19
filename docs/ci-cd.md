# CI/CD
This document describes the checks that run on every change, the equivalent commands to run locally before pushing, and how workflows and images are released.


## GitHub Actions
Six workflows run on push and pull request to `main`, each gated on the paths it checks. Every check is a Python 3.8 stdlib script - plus `pyyaml` for the Dockstore check and a Java 17 / womtool download for validation - so each can be run locally without CI.

| Workflow | Trigger paths | What it does |
|---|---|---|
| [`wdl-validation.yml`](../.github/workflows/wdl-validation.yml) | `wdl/**` | Downloads womtool 87 and runs `womtool validate` over every `.wdl` file under `wdl/`, failing on the first file that does not validate. |
| [`wdl-style-check.yml`](../.github/workflows/wdl-style-check.yml) | `wdl/**`, `.github/scripts/check_wdl_style.py` | Runs [`check_wdl_style.py`](../.github/scripts/check_wdl_style.py), which enforces the mechanically checkable WDL rules in [conventions.md](conventions.md) - indentation, whitespace, task section order, naming (including call aliases), unused `RuntimeAttr?` inputs, and alphabetical task order in task libraries. `archive/` is not checked. |
| [`markdown-style-check.yml`](../.github/workflows/markdown-style-check.yml) | `**.md`, `.github/scripts/check_markdown_style.py` | Runs [`check_markdown_style.py`](../.github/scripts/check_markdown_style.py), which enforces the Markdown rules in [conventions.md](conventions.md) - whitespace, blank-line spacing around headings, thematic dividers, and nested list indentation. Only `README.md`, `AGENTS.md` and `docs/` are checked; `archive/` and the `AGENTS.md` symlinks are skipped. |
| [`python-linting.yaml`](../.github/workflows/python-linting.yaml) | `scripts/**` | Runs `flake8` over `scripts/`. |
| [`dockstore-sync.yml`](../.github/workflows/dockstore-sync.yml) | `wdl/**`, `.dockstore.yml` | Runs [`check_dockstore_sync.py`](../.github/scripts/check_dockstore_sync.py), which fails if any active workflow in `wdl/annotation`, `wdl/annotation_utils` or `wdl/tools` is missing a `.dockstore.yml` entry (or vice versa), aside from the allowlisted external-repo workflows `AnnotateAF` and `QcAnnotations`. |
| [`agents-sync-check.yml`](../.github/workflows/agents-sync-check.yml) | `AGENTS.md`, `.claude/CLAUDE.md`, `.github/copilot-instructions.md` | Runs [`check_agents_sync.py`](../.github/scripts/check_agents_sync.py), which fails if either instruction mirror has drifted from the canonical `AGENTS.md`. |


### Style checker escape hatch
`check_wdl_style.py` carries an empty `SET_LINE_ALLOWLIST`, keyed by `(path, task name)`, for tasks whose command block genuinely cannot open with `set -euo pipefail`. Add an entry with a comment giving the reason rather than weakening the rule for the whole repository.


## Local Checks
Run the checks that match what you changed before pushing. These mirror the instructions in [`AGENTS.md`](../AGENTS.md).
- Changed anything under `wdl/`: validate syntax with `find wdl -type f -name "*.wdl" -exec java -jar womtool.jar validate {} \;` and style with `python .github/scripts/check_wdl_style.py`.
- Changed any `.md` file: run `python .github/scripts/check_markdown_style.py`.
- Changed anything under `scripts/`: run `flake8 scripts/`. Configuration lives in `.flake8` - maximum line length 130, `E203` and `W503` ignored. `pyproject.toml` pins Black to `line-length = 88` for optional local formatting, but Black is not enforced in CI.
- Added, renamed or deleted a workflow under `wdl/annotation`, `wdl/annotation_utils` or `wdl/tools`: run `python .github/scripts/check_dockstore_sync.py`, which needs `pip install pyyaml`.
- Changed `AGENTS.md`: run `python .github/scripts/check_agents_sync.py`. `.claude/CLAUDE.md` and `.github/copilot-instructions.md` are symlinks to `../AGENTS.md`, so only the canonical file should ever be edited.
- Optional deeper audit: `miniwdl check --strict <file>.wdl` reports unused declarations and name collisions. It is deliberately not part of CI, because it also flags the index localization inputs and the sub-workflow namespace collisions that [conventions.md](conventions.md) documents as intentional.


## Dockstore Registration
Terra imports workflows from Dockstore, which syncs from `main` on push. There is no release tagging step.
- Every directly-run workflow needs an entry in [`.dockstore.yml`](../.dockstore.yml), placed under its matching `# Annotation Workflows`, `# Annotation Utilities` or `# Tools` comment block.
- An entry sets `subclass: WDL`, `name` to the file stem, and `primaryDescriptorPath` to `/wdl/<directory>/<Name>.wdl`, with filters `branches: [main]` and `tags: /.*/`.
- Task libraries and the sub-workflows under `wdl/utils/` are never registered, since they are imported rather than run directly.


## Docker Images
Images are not built by CI. They are built and pushed locally - see [Building and pushing](dockers.md#building-and-pushing) in the Docker images document for the `build_docker.sh` flow.

A [`docker-build-push.yml`](../archive/.github/workflows/docker-build-push.yml) workflow and its [`build_changed_dockers.sh`](../archive/.github/scripts/build_changed_dockers.sh) helper were designed to build changed images on push, but are kept in `archive/` and are intentionally inactive: they require a `GCP_SA_KEY` service-account secret that has not been configured for the repository.
