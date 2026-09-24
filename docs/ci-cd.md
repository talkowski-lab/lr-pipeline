# CI/CD
This document describes the checks that run on every change, the equivalent commands to run locally before pushing, and how workflows and images are released.


## GitHub Actions
Seven workflows run on push and pull request to `main`, each gated on the paths it checks. Every check is a Python 3.8 stdlib script - plus `pyyaml` for the Dockstore check and a Java 17 / womtool download for validation - so each can be run locally without CI. Six are read-only checkers; `workflows-doc-sync.yml` is the one workflow that writes to the repository.

| Workflow | Trigger paths | What it does |
|---|---|---|
| [`wdl-validation.yml`](../.github/workflows/wdl-validation.yml) | `wdl/**` | Downloads womtool 87 and runs `womtool validate` over every `.wdl` file under `wdl/`, failing on the first file that does not validate. |
| [`wdl-style-check.yml`](../.github/workflows/wdl-style-check.yml) | `wdl/**`, `.github/scripts/check_wdl_style.py`, `.github/scripts/wdl_meta.py` | Runs [`check_wdl_style.py`](../.github/scripts/check_wdl_style.py), which enforces the mechanically checkable WDL rules in [conventions.md](conventions.md) - indentation, whitespace, task section order, naming (including call aliases), unused `RuntimeAttr?` inputs, alphabetical task order in task libraries, and the `meta` and `parameter_meta` rules that [workflows.md](workflows.md) is generated from. `archive/` is not checked. |
| [`workflows-doc-sync.yml`](../.github/workflows/workflows-doc-sync.yml) | `wdl/**`, `archive/wdl/**`, `.dockstore.yml`, `.github/scripts/wdl_meta.py`, `.github/scripts/generate_workflows_doc.py` | On a pull request, checks the WDL style rules and confirms [workflows.md](workflows.md) and its archive counterpart can be generated. On push to `main`, regenerates both with [`generate_workflows_doc.py`](../.github/scripts/generate_workflows_doc.py) and commits the result when it changed. This is the only job granted `contents: write`; because the push uses `GITHUB_TOKEN`, it does not trigger another run. |
| [`markdown-style-check.yml`](../.github/workflows/markdown-style-check.yml) | `**.md`, `.github/scripts/check_markdown_style.py` | Runs [`check_markdown_style.py`](../.github/scripts/check_markdown_style.py), which enforces the Markdown rules in [conventions.md](conventions.md) - whitespace, blank-line spacing around headings, thematic dividers, and nested list indentation. Only `README.md`, `AGENTS.md` and `docs/` are checked; `archive/` and the `AGENTS.md` symlinks are skipped. |
| [`python-linting.yaml`](../.github/workflows/python-linting.yaml) | `scripts/**`, `.github/scripts/**` | Runs `flake8` over `scripts/` and `.github/scripts/`. |
| [`dockstore-sync.yml`](../.github/workflows/dockstore-sync.yml) | `wdl/**`, `.dockstore.yml` | Runs [`check_dockstore_sync.py`](../.github/scripts/check_dockstore_sync.py), which fails if any active workflow in `wdl/annotation`, `wdl/annotation_utils` or `wdl/tools` is missing a `.dockstore.yml` entry (or vice versa), or if an entry does not list `main` under `filters.branches`. On push to `main` only, it runs again with `--main-only`, which fails if a feature branch is still listed. |
| [`agents-sync-check.yml`](../.github/workflows/agents-sync-check.yml) | `AGENTS.md`, `.claude/CLAUDE.md`, `.github/copilot-instructions.md` | Runs [`check_agents_sync.py`](../.github/scripts/check_agents_sync.py), which fails if either instruction mirror has drifted from the canonical `AGENTS.md`. |


### Style checker escape hatch
`check_wdl_style.py` carries an empty `SET_LINE_ALLOWLIST`, keyed by `(path, task name)`, for tasks whose command block genuinely cannot open with `set -euo pipefail`. Add an entry with a comment giving the reason rather than weakening the rule for the whole repository.


## Local Checks
Run the checks that match what you changed before pushing. These mirror the instructions in [`AGENTS.md`](../AGENTS.md).
- Changed anything under `wdl/`: validate syntax with `find wdl -type f -name "*.wdl" -exec java -jar womtool.jar validate {} \;` and style with `python .github/scripts/check_wdl_style.py`.
- Changed any `.md` file: run `python .github/scripts/check_markdown_style.py`. `docs/workflows.md` is generated, so edit the `meta` and `parameter_meta` blocks of the workflow instead; `python .github/scripts/generate_workflows_doc.py` previews the result locally, and CI regenerates and commits it on push to `main`. `archive/docs/workflows.md` is generated the same way from the retired workflows, with `--site archive`.
- Changed anything under `scripts/` or `.github/scripts/`: run `flake8 scripts/ .github/scripts/`. Configuration lives in `.flake8` - maximum line length 130, `E203` and `W503` ignored. `pyproject.toml` pins Black to `line-length = 88` for optional local formatting, but Black is not enforced in CI.
- Added, renamed or deleted a workflow under `wdl/annotation`, `wdl/annotation_utils` or `wdl/tools`, or edited `.dockstore.yml`: run `python .github/scripts/check_dockstore_sync.py`, which needs `pip install pyyaml`. Before merging a feature branch into `main`, run it again with `--main-only` to confirm no branch filter was left behind.
- Changed `AGENTS.md`: run `python .github/scripts/check_agents_sync.py`. `.claude/CLAUDE.md` and `.github/copilot-instructions.md` are symlinks to `../AGENTS.md`, so only the canonical file should ever be edited.
- Optional deeper audit: `miniwdl check --strict <file>.wdl` reports unused declarations and name collisions. It is deliberately not part of CI, because it also flags the index localization inputs and the sub-workflow namespace collisions that [conventions.md](conventions.md) documents as intentional.


## Dockstore Registration
Terra imports workflows from Dockstore, which syncs on push from every branch listed in an entry's `filters.branches`. There is no release tagging step.
- Every directly-run workflow needs an entry in [`.dockstore.yml`](../.dockstore.yml), placed under its matching `# Annotation Workflows`, `# Annotation Utilities` or `# Tools` comment block.
- An entry sets `subclass: WDL`, `name` to the file stem, and `primaryDescriptorPath` to `/wdl/<directory>/<Name>.wdl`, with filters `branches: [main]` - plus any feature branch currently under test - and `tags: /.*/`.
- `main` is what Terra runs in production and must be listed in every entry at all times.
- Task libraries and the sub-workflows under `wdl/utils/` are never registered, since they are imported rather than run directly.


### Feature branch test versions
Development happens on a `kj-<kebab-topic>` branch rather than on `main`, in its own worktree created by [`new_worktree.sh`](../.github/scripts/new_worktree.sh) - see [Branches](conventions.md#branches). To test a changed workflow from that branch, add the branch name under `filters.branches` in that workflow's entry and push; Dockstore reads the `.dockstore.yml` on the pushed branch, so the version appears under the branch name and can be imported into Terra. Only the workflows being tested should list the branch.

Once the changes are validated and a merge has been explicitly approved, merge locally - no pull request - with [`merge_branch.sh`](../.github/scripts/merge_branch.sh), run from the main checkout, which must be on `main`, and with the branch's worktree clean.
```bash
.github/scripts/merge_branch.sh kj-topic # merges a named branch from its worktree
.github/scripts/merge_branch.sh          # run inside a worktree, merges its branch
```
The script performs the whole sequence, and each step is what to do by hand if it is ever run one piece at a time.
1. Removes the branch from every `filters.branches` list it was added to and commits that, after `check_dockstore_sync.py --main-only` confirms nothing else was left behind. Doing this before the merge keeps the filter off `main` entirely.
2. Rebases the branch, in its worktree, onto `origin/main`. A conflict stops the script with the rebase still in progress and nothing merged, because conflicts are resolved by a human, never automatically. Resolve them, `git rebase --continue`, and run the script again, or `git rebase --abort` to back out.
3. Fast-forwards `main` in the main checkout and pushes it. The merge is `--ff-only`, so it fails rather than quietly creating a merge commit if `main` moved in between.
4. Removes the worktree, then deletes the branch on the remote and locally. Deleting the remote branch is what removes its Dockstore version. The rebased branch is never force-pushed, since it is deleted moments later anyway.
5. Deletes every `<image>:<branch>` tag the branch's builds pushed to Artifact Registry - see [Building and Pushing](dockers.md#building-and-pushing). This runs last, so a `gcloud` failure leaves the git state already final. Any image the branch changed should then be rebuilt from the main checkout, which is what moves `:latest` onto the merged code.

The repository also has `delete_branch_on_merge` enabled, which covers the occasional pull request; a local merge is not a merged pull request as far as GitHub is concerned, so the script deletes the branch explicitly. There is deliberately no CI job that bulk-deletes merged branches, since the repository carries long-lived collaborator branches that such a job would remove.


## Docker Images
Images are not built by CI, and are tagged by the branch they are built on - `kj_V<N>` and `:latest` from `main`, `:<branch>` from a feature branch, per [Building and Pushing](dockers.md#building-and-pushing). They are built and pushed locally - see [Building and pushing](dockers.md#building-and-pushing) in the Docker images document for the `build_docker.sh` flow.

A [`docker-build-push.yml`](../archive/.github/workflows/docker-build-push.yml) workflow and its [`build_changed_dockers.sh`](../archive/.github/scripts/build_changed_dockers.sh) helper were designed to build changed images on push, but are kept in `archive/` and are intentionally inactive: they require a `GCP_SA_KEY` service-account secret that has not been configured for the repository.
