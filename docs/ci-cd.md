# CI/CD
This document describes the checks that run on every change, the equivalent commands to run locally before pushing, and how workflows and images are released.


## GitHub Actions
Ten workflows run on push and pull request to `main`, each gated on the paths it checks; `dockstore-published.yml` also runs weekly. The style and sync checks are Python stdlib scripts - plus `pyyaml` for the Dockstore check and a Java 17 / womtool download for validation - so each can be run locally without CI. Every action is pinned by commit SHA and every workflow runs with `permissions: contents: read`; `workflows-doc-sync.yml` is the one workflow with a job granted `contents: write`. Each workflow has a concurrency group that cancels a superseded pull request run, but never a run on `main`. Python is 3.10 throughout. Every job has a `timeout-minutes` limit, each workflow also re-runs when its own file changes, and `tests.yml` and `wdl-validation.yml` can be re-run manually with `workflow_dispatch`. Dependabot ([`dependabot.yml`](../.github/dependabot.yml)) proposes weekly updates for the GitHub Actions, pip dependencies and Docker base images, ignoring the private registry base images.

| Workflow | Trigger paths | What it does |
|---|---|---|
| [`wdl-validation.yml`](../.github/workflows/wdl-validation.yml) | `wdl/**` | Downloads womtool 87 (cached between runs) and runs `womtool validate` over every `.wdl` file under `wdl/` in parallel, failing if any file does not validate. |
| [`wdl-style-check.yml`](../.github/workflows/wdl-style-check.yml) | `wdl/**`, `.github/scripts/check_wdl_style.py`, `.github/scripts/wdl_meta.py` | Runs [`check_wdl_style.py`](../.github/scripts/check_wdl_style.py), which enforces the mechanically checkable WDL rules in [conventions.md](conventions.md) - indentation, whitespace, task section order, naming (including call aliases), unused `RuntimeAttr?` inputs, undefined values entering a nested scatter, alphabetical task order in task libraries, and the `meta` and `parameter_meta` rules that [workflows.md](workflows.md) is generated from. `archive/` is not checked. |
| [`workflows-doc-sync.yml`](../.github/workflows/workflows-doc-sync.yml) | `wdl/**`, `archive/wdl/**`, `.dockstore.yml`, `.github/scripts/wdl_meta.py`, `.github/scripts/generate_workflows_doc.py` | On a pull request, checks the WDL style rules and confirms [workflows.md](workflows.md) and its archive counterpart can be generated. On push to `main`, regenerates both with [`generate_workflows_doc.py`](../.github/scripts/generate_workflows_doc.py) and commits the result when it changed. It pushes over SSH with a write deploy key stored as the `DOCS_DEPLOY_KEY` secret, the one actor allowed to bypass the pull-request ruleset - see [Branch Protection](#branch-protection); the commit touches only the generated documents, which are outside the workflow's trigger paths, so it cannot loop. |
| [`markdown-style-check.yml`](../.github/workflows/markdown-style-check.yml) | `**.md`, `.github/scripts/check_markdown_style.py` | Runs [`check_markdown_style.py`](../.github/scripts/check_markdown_style.py), which enforces the Markdown rules in [conventions.md](conventions.md) - whitespace, blank-line spacing around headings, thematic dividers, and nested list indentation. `README.md`, `AGENTS.md`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `SECURITY.md`, `.github/PULL_REQUEST_TEMPLATE.md`, `docs/` and `tests/` are checked; `archive/` and the `AGENTS.md` symlinks are skipped. |
| [`python-lint.yml`](../.github/workflows/python-lint.yml) | `scripts/**`, `.github/scripts/**`, `tests/**`, `pyproject.toml` | Runs `ruff check` and `ruff format --check` over `scripts`, `.github/scripts` and `tests`, with ruff pinned to the version in the workflow. Configuration lives in `pyproject.toml`. |
| [`tests.yml`](../.github/workflows/tests.yml) | `wdl/**`, `scripts/**`, `tests/**`, `.github/scripts/**`, `dockerfiles/Dockerfile.utils`, `dockerfiles/versions.env`, `dockerfiles/build_args.sh`, `pyproject.toml` | Two jobs. `unit` runs `pytest tests/unit`, the fast checks of script functions and the `.github/scripts` checkers. `wdl` builds the utils image from `Dockerfile.utils` with a GitHub Actions layer cache - the Artifact Registry is private, so CI cannot pull it - and runs the `miniwdl` smoke tests in `tests/wdl/` against it with Docker. See [Tests](#tests). |
| [`shell-lint.yml`](../.github/workflows/shell-lint.yml) | `**/*.sh`, `.githooks/**`, `dockerfiles/**`, `.hadolint.yaml` | Two jobs. `shellcheck` runs `shellcheck -S warning` over `dockerfiles/*.sh`, `.github/scripts/*.sh`, `.github/hooks/*.sh`, `tests/*.sh` and `.githooks/pre-commit`. `hadolint` lints every `dockerfiles/Dockerfile.*` using [`.hadolint.yaml`](../.hadolint.yaml), which fails on warnings and carries a justified comment for each ignored rule. |
| [`docker-build.yml`](../.github/workflows/docker-build.yml) | `dockerfiles/**`, manual dispatch with an `images` input | Builds `utils` as `utils:ci`, then every changed Dockerfile on top of it with `--build-arg BASE_IMAGE=utils:ci` and the arguments from [`build_args.sh`](../dockerfiles/build_args.sh), so a broken pin or install fails the pull request. A change to `versions.env` rebuilds the images whose keys changed. Never pushes. |
| [`dockstore-sync.yml`](../.github/workflows/dockstore-sync.yml) | `wdl/**`, `.dockstore.yml` | Runs [`check_dockstore_sync.py`](../.github/scripts/check_dockstore_sync.py), which fails if any active workflow in `wdl/annotation`, `wdl/annotation_utils` or `wdl/tools` is missing a `.dockstore.yml` entry (or vice versa), if an entry does not set `publish: true`, or if an entry does not list `main` under `filters.branches`. On push to `main` only, it runs again with `--main-only`, which fails if a feature branch is still listed. |
| [`dockstore-published.yml`](../.github/workflows/dockstore-published.yml) | `.dockstore.yml`, weekly schedule, manual dispatch | Runs [`check_dockstore_published.py`](../.github/scripts/check_dockstore_published.py), which asks the public Dockstore API whether each `.dockstore.yml` entry is published and fails listing those that are not, since Terra cannot import an unpublished workflow. Network errors and 5xx responses are retried. On push to `main` it first waits five minutes for Dockstore to sync, then `--added-since <before-sha>` reports entries added by that push as warnings in case the sync is still running; the weekly and manual runs are strict. |
| [`agents-sync-check.yml`](../.github/workflows/agents-sync-check.yml) | `AGENTS.md`, `.claude/CLAUDE.md`, `.github/copilot-instructions.md` | Runs [`check_agents_sync.py`](../.github/scripts/check_agents_sync.py), which fails if either instruction mirror has drifted from the canonical `AGENTS.md`. |


### Style checker escape hatch
`check_wdl_style.py` carries an empty `SET_LINE_ALLOWLIST`, keyed by `(path, task name)`, for tasks whose command block genuinely cannot open with `set -euo pipefail`. Add an entry with a comment giving the reason rather than weakening the rule for the whole repository.


## Local Checks
Run the checks that match what you changed before pushing. These mirror the instructions in [`AGENTS.md`](../AGENTS.md).
- Changed anything under `wdl/`: validate syntax with `find wdl -type f -name '*.wdl' -print0 | xargs -0 -n1 java -jar womtool.jar validate` and style with `python .github/scripts/check_wdl_style.py`. Fetch the jar once into the main checkout's root with `curl -fsSL -o womtool.jar https://github.com/broadinstitute/cromwell/releases/download/87/womtool-87.jar`; it is gitignored and [`.worktreelinks`](../.worktreelinks) shares it into every worktree.
- Changed any `.md` file: run `python .github/scripts/check_markdown_style.py`. `docs/workflows.md` and `archive/docs/workflows.md` are generated rather than edited - see [Workflow documentation](conventions.md#workflow-documentation); `python .github/scripts/generate_workflows_doc.py` (with `--site archive` for the archive) previews the result locally.
- Changed anything under `scripts/`, `.github/scripts/` or `tests/`: run `ruff check scripts .github/scripts tests`, `ruff format --check scripts .github/scripts tests` and `pytest tests/unit`. `ruff format` applies the formatting. Configuration lives in `pyproject.toml` - line length 130, Python 3.9 target, rule sets `E`, `F`, `W`, `I` and `B`. Install the tooling with `pip install -e ".[dev]"`.
- Changed a workflow covered by `tests/wdl/`, or `dockerfiles/Dockerfile.utils`: run the smoke tests described under [Tests](#tests).
- Changed a shell script or a Dockerfile: run `shellcheck -S warning` on the script, or `hadolint` on the Dockerfile from the repository root so `.hadolint.yaml` is picked up.
- Added, renamed or deleted a workflow under `wdl/annotation`, `wdl/annotation_utils` or `wdl/tools`, or edited `.dockstore.yml`: run `python .github/scripts/check_dockstore_sync.py`, which needs `pip install pyyaml`. Before merging a feature branch into `main`, run it again with `--main-only` to confirm no branch filter was left behind.
- Changed `AGENTS.md`: run `python .github/scripts/check_agents_sync.py`. `.claude/CLAUDE.md` and `.github/copilot-instructions.md` are symlinks to `../AGENTS.md`, so only the canonical file should ever be edited.
- Optional deeper audit: `miniwdl check --strict <file>.wdl` reports unused declarations and name collisions. It is deliberately not part of CI, because it also flags the index localization inputs and the sub-workflow namespace collisions that [conventions.md](conventions.md) documents as intentional.


## Tests
Two layers of tests live under [`tests/`](../tests/), both configured in `pyproject.toml`; install their dependencies with `pip install -e ".[test]"`.
- `tests/unit/` - fast `pytest` tests of the functions in `scripts/` and of the `.github/scripts/` checkers. They need no containers: `pytest tests/unit`.
- `tests/wdl/` - `miniwdl` smoke tests, marked `wdl`, that run several utils-image workflows on `tests/data/hprc_hgsvc_v10.chr22.subset.vcf.gz`, an open-access subset of the HPRC/HGSVC V10 callset with 385 records on chr22 and 5 samples. They need a local utils image:
  ```bash
  dockerfiles/build_docker.sh --no-push utils
  MINIWDL__SCHEDULER__CONTAINER_BACKEND=podman UTILS_DOCKER=localhost/utils:<branch> pytest -m wdl tests/wdl
  ```
  CI uses Docker and builds the utils image itself from `Dockerfile.utils`, because the Artifact Registry is private.

See [`tests/README.md`](../tests/README.md) for what each test covers.


## Dockstore Registration
Terra imports workflows from Dockstore, which syncs on push from every branch listed in an entry's `filters.branches`. There is no release tagging step.
- Every directly-run workflow needs an entry in [`.dockstore.yml`](../.dockstore.yml), placed under its matching `# Annotation Workflows`, `# Annotation Utilities` or `# Tools` comment block.
- An entry sets `subclass: WDL`, `name` to the file stem, `primaryDescriptorPath` to `/wdl/<directory>/<Name>.wdl` and `publish: true`, with filters `branches: [main]` - plus any feature branch currently under test - and `tags: /.*/`.
- `main` is what Terra runs in production and must be listed in every entry at all times.
- `publish: true` makes Dockstore publish the entry when it next syncs, so no manual step is needed; [`check_dockstore_published.py`](../.github/scripts/check_dockstore_published.py) reports any entry the public API still shows as unpublished.
- Task libraries and the sub-workflows under `wdl/utils/` are never registered, since they are imported rather than run directly.


### Feature branch test versions
Development happens on a branch rather than on `main`, named with any descriptive name other than `main`, in its own worktree created by [`new_worktree.sh`](../.github/scripts/new_worktree.sh) - see [Branches](conventions.md#branches). To test a workflow the branch changes - directly, or through a task library or `wdl/utils/` sub-workflow it imports - add the branch name under `filters.branches` in that workflow's entry, after `main`, and push; Dockstore reads the `.dockstore.yml` on the pushed branch, so the version appears under the branch name and can be imported into Terra. Only the workflows being tested should list the branch.

This is the maintainer workflow; external contributors use a fork and a pull request instead - see [CONTRIBUTING.md](../CONTRIBUTING.md). Once the changes are validated and a merge has been explicitly approved, the maintainer merges locally - no pull request - with [`merge_branch.sh`](../.github/scripts/merge_branch.sh), run from the main checkout, which must be on `main`, and with the branch's worktree clean.
```bash
.github/scripts/merge_branch.sh my-topic # merges a named branch from its worktree
.github/scripts/merge_branch.sh          # run inside a worktree, merges its branch
```
The script performs the whole sequence, and each step is what to do by hand if it is ever run one piece at a time.
1. Removes the branch from every `filters.branches` list it was added to and commits that, after `check_dockstore_sync.py --main-only` confirms nothing else was left behind. Doing this before the merge keeps the filter off `main` entirely.
2. Rebases the branch, in its worktree, onto `origin/main`. A conflict stops the script with the rebase still in progress and nothing merged, because conflicts are resolved by a human, never automatically. Resolve them, `git rebase --continue`, and run the script again, or `git rebase --abort` to back out.
3. Fast-forwards `main` in the main checkout and pushes it. The merge is `--ff-only`, so it fails rather than quietly creating a merge commit if `main` moved in between.
4. Removes the worktree, then deletes the branch on the remote and locally. Deleting the remote branch is what removes its Dockstore version. The rebased branch is never force-pushed, since it is deleted moments later anyway.
5. Deletes every `<image>:<branch>` tag the branch's builds pushed to Artifact Registry - see [Building and Pushing](dockers.md#building-and-pushing). This runs last, so a `gcloud` failure leaves the git state already final. Any image the branch changed should then be rebuilt from the main checkout, which is what moves `:latest` onto the merged code.

The merge is a push to `main`, which the pull request ruleset in [Branch Protection](#branch-protection) allows for repository administrators. The repository also has `delete_branch_on_merge` enabled, which covers the occasional pull request; a local merge is not a merged pull request as far as GitHub is concerned, so the script deletes the branch explicitly. There is deliberately no CI job that bulk-deletes merged branches, since the repository carries long-lived collaborator branches that such a job would remove.


## Branch Protection
Two GitHub rulesets protect `main`.
- `main: no deletion or force push` applies to everyone, with no bypass.
- `main: pull request required` requires every change to arrive through a pull request, except from repository administrators and deploy keys. Those two exceptions keep the maintainer's fast-forward push from `merge_branch.sh` and the `workflows-doc-sync.yml` bot commit working; GitHub does not accept the GitHub Actions integration as a bypass actor on a repository ruleset, so that job pushes with the `workflows-doc-sync` deploy key, whose private half is the `DOCS_DEPLOY_KEY` repository secret. Rotating the key means replacing both.

No status checks are required, because every workflow is path-filtered: a required check that never runs on a pull request touching other paths would block it indefinitely. Reviewers check the results of the workflows that did run. Private vulnerability reporting is enabled - see [SECURITY.md](../SECURITY.md).


## Docker Images
CI builds images only to check that they build: `docker-build.yml` builds the changed Dockerfiles from the pull request's tree and discards them, so a broken pin or install is caught before release. It uses Docker, builds `utils` first as `utils:ci` (sharing the `tests.yml` layer cache) and builds every other image `FROM` it through `BASE_IMAGE`, with `trgt:ci` first for `trgtlps`. Images that are released are built locally with `build_docker.sh` and pushed to Artifact Registry by the maintainer - see [Building and Pushing](dockers.md#building-and-pushing). Contributors without registry access build with `--no-push`.

[`.hadolint.yaml`](../.hadolint.yaml) configures `hadolint` for the Dockerfiles, run in CI by `shell-lint.yml`; each ignored rule carries a comment justifying it.

A [`docker-build-push.yml`](../archive/.github/workflows/docker-build-push.yml) workflow and its [`build_changed_dockers.sh`](../archive/.github/scripts/build_changed_dockers.sh) helper were designed to build changed images on push, but are kept in `archive/` and are intentionally inactive: they require a `GCP_SA_KEY` service-account secret that has not been configured for the repository.
