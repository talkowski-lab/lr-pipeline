# Contributing
Thanks for helping improve lr-pipeline. This page covers how to propose a change as an outside contributor, the checks to run first, and, for the maintainers, how branches reach `main`. All participation falls under the [Code of Conduct](CODE_OF_CONDUCT.md); report security problems as described in [SECURITY.md](SECURITY.md), not in public issues.


## Proposing a Change
1. Fork the repository and create a branch with any descriptive name other than `main` - e.g. `annotate-dbvar-mem`.
2. Make the change following [Conventions](docs/conventions.md) - WDL, Python and Markdown style, and the comment rules.
3. Run the local checks below that match what you changed.
4. Open a pull request against `main` and fill in the template. A pull request is required: `main` is protected by a ruleset - see [Branch Protection](docs/ci-cd.md#branch-protection).

Do not edit `docs/workflows.md` or `archive/docs/workflows.md`. They are generated from each workflow's `meta` and `parameter_meta` blocks, and CI regenerates and commits them after the merge.


## Local Checks
Install the Python tooling with `pip install -e ".[dev]"`, and run the checks that match what you changed. [CI/CD](docs/ci-cd.md#local-checks) explains each one.

| Changed | Command |
|---|---|
| `wdl/` | `find wdl -type f -name '*.wdl' -print0 \| xargs -0 -n1 java -jar "${WOMTOOL_JAR:-womtool.jar}" validate` |
| `wdl/` | `python .github/scripts/check_wdl_style.py` |
| `wdl/utils/`, `archive/wdl/` | `find archive/wdl -type f -name '*.wdl' \| grep -vxFf archive/wdl/womtool-excluded.txt \| xargs -n1 java -jar "${WOMTOOL_JAR:-womtool.jar}" validate` |
| `scripts/`, `.github/scripts/`, `tests/` | `ruff check scripts .github/scripts tests` and `ruff format --check scripts .github/scripts tests` |
| `scripts/`, `.github/scripts/`, `tests/` | `pytest tests/unit` |
| A workflow in `tests/wdl/`, or `Dockerfile.utils` | `dockerfiles/build_docker.sh --no-push utils`, then `MINIWDL__SCHEDULER__CONTAINER_BACKEND=podman UTILS_DOCKER=localhost/utils:<branch> pytest -m wdl tests/wdl` |
| Shell scripts, `dockerfiles/` | `shellcheck -S warning <script>` and `hadolint dockerfiles/Dockerfile.<image>` |
| Any `.md` file | `python .github/scripts/check_markdown_style.py` |
| `.dockstore.yml` or a directly-run workflow | `python .github/scripts/check_dockstore_sync.py` |
| `AGENTS.md` | `python .github/scripts/check_agents_sync.py` |

Womtool is a single jar: set `WOMTOOL_JAR` to its path, or fetch it into the repository root with `curl -fsSL -o womtool.jar https://github.com/broadinstitute/cromwell/releases/download/87/womtool-87.jar`. It is gitignored.


## Adding a Workflow
- Implement the workflow in `wdl/`, adding a script under `scripts/` only if inline Python in the workflow is not enough.
- Add or update a Dockerfile if new dependencies are needed.
- Give every directly-run workflow an entry in [`.dockstore.yml`](.dockstore.yml), matching the existing entries and listing `main` - see [Dockstore Registration](docs/ci-cd.md#dockstore-registration).
- Write a `meta` block and a `parameter_meta` block covering every input and output - see [Workflow documentation](docs/conventions.md#workflow-documentation).
- Pass Docker images and reference files in as inputs rather than hardcoding them, and update any other affected page in `docs/`.


## Docker Images
Images are built from [`dockerfiles/`](dockerfiles/) with `dockerfiles/build_docker.sh`. Pushing needs access to the maintainers' registry, so contributors build locally with `dockerfiles/build_docker.sh --no-push <image-name>`, which produces `localhost/<image-name>:<branch>`. Pin new tool versions or upstream commits in `dockerfiles/versions.env`. See [Docker images](docs/dockers.md#building-and-pushing).


## Maintainer Workflow
Maintainers work from a worktree per branch rather than from a fork, and merge locally rather than through a pull request.
- `.github/scripts/new_worktree.sh <branch>` cuts the branch from `origin/main`, pushes it and creates `../lr-pipeline-worktrees/<branch>`.
- Every push to the branch runs the CI workflows whose paths it touches, so check those runs before merging.
- To test a workflow from the branch in Terra, list the branch under `filters.branches` in its `.dockstore.yml` entry.
- `.github/scripts/merge_branch.sh <branch>` strips the branch filters, rebases, fast-forwards `main`, pushes, and removes the worktree, branch and branch image tags. Repository administrators are exempt from the pull request ruleset, which is what allows this push.

The full sequence is in [Feature branch test versions](docs/ci-cd.md#feature-branch-test-versions).


## Agent-Assisted Development
[AGENTS.md](AGENTS.md) is the canonical instruction file for coding agents (Claude Code, Codex and GitHub Copilot); `.claude/CLAUDE.md` and `.github/copilot-instructions.md` are symlinks to it, so edit only `AGENTS.md`. The hooks in `.github/hooks/` and `.githooks/pre-commit` keep agent sessions out of the main checkout - see [Branches](docs/conventions.md#branches). Agent-written changes meet the same checks and review as any other.
