# Tests
This directory holds the test suite for the scripts and workflows in this repository.


## Layout
- `unit/` - fast pytest tests of the Python under `scripts/` and `.github/scripts/`, plus repository invariants: Dockerfile pins resolve and are used, and every `/opt/scripts` path a workflow calls exists. They need no container and no Hail.
- `unit/test_smoke_coverage.py` - the coverage rule for the smoke tests; see [Coverage rule](#coverage-rule).
- `wdl/` - smoke tests that run workflows under miniwdl in the `utils` image, marked `wdl`. `test_annotation.py` and `test_annotation_utils.py` cover the directly-run workflows in `wdl/annotation` and `wdl/annotation_utils`; each test computes its expected outputs independently in Python from the fixture records.
- `data/` - test fixtures: a 385-record chr22 subset VCF and a small chr22 window FASTA; see [data/README.md](data/README.md) for how they were built.
- `conftest.py` - shared fixtures and helpers, including `run_workflow`, which runs a WDL under miniwdl and returns its outputs.


## Setup
Install the test dependencies in a Python 3.9 or newer environment:

```bash
pip install -e ".[test]"
```


## Unit tests
```bash
pytest tests/unit
```


## WDL smoke tests
The smoke tests need the `utils` image and a container backend. Build the image locally without pushing it:

```bash
dockerfiles/build_docker.sh --no-push utils
```

Then run the tests against it, using podman as the miniwdl container backend. The image tag is the current branch name:

```bash
MINIWDL__SCHEDULER__CONTAINER_BACKEND=podman \
UTILS_DOCKER=localhost/utils:<branch> \
pytest -m wdl tests/wdl
```

- `UTILS_DOCKER` is the image passed to every `*_docker` workflow input. The `wdl` tests are skipped when it is unset.
- With the podman backend the tests set `MINIWDL__PODMAN__EXE` to `["podman"]` unless it is already set, since miniwdl otherwise runs podman through `sudo`.
- Each test runs in its own directory under `tests/.runs/`, which is git-ignored. Inspect it for task logs and outputs after a failure.
- Plain `pytest` collects both suites, and skips the `wdl` tests without `UTILS_DOCKER`.

The full suite takes about six minutes.


## Coverage rule
Every directly-run workflow in `wdl/annotation` and `wdl/annotation_utils` needs a smoke test in `tests/wdl/`. `tests/unit/test_smoke_coverage.py` enforces this: it finds each workflow, collects the `.wdl` files the tests in `tests/wdl/` run, and fails when a workflow is neither tested nor listed in `EXEMPT`, or when a workflow is both.

To add an exemption, add the workflow name to the `EXEMPT` dict in that file with a one-line reason. Use it only for a workflow that needs an image other than `utils` or a large external resource. Remove the entry when you add the test.


## Continuous integration
`.github/workflows/tests.yml` runs on pushes to any branch and on pull requests that touch the workflows, scripts, tests, the `utils` image or `pyproject.toml`.
- The `unit` job installs the test dependencies on Python 3.10 and runs `pytest tests/unit`.
- The `wdl` job builds `dockerfiles/Dockerfile.utils` with the build arguments from `dockerfiles/build_args.sh`, loads it into Docker as `utils:ci`, and runs `pytest -m wdl tests/wdl` with the default miniwdl Docker backend. miniwdl uses an image already present in the local Docker daemon rather than pulling it.
