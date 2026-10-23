# Tests
This directory holds the test suite for the scripts and workflows in this repository.


## Layout
- `unit/` - fast pytest tests of the Python under `scripts/` and `.github/scripts/`, plus repository invariants: Dockerfile pins resolve and are used, and every `/opt/scripts` path a workflow calls exists. They need no container and no Hail.
- `wdl/` - smoke tests that run workflows under miniwdl in the `utils` image, marked `wdl`.
- `data/` - test fixtures; see [data/README.md](data/README.md) for how the subset VCF was built.
- `conftest.py` - shared fixtures, including `run_workflow`, which runs a WDL under miniwdl and returns its outputs.


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


## Continuous integration
`.github/workflows/tests.yml` runs on pushes and pull requests to `main` that touch the workflows, scripts, tests, the `utils` image or `pyproject.toml`.
- The `unit` job installs the test dependencies on Python 3.10 and runs `pytest tests/unit`.
- The `wdl` job builds `dockerfiles/Dockerfile.utils` with the build arguments from `dockerfiles/build_args.sh`, loads it into Docker as `utils:ci`, and runs `pytest -m wdl tests/wdl` with the default miniwdl Docker backend. miniwdl uses an image already present in the local Docker daemon rather than pulling it.
