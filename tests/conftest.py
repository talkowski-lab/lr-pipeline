import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RUNS_DIR = REPO_ROOT / "tests" / ".runs"
SUBSET_VCF = REPO_ROOT / "tests" / "data" / "hprc_hgsvc_v10.chr22.subset.vcf.gz"

# The checker scripts import one another by bare module name
sys.path.insert(0, str(REPO_ROOT / ".github" / "scripts"))


def load_script(rel_path):
    """Import a script under the repo root by path, without executing its __main__ block."""
    path = REPO_ROOT / rel_path
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def repo_root():
    return REPO_ROOT


@pytest.fixture(scope="session")
def subset_vcf():
    return SUBSET_VCF


@pytest.fixture(scope="session")
def utils_docker():
    image = os.environ.get("UTILS_DOCKER")
    if not image:
        pytest.skip("UTILS_DOCKER is not set")
    return image


class WorkflowError(Exception):
    pass


@pytest.fixture
def run_dir(request):
    """Return a clean per-test directory under tests/.runs/, which is git-ignored."""
    path = RUNS_DIR / re.sub(r"[^\w.-]", "_", request.node.name)
    shutil.rmtree(path, ignore_errors=True)
    path.mkdir(parents=True)
    return path


@pytest.fixture
def run_workflow(run_dir):
    """Return a function that runs a WDL under miniwdl and returns its outputs dict."""

    runs = []

    def run(wdl_rel_path, inputs):
        runs.append(None)
        inputs_path = run_dir / f"inputs{len(runs)}.json"
        inputs_path.write_text(json.dumps(inputs, indent=2))
        # A trailing "/." makes miniwdl run directly in this directory instead of a timestamped subdirectory
        command = [
            "miniwdl",
            "run",
            str(REPO_ROOT / wdl_rel_path),
            "-i",
            str(inputs_path),
            "--dir",
            f"{run_dir}/run{len(runs)}/.",
            "--no-cache",
        ]
        result = subprocess.run(command, capture_output=True, text=True, env=_miniwdl_env())
        if result.returncode != 0:
            raise WorkflowError(f"miniwdl exited {result.returncode}\n{result.stderr[-4000:]}")
        return json.loads(result.stdout)["outputs"]

    return run


def _miniwdl_env():
    env = dict(os.environ)
    # miniwdl runs podman through sudo and always pulls, so use a wrapper that is rootless and keeps locally built images
    if env.get("MINIWDL__SCHEDULER__CONTAINER_BACKEND") == "podman":
        env.setdefault("MINIWDL__PODMAN__EXE", json.dumps([str(REPO_ROOT / "tests" / "podman_local.sh")]))
    return env
