import re
import subprocess

import pytest
from conftest import REPO_ROOT

DOCKERFILES = sorted((REPO_ROOT / "dockerfiles").glob("Dockerfile.*"))
WDL_FILES = sorted((REPO_ROOT / "wdl").rglob("*.wdl"))

# PALMER builds from the maintainer's fork at its default branch by choice
UNPINNED_CLONES = {"palmer"}


def image_name(dockerfile):
    return dockerfile.name.removeprefix("Dockerfile.")


def instructions(dockerfile):
    return re.sub(r"\\\n", " ", dockerfile.read_text()).splitlines()


@pytest.mark.parametrize("dockerfile", DOCKERFILES, ids=image_name)
def test_build_args_resolve(dockerfile):
    result = subprocess.run(
        ["bash", str(REPO_ROOT / "dockerfiles" / "build_args.sh"), image_name(dockerfile)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_versions_env_keys_are_used():
    declared = {(image_name(d), arg) for d in DOCKERFILES for arg in re.findall(r"^ARG ([A-Z0-9_]+)$", d.read_text(), re.M)}
    keys = re.findall(r"^([a-z0-9-]+)__([A-Z0-9_]+)=", (REPO_ROOT / "dockerfiles" / "versions.env").read_text(), re.M)
    assert sorted(set(keys) - declared) == []


@pytest.mark.parametrize("dockerfile", [d for d in DOCKERFILES if image_name(d) not in UNPINNED_CLONES], ids=image_name)
def test_git_clones_are_pinned(dockerfile):
    unpinned = [
        line for line in instructions(dockerfile) if "git clone" in line and "checkout" not in line and "--branch" not in line
    ]
    assert unpinned == []


def test_wdl_script_paths_exist():
    referenced = {path for wdl in WDL_FILES for path in re.findall(r"/opt/scripts/([\w./-]+)", wdl.read_text())}
    assert referenced
    assert sorted(path for path in referenced if not (REPO_ROOT / "scripts" / path).is_file()) == []


def test_wdl_never_fetches_from_github():
    offenders = [str(wdl.relative_to(REPO_ROOT)) for wdl in WDL_FILES if "raw.githubusercontent.com" in wdl.read_text()]
    assert offenders == []
