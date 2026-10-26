import json
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


PIN_TOKEN = re.compile(
    r"^[A-Za-z0-9_.\[\]-]+==(\$\{[A-Z0-9_]+\}|\d[\w.+!-]*)$"
    r"|^git\+https?://\S+@(\$\{[A-Z0-9_]+_COMMIT\}|[0-9a-f]{40})$"
)
PIP_INSTALL = re.compile(r"(?:pip3?|-m pip) install\b|python3? \S*get-pip\.py")


def unpinned_pip_packages(dockerfile):
    unpinned = []
    for line in instructions(dockerfile):
        for segment in re.split(r"&&|;|\|", line):
            match = PIP_INSTALL.search(segment)
            if not match:
                continue
            for token in segment[match.end() :].split():
                token = token.strip("\"'")
                if token and not token.startswith("-") and not PIN_TOKEN.match(token):
                    unpinned.append(token)
    return unpinned


@pytest.mark.parametrize("dockerfile", [d for d in DOCKERFILES if image_name(d) not in UNPINNED_CLONES], ids=image_name)
def test_pip_installs_are_pinned(dockerfile):
    assert unpinned_pip_packages(dockerfile) == []


@pytest.mark.parametrize("dockerfile", DOCKERFILES, ids=image_name)
def test_no_apt_upgrade_or_masked_failures(dockerfile):
    offenders = [line for line in instructions(dockerfile) if re.search(r"apt(-get)? (dist-)?upgrade|\|\| true", line)]
    assert offenders == []


@pytest.mark.parametrize("dockerfile", DOCKERFILES, ids=image_name)
def test_derived_dockerfiles_use_base_image_arg(dockerfile):
    text = dockerfile.read_text()
    assert not re.search(r"^FROM us-central1-docker\.pkg\.dev/", text, re.M)
    if re.search(r"^ARG BASE_IMAGE=", text, re.M):
        assert re.search(r"^ARG BASE_IMAGE=us-central1-docker\.pkg\.dev/\S+:latest\n+FROM \$\{BASE_IMAGE\}$", text, re.M)


NOTEBOOKS = [
    REPO_ROOT / path for path in subprocess.check_output(["git", "ls-files", "*.ipynb"], cwd=REPO_ROOT, text=True).split()
]


@pytest.mark.parametrize("notebook", NOTEBOOKS, ids=lambda path: path.name)
def test_notebooks_hold_only_code(notebook):
    cells = json.loads(notebook.read_text())["cells"]
    assert [cell for cell in cells if cell["metadata"] or cell.get("outputs") or cell.get("execution_count") is not None] == []
