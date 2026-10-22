import sys

import pytest
from conftest import load_script

sync = load_script(".github/scripts/check_dockstore_sync.py")

ENTRY = """\
  - subclass: WDL
    name: {name}
    primaryDescriptorPath: {path}
    filters:
      branches:
{branches}
"""


def entry(name, path=None, branches=("main",)):
    path = path or f"/wdl/tools/{name}.wdl"
    return ENTRY.format(name=name, path=path, branches="\n".join(f"        - {b}" for b in branches))


@pytest.fixture
def fake_repo(tmp_path, monkeypatch):
    (tmp_path / "wdl" / "tools").mkdir(parents=True)
    (tmp_path / "wdl" / "annotation").mkdir()
    (tmp_path / "wdl" / "annotation_utils").mkdir()
    monkeypatch.setattr(sync, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(sync, "DOCKSTORE_YML", tmp_path / ".dockstore.yml")

    def build(workflows, yml_entries):
        for name in workflows:
            (tmp_path / "wdl" / "tools" / f"{name}.wdl").write_text("version 1.0\n")
        (tmp_path / ".dockstore.yml").write_text("version: 1.2\nworkflows:\n" + "".join(yml_entries))

    return build


def run(monkeypatch, capsys, *args):
    monkeypatch.setattr(sys, "argv", ["check_dockstore_sync.py", *args])
    try:
        sync.main()
        code = 0
    except SystemExit as exc:
        code = exc.code
    return code, capsys.readouterr().out


def test_check_branch_filters_requires_branches_list():
    assert sync.check_branch_filters("A", {}, False)
    assert sync.check_branch_filters("A", {"filters": {"branches": "main"}}, False)


def test_check_branch_filters_requires_main():
    errors = sync.check_branch_filters("A", {"filters": {"branches": ["kj-topic"]}}, False)
    assert len(errors) == 1 and "'main'" in errors[0]


def test_check_branch_filters_allows_feature_branch_unless_main_only():
    branch_entry = {"filters": {"branches": ["main", "kj-topic"]}}
    assert sync.check_branch_filters("A", branch_entry, False) == []
    errors = sync.check_branch_filters("A", branch_entry, True)
    assert len(errors) == 1 and "kj-topic" in errors[0]


def test_check_branch_filters_main_only_accepts_main_alone():
    assert sync.check_branch_filters("A", {"filters": {"branches": ["main"]}}, True) == []


def test_in_sync_repo_passes(fake_repo, monkeypatch, capsys):
    fake_repo(["A", "B"], [entry("A"), entry("B")])
    code, out = run(monkeypatch, capsys)
    assert code == 0
    assert "2 active workflow files" in out


def test_workflow_without_entry_fails(fake_repo, monkeypatch, capsys):
    fake_repo(["A", "B"], [entry("A")])
    code, out = run(monkeypatch, capsys)
    assert code == 1
    assert "no .dockstore.yml entry found" in out and "B.wdl" in out


def test_entry_without_file_fails(fake_repo, monkeypatch, capsys):
    fake_repo(["A"], [entry("A"), entry("Ghost")])
    code, out = run(monkeypatch, capsys)
    assert code == 1
    assert "missing file" in out and "no matching active workflow file" in out


def test_wrong_descriptor_path_fails(fake_repo, monkeypatch, capsys):
    fake_repo(["A"], [entry("A", path="/wdl/annotation/A.wdl")])
    code, out = run(monkeypatch, capsys)
    assert code == 1


def test_duplicate_entry_fails(fake_repo, monkeypatch, capsys):
    fake_repo(["A"], [entry("A"), entry("A")])
    code, out = run(monkeypatch, capsys)
    assert code == 1 and "duplicate" in out


def test_feature_branch_passes_by_default_but_fails_with_main_only(fake_repo, monkeypatch, capsys):
    fake_repo(["A"], [entry("A", branches=("main", "kj-topic"))])
    assert run(monkeypatch, capsys)[0] == 0
    code, out = run(monkeypatch, capsys, "--main-only")
    assert code == 1 and "kj-topic" in out


def test_repository_dockstore_is_in_sync(repo_root):
    entries = sync.load_dockstore_entries()
    names = {e["name"] for e in entries}
    assert names == {p.stem for p in sync.find_active_workflow_files()}
    assert all("main" in e["filters"]["branches"] for e in entries)
