import pytest
from conftest import load_script

agents = load_script(".github/scripts/check_agents_sync.py")


def test_repository_mirrors_match_agents_md(capsys):
    agents.main()
    assert "passed" in capsys.readouterr().out


def test_drifted_mirror_fails(tmp_path, monkeypatch, capsys):
    (tmp_path / "AGENTS.md").write_text("canonical\n")
    good, drifted = tmp_path / "good.md", tmp_path / "drifted.md"
    good.write_text("canonical\n")
    drifted.write_text("changed\n")
    monkeypatch.setattr(agents, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(agents, "CANONICAL", tmp_path / "AGENTS.md")
    monkeypatch.setattr(agents, "MIRRORS", [good, drifted, tmp_path / "missing.md"])
    with pytest.raises(SystemExit) as exc:
        agents.main()
    assert exc.value.code == 1
    out = capsys.readouterr().out
    assert "drifted.md content has drifted" in out
    assert "missing.md is missing" in out
    assert "good.md" not in out
