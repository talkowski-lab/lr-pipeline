import subprocess

import pytest
from conftest import load_script

md = load_script(".github/scripts/check_markdown_style.py")

GOOD = """\
# Title
Intro.


## Section
- one
  - nested
- two


### First
Text.

### Second
Text.
"""


def codes(text):
    return [error[2] for error in md.Checker("x.md", text).run()]


def test_clean_document_passes():
    assert codes(GOOD) == []


@pytest.mark.parametrize(
    "text,code",
    [
        (GOOD.replace("- one", "-\tone"), "M001"),
        (GOOD.replace("Intro.", "Intro. "), "M002"),
        (GOOD.rstrip("\n"), "M003"),
        (GOOD + "\n", "M003"),
        (GOOD.replace("## Section\n", "## Section\n\n"), "M004"),
        (GOOD.replace("\n\n\n## Section", "\n\n\n\n## Section"), "M005"),
        (GOOD.replace("Intro.", "Intro.\n\n---"), "M006"),
        (GOOD.replace("  - nested", "   - nested"), "M007"),
    ],
)
def test_violations_are_reported(text, code):
    assert code in codes(text)


def test_fenced_code_is_exempt_from_heading_and_divider_rules():
    text = "# Title\n```\n## not a heading\n---\n```\n"
    assert codes(text) == []


def test_fenced_code_still_checks_trailing_whitespace():
    assert "M002" in codes("# Title\n```\ncode \n```\n")


def test_repository_docs_pass(repo_root):
    errors = []
    for path in md.collect_files(["README.md", "docs"]):
        errors.extend(md.Checker(str(path), path.read_text()).run())
    assert errors == []


@pytest.mark.parametrize(
    "texts,orphans",
    [
        # The root README needs no referrer; a file mentioning itself does not count
        ({"README.md": "", "docs/a.md": "see docs/a.md"}, ["docs/a.md"]),
        ({"README.md": "[A](docs/a.md#top)", "docs/a.md": ""}, []),
        # Relative links resolve from the mentioning file's directory
        ({"README.md": "", "docs/a.md": "[B](../tests/b.md)", "tests/b.md": "[A](../docs/a.md)"}, []),
        ({"README.md": "tests/README.md", "tests/README.md": "[data](data/README.md)", "tests/data/README.md": ""}, []),
        # A unique basename counts from anywhere, including non-Markdown files
        ({"README.md": "", "docs/a.md": "", "scripts/x.py": "# documented in a.md"}, []),
        # A shared basename alone does not identify the file
        ({"README.md": "see README.md", "tests/README.md": ""}, ["tests/README.md"]),
        (
            {"README.md": "workflows.md", "docs/workflows.md": "", "archive/docs/workflows.md": ""},
            ["archive/docs/workflows.md", "docs/workflows.md"],
        ),
        (
            {
                "README.md": "",
                "CONTRIBUTING.md": "`archive/docs/workflows.md`",
                "archive/docs/workflows.md": "docs/workflows.md",
                "docs/workflows.md": "CONTRIBUTING.md",
            },
            [],
        ),
    ],
)
def test_find_unreferenced(texts, orphans):
    assert md.find_unreferenced(texts) == orphans


def test_check_references_reads_tracked_files(tmp_path):
    (tmp_path / "docs").mkdir()
    (tmp_path / "README.md").write_text("[A](docs/a.md)\n")
    (tmp_path / "docs" / "a.md").write_text("")
    (tmp_path / "docs" / "orphan.md").write_text("")
    (tmp_path / "docs" / "untracked.md").write_text("")
    (tmp_path / "binary.dat").write_bytes(b"\xff\xfe orphan.md")
    (tmp_path / "link.md").symlink_to("README.md")
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    subprocess.run(["git", "add", "README.md", "docs/a.md", "docs/orphan.md", "binary.dat", "link.md"], cwd=tmp_path, check=True)
    assert [(error[0], error[2]) for error in md.check_references(tmp_path)] == [("docs/orphan.md", "M008")]


def test_repository_markdown_is_referenced(repo_root):
    assert md.check_references(repo_root) == []
