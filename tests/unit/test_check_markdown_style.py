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
