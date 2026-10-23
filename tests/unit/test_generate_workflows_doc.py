import re

import pytest
from conftest import load_script

gen = load_script(".github/scripts/generate_workflows_doc.py")
wdl_meta = gen.wdl_meta


def test_built_document_contains_every_active_workflow():
    text = gen.build(gen.ACTIVE_SITE)
    names = [wdl_meta.parse_workflow(path).name for path in wdl_meta.find_workflow_files()]
    assert names
    for name in names:
        assert f"### [{name}](" in text, name


def test_built_document_has_expected_sections_in_order():
    text = gen.build(gen.ACTIVE_SITE)
    headings = [line for line in text.splitlines() if line.startswith("## ")]
    assert headings == ["## Contents", "## Annotations", "## Annotation Utilities", "## Tools", "## Sub-workflows"]


def test_dockstore_order_groups_match_sections():
    groups = gen.dockstore_order()
    for _, group, _ in gen.SECTIONS:
        assert groups[group]


def test_render_default_unwraps_quoted_strings():
    assert gen.render_default('"abc"') == " (default `abc`)"
    assert gen.render_default('""') == " (default empty)"
    assert gen.render_default("5") == " (default `5`)"


def test_render_workflow_lists_documented_inputs_outputs_and_standard_bullets(repo_root):
    path = repo_root / "wdl" / "annotation_utils" / "StripGenotypes.wdl"
    workflow = wdl_meta.parse_workflow(path)
    lines = gen.render_workflow(workflow, path, gen.ACTIVE_SITE)
    assert lines[0] == "### [StripGenotypes](../wdl/annotation_utils/StripGenotypes.wdl)"
    assert "Inputs:" in lines and "Outputs:" in lines
    assert any(line.startswith("- `String prefix`: Prefix for output file names.") for line in lines)
    assert any("Optional per-task runtime overrides" in line for line in lines)
    for decl in workflow.outputs:
        assert any(line.startswith(f"- `{decl.type} {decl.name}`:") for line in lines)


def github_anchors(text):
    """Anchors GitHub gives the headings of a rendered document, in order, with duplicates suffixed."""
    seen = {}
    anchors = set()
    in_fence = False
    for line in text.splitlines():
        if line.startswith("```"):
            in_fence = not in_fence
        match = re.match(r"^#{1,6}\s+(.*?)\s*$", line)
        if in_fence or not match:
            continue
        heading = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", match.group(1))
        slug = re.sub(r"[^\w\- ]", "", heading.lower()).replace(" ", "-")
        count = seen.get(slug, 0)
        seen[slug] = count + 1
        anchors.add(slug if count == 0 else f"{slug}-{count}")
    return anchors


@pytest.mark.parametrize("site", [gen.ACTIVE_SITE, gen.ARCHIVE_SITE], ids=["active", "archive"])
def test_contents_links_every_workflow_and_resolves_to_headings(site):
    text = gen.build(site)
    lines = text.splitlines()
    start = lines.index("## Contents")
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("## "))
    links = re.findall(r"^\s*- \[([^\]]+)\]\(#([^)]+)\)$", "\n".join(lines[start + 1 : end]), re.M)
    assert links
    anchors = github_anchors(text)
    for _, anchor in links:
        assert anchor in anchors, anchor
    headings = {m.group(1) for m in re.finditer(r"^#{2,3} (?:\[([^\]]+)\]\([^)]*\)|(.+))$", text, re.M) if m.group(1)}
    assert headings <= {name for name, _ in links}


def test_anchor_slugs_follow_github_rules():
    assert gen.anchor_slugs(["Depth-based CNV pipeline", "Sub-workflows", "Name_1", "Name_1"]) == [
        "depth-based-cnv-pipeline",
        "sub-workflows",
        "name_1",
        "name_1-1",
    ]
