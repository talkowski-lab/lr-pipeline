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
    assert headings == ["## Annotations", "## Annotation Utilities", "## Tools", "## Sub-workflows"]


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
