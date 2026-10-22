from conftest import load_script

wdl_meta = load_script(".github/scripts/wdl_meta.py")

WORKFLOW = """\
version 1.0

task Ignored {
    input {
        File not_a_workflow_input
    }
}

workflow Example {
    meta {
        description: [
            "First paragraph.",
            "Second paragraph."
        ]
    }

    parameter_meta {
        vcf: "Input VCF."
        sample_count: "Number of samples."
        out: "Output VCF."
    }

    input {
        File vcf
        Int sample_count = 5
        String prefix
        String utils_docker
        RuntimeAttr? runtime_attr_x
    }

    output {
        File out = vcf
    }
}
"""


def test_parse_workflow_reads_name_and_description():
    workflow = wdl_meta.parse_workflow("Example.wdl", WORKFLOW)
    assert workflow.name == "Example"
    assert workflow.description == ["First paragraph.", "Second paragraph."]
    assert workflow.description_is_array
    assert workflow.errors == []


def test_parse_workflow_reads_declarations_and_skips_task_bodies():
    workflow = wdl_meta.parse_workflow("Example.wdl", WORKFLOW)
    assert [d.name for d in workflow.inputs] == ["vcf", "sample_count", "prefix", "utils_docker", "runtime_attr_x"]
    assert [d.name for d in workflow.outputs] == ["out"]
    assert workflow.declaration("sample_count").default == "5"
    assert workflow.declaration("runtime_attr_x").optional
    assert workflow.declaration("missing") is None


def test_parse_workflow_reads_parameter_meta_in_order():
    workflow = wdl_meta.parse_workflow("Example.wdl", WORKFLOW)
    assert workflow.param_meta_keys == ["vcf", "sample_count", "out"]
    assert workflow.param_meta[0][1] == "Input VCF."
    assert [name for name, _ in workflow.blocks] == ["meta", "parameter_meta", "input", "output"]


def test_is_exempt_covers_prefix_docker_and_runtime_attr():
    workflow = wdl_meta.parse_workflow("Example.wdl", WORKFLOW)
    exempt = [d.name for d in workflow.inputs if wdl_meta.is_exempt(d)]
    assert exempt == ["prefix", "utils_docker", "runtime_attr_x"]


def test_parse_workflow_returns_none_without_workflow():
    assert wdl_meta.parse_workflow("Lib.wdl", "version 1.0\n\ntask A {\n}\n") is None


def test_unparseable_input_line_is_reported_not_dropped():
    workflow = wdl_meta.parse_workflow("Example.wdl", WORKFLOW.replace("        File vcf\n", "        File\n"))
    assert [code for _, code, _ in workflow.errors] == ["W025"]


def test_unescape_restores_quotes_and_backslashes():
    assert wdl_meta.unescape('say \\"hi\\" \\\\') == 'say "hi" \\'


def test_find_workflow_files_only_returns_files_declaring_a_workflow(repo_root):
    files = wdl_meta.find_workflow_files(["wdl/annotation_utils", "wdl/utils"])
    assert files
    assert all("workflow " in path.read_text() for path in files)
    assert repo_root / "wdl/utils/Helpers.wdl" not in files
