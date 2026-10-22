import pytest
from conftest import load_script

style = load_script(".github/scripts/check_wdl_style.py")

GOOD_TASK = """\
version 1.0

task AlphaTask {
    input {
        String docker
    }

    command <<<
        set -euo pipefail

        echo hi
    >>>

    output {
        String message = "hi"
    }

    runtime {
        docker: docker
    }
}
"""

GOOD_WORKFLOW = """\
version 1.0

workflow Example {
    meta {
        description: [
            "Does one thing."
        ]
    }

    parameter_meta {
        vcf: "Input VCF."
        out: "Output file."
    }

    input {
        File vcf
        String prefix
        String utils_docker
    }

    output {
        File out = vcf
    }
}
"""


def codes(text, path="wdl/utils/Snippet.wdl"):
    return [error[2] for error in style.Checker(path, text).run()]


def test_clean_task_passes():
    assert codes(GOOD_TASK) == []


def test_clean_workflow_with_full_documentation_passes():
    assert codes(GOOD_WORKFLOW, "wdl/annotation_utils/Example.wdl") == []


@pytest.mark.parametrize(
    "mutate,code",
    [
        (lambda t: t.replace("    input {", "\tinput {"), "W001"),
        (lambda t: t.replace("echo hi", "echo hi  "), "W002"),
        (lambda t: t.replace("\n    output", "\n\n    output"), "W003"),
        (lambda t: t.rstrip("\n"), "W004"),
        (lambda t: t.replace("        set -euo", "         set -euo"), "W006"),
        (lambda t: t.replace("set -euo pipefail", "set -euxo pipefail"), "W007"),
        (lambda t: t.replace("version 1.0", "version 1.1"), "W011"),
        (lambda t: t.replace("AlphaTask", "alpha_task"), "W012"),
        (lambda t: t.replace("command <<<", "command {").replace(">>>", "}"), "W006"),
        (lambda t: t.replace("version 1.0\n", 'version 1.0\n\nimport "A.wdl" as B\n'), "W010"),
    ],
)
def test_task_violations_are_reported(mutate, code):
    assert code in codes(mutate(GOOD_TASK))


def test_task_blocks_out_of_order_are_reported():
    text = GOOD_TASK.replace('    output {\n        String message = "hi"\n    }\n\n', "")
    text = text.replace("    command <<<", "    runtime {\n        docker: docker\n    }\n\n    command <<<", 1)
    assert "W008" in codes(text)


def test_missing_runtime_block_is_reported():
    text = GOOD_TASK.replace("\n    runtime {\n        docker: docker\n    }\n", "")
    assert "W008" in codes(text)


def test_task_library_must_declare_tasks_alphabetically():
    second = GOOD_TASK.split("\n", 2)[2].replace("AlphaTask", "BetaTask")
    assert "W017" not in codes(GOOD_TASK + "\n" + second)
    assert "W017" in codes(second + "\n" + GOOD_TASK.split("\n", 2)[2])


def test_workflow_name_must_match_file_stem():
    assert "W013" in codes(GOOD_WORKFLOW, "wdl/annotation_utils/Other.wdl")


def test_undocumented_input_is_reported():
    text = GOOD_WORKFLOW.replace('        vcf: "Input VCF."\n', "")
    assert "W022" in codes(text, "wdl/annotation_utils/Example.wdl")


def test_markdown_in_description_is_reported():
    text = GOOD_WORKFLOW.replace("Does one thing.", "Does **one** thing.")
    assert "W024" in codes(text, "wdl/annotation_utils/Example.wdl")


def test_unused_runtime_attr_is_reported():
    text = GOOD_WORKFLOW.replace(
        "        String utils_docker\n", "        String utils_docker\n        RuntimeAttr? runtime_attr_x\n"
    )
    assert "W016" in codes(text, "wdl/annotation_utils/Example.wdl")


def test_provenance_header_before_version_is_allowed():
    assert codes("# Derived from upstream/repo path.wdl\n\n" + GOOD_TASK) == []


def test_multi_line_comment_inside_file_is_reported():
    text = GOOD_TASK.replace("    output {", "    # first\n    # second\n    output {")
    assert "W027" in codes(text)


def test_repository_wdl_passes_style_check(repo_root):
    errors = []
    for path in style.collect_files(["wdl"]):
        errors.extend(style.Checker(str(path.relative_to(repo_root)), path.read_text()).run())
    assert errors == []
