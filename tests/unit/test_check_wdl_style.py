import pytest
from conftest import REPO_ROOT, load_script

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


PROVENANCE = "# Derived from wdl/Path.wdl of owner/repo (branch dev)\n# https://github.com/owner/repo/blob/main/wdl/Path.wdl\n"


def test_provenance_header_before_version_is_allowed():
    assert codes(PROVENANCE + "\n" + GOOD_TASK) == []
    assert codes(PROVENANCE * 2 + "\n" + GOOD_TASK) == []


def test_single_provenance_line_above_a_library_task_is_allowed():
    library = GOOD_TASK.replace("task AlphaTask", "# Derived from wdl/Path.wdl of owner/repo\ntask AlphaTask")
    assert codes(library) == []


@pytest.mark.parametrize(
    "header",
    [
        "# Derived from owner/repo (branch dev) wdl/Path.wdl\n# https://github.com/owner/repo/blob/main/wdl/Path.wdl\n",
        "# Derived from wdl/Path.wdl of owner/repo\n",
        PROVENANCE + "# Fork of other/repo\n",
    ],
)
def test_malformed_provenance_header_is_reported(header):
    assert "W035" in codes(header + "\n" + GOOD_TASK)


def test_script_provenance_lines_follow_the_standard_form():
    bad = []
    for path in sorted((REPO_ROOT / "scripts").rglob("*")):
        if not path.is_file() or path.suffix not in {".py", ".R", ".sh"}:
            continue
        lines = path.read_text().splitlines()
        for index, line in enumerate(lines):
            if style.PROVENANCE_PREFIX_RE.match(line) and not (
                style.PROVENANCE_RE.match(line) and index + 1 < len(lines) and style.PROVENANCE_URL_RE.match(lines[index + 1])
            ):
                bad.append(f"{path.relative_to(REPO_ROOT)}:{index + 1}")
    assert bad == []


def test_multi_line_comment_inside_file_is_reported():
    text = GOOD_TASK.replace("    output {", "    # first\n    # second\n    output {")
    assert "W027" in codes(text)


RUNTIME_TASK = GOOD_TASK.replace(
    "        docker: docker\n",
    "        docker: docker\n        memory: MEMORY\n        preemptible: PREEMPT\n",
)


def default_attr(preemptible="1", retries="0"):
    text = (
        "    RuntimeAttr default_attr = object {\n"
        "        mem_gb: 1,\n"
        "        preemptible_tries: PREEMPT,\n"
        "        max_retries: RETRIES\n"
        "    }\n"
    )
    return text.replace("PREEMPT", preemptible).replace("RETRIES", retries)


@pytest.mark.parametrize("value", ["1", "preemptible_default"])
def test_default_attr_allowed_preemptible_values_pass(value):
    assert "W030" not in codes(GOOD_TASK.replace("    runtime {", default_attr(preemptible=value) + "    runtime {"))


@pytest.mark.parametrize("preemptible,retries", [("0", "0"), ("3", "0"), ("1", "1"), ("4", "2")])
def test_default_attr_literal_retry_values_are_reported(preemptible, retries):
    text = GOOD_TASK.replace("    runtime {", default_attr(preemptible, retries) + "    runtime {")
    assert "W030" in codes(text)


def test_default_attr_identifier_retries_pass():
    text = GOOD_TASK.replace("    runtime {", default_attr(retries="default_retries") + "    runtime {")
    assert "W030" not in codes(text)


def test_retry_literals_inside_command_are_ignored():
    assert "W030" not in codes(GOOD_TASK.replace("echo hi", "echo max_retries: 5"))


def test_standard_memory_expression_passes():
    standard = 'select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"'
    assert "W031" not in codes(RUNTIME_TASK.replace("MEMORY", standard).replace("PREEMPT", "1"))


@pytest.mark.parametrize(
    "memory",
    [
        '"1 GiB"',
        '4 + " GiB"',
        'select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"',
        'select_first([default_attr.mem_gb, runtime_attr.mem_gb]) + " GiB"',
    ],
)
def test_nonstandard_memory_is_reported(memory):
    assert "W031" in codes(RUNTIME_TASK.replace("MEMORY", memory).replace("PREEMPT", "1"))


def test_allowlisted_memory_task_passes():
    text = RUNTIME_TASK.replace("MEMORY", '"10 GiB"').replace("PREEMPT", "1").replace("AlphaTask", "GermlineCNVCallerCaseMode")
    assert "W031" not in codes(text, "wdl/utils/LRCNVs.wdl")
    assert "W031" in codes(text, "wdl/utils/Other.wdl")


def test_task_names_must_be_unique_across_files():
    sources = [("wdl/a/One.wdl", GOOD_TASK), ("wdl/a/Two.wdl", GOOD_TASK.replace("AlphaTask", "BetaTask"))]
    assert style.duplicate_task_errors(sources) == []
    duplicate = style.duplicate_task_errors(sources + [("wdl/a/Three.wdl", GOOD_TASK)])
    assert [error[:3] for error in duplicate] == [("wdl/a/Three.wdl", 3, "W032")]


def test_task_name_in_command_text_is_not_a_declaration():
    text = GOOD_TASK.replace("AlphaTask", "BetaTask").replace("echo hi", "\ntask AlphaTask")
    assert style.duplicate_task_errors([("wdl/a/One.wdl", GOOD_TASK), ("wdl/a/Two.wdl", text)]) == []


@pytest.mark.parametrize(
    "name",
    [
        "PALMER_vcf",
        "chrY_coverage_cutoff",
        "bamOut",
        "ref_fasta",
        "ref_fasta_index",
        "ref_fa_idx",
        "vcf_index",
        "vcf_tbi",
        "combined_tbi",
    ],
)
def test_bad_input_names_are_reported(name):
    assert "W033" in codes(GOOD_TASK.replace("String docker", "String " + name))


@pytest.mark.parametrize("name", ["ref_fa", "ref_fai", "vcf_idx", "more_pl", "chry_coverage_cutoff", "bam", "base_vcf_number"])
def test_good_input_names_pass(name):
    assert "W033" not in codes(GOOD_TASK.replace("String docker", "String " + name))


def test_bad_output_name_is_reported():
    assert "W033" in codes(GOOD_TASK.replace("String message", "String statusOut"))


def test_struct_members_are_exempt_from_naming():
    struct = "version 1.0\n\nstruct Bundle {\n    File fasta\n    File PAR_bed\n}\n"
    assert "W033" not in codes(struct)


def test_workflow_docker_input_must_end_in_docker():
    workflow = "wdl/annotation_utils/Example.wdl"
    assert "W033" in codes(GOOD_WORKFLOW.replace("String utils_docker", "String docker"), workflow)
    assert "W033" in codes(GOOD_WORKFLOW.replace("String utils_docker", "String dockerimage"), workflow)
    assert "W033" not in codes(GOOD_WORKFLOW, workflow)


def test_task_docker_input_is_not_renamed():
    assert "W033" not in codes(GOOD_TASK)


def test_repository_wdl_passes_style_check(repo_root):
    errors = []
    sources = []
    for path in style.collect_files(["wdl"]):
        rel, text = str(path.relative_to(repo_root)), path.read_text()
        sources.append((rel, text))
        errors.extend(style.Checker(rel, text).run())
    errors.extend(style.duplicate_task_errors(sources))
    assert errors == []


NESTED_WORKFLOW = """\
version 1.0

workflow Nested {
    input {
        Array[String] contigs
        Int records = 10
        RuntimeAttr? runtime_attr_inner
    }

    Array[RuntimeAttr?] runtime_attr_inner_wrapped = [runtime_attr_inner]

    scatter (contig in contigs) {
        if (records > 0) {
            call Shard { input: n = records }
        }
        INNER
    }
}
"""


def nested_hits(inner):
    text = NESTED_WORKFLOW.replace("INNER", inner)
    return [name for _, name in style.Checker("wdl/annotation/Nested.wdl", text).nested_scatter_optionals()]


@pytest.mark.parametrize(
    "inner,expected",
    [
        ("scatter (i in range(2)) {\n call Inner { input: attr = runtime_attr_inner }\n }", ["runtime_attr_inner"]),
        ('scatter (i in range(2)) {\n call Inner { input: s = "~{runtime_attr_inner}" }\n }', ["runtime_attr_inner"]),
        ("scatter (i in range(2)) {\n call Inner { input: s = Shard.out }\n }", ["Shard"]),
        ("scatter (i in range(2)) {\n call Inner { input: attr = runtime_attr_inner_wrapped[0] }\n }", []),
        ("scatter (i in range(2)) {\n call Inner { input: n = records, s = contig, i = i }\n }", []),
        ("call Inner { input: attr = runtime_attr_inner }", []),
    ],
)
def test_undefined_optionals_entering_nested_scatters_are_reported(inner, expected):
    assert nested_hits(inner) == expected


def test_conditional_call_used_inside_its_own_if_passes():
    inner = (
        "if (records > 0) {\n call Shard2 { input: n = 1 }\n"
        " scatter (i in Shard2.out) {\n call Inner { input: s = Shard2.out }\n }\n }"
    )
    assert nested_hits(inner) == []
