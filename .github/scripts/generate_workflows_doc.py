#!/usr/bin/env python3
"""Generate a workflow document from the meta and parameter_meta blocks of each workflow.

Every section is derived from the WDL: the description paragraphs come from
`meta.description`, the input and output bullets from `parameter_meta`, and each
bullet's type, optionality and default from the declaration itself. Nothing in the
document is hand-written except the constants in this file.

Two documents are generated. `docs/workflows.md` covers the active pipeline, ordering
its sections by `.dockstore.yml` so the document and the Dockstore registration cannot
drift apart; the sub-workflows under `wdl/utils` are not registered, so their order and
grouping headings live in SUBWORKFLOW_GROUPS. `archive/docs/workflows.md` covers the
retired workflows, which have no Dockstore entries and are listed alphabetically.

Run with no arguments to write the file. Run with --check to verify the committed copy
is current, or --dry-run to confirm the document can be built without writing it.

Usage: generate_workflows_doc.py [--site {active,archive}] [--check | --dry-run]
"""
import argparse
import difflib
import os
import re
import sys
from pathlib import Path

import wdl_meta

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCKSTORE_YML = REPO_ROOT / ".dockstore.yml"

TITLE = "Workflows"
INTRO = (
    "This document describes each WDL workflow in the pipeline, including its purpose, inputs and outputs. "
    "Annotations, annotation utilities and tools are run directly and are registered in `.dockstore.yml`; "
    "the sub-workflows in the final section are imported building blocks that are never run on their own."
)
GENERATED_NOTE = (
    "This file is generated from the `meta` and `parameter_meta` blocks of each workflow by "
    "[`generate_workflows_doc.py`](../.github/scripts/generate_workflows_doc.py). Edit those blocks rather than "
    "this document. Inputs described as `From references.` are the shared reference files listed in "
    "[references](references.md)."
)

# Registered sections, paired with the comment block that orders them in .dockstore.yml
SECTIONS = [
    ("Annotations", "Annotation Workflows", "wdl/annotation"),
    ("Annotation Utilities", "Annotation Utilities", "wdl/annotation_utils"),
    ("Tools", "Tools", "wdl/tools"),
]

SUBWORKFLOWS_TITLE = "Sub-workflows"
SUBWORKFLOWS_INTRO = (
    "These workflows live in `wdl/utils/` and are building blocks rather than entry points. They are never "
    'registered in `.dockstore.yml` and are not run directly; a workflow imports one with `import "../utils/<Name>.wdl"` '
    "and calls it as `<Name>.<Name>`."
)
SUBWORKFLOW_GROUPS = [
    (
        "Depth-based CNV pipeline",
        "`LRCNVs`, `DepthPreprocessing`, `DepthClustering` and `GenotypeDepth` are called in that order by "
        "`LongReadCNVs`, which supplies their shared inputs.",
        ["LRCNVs", "DepthPreprocessing", "DepthClustering", "GenotypeDepth"],
    ),
    (
        "Callset matching and sharding",
        "`ExactMatch`, `TruvariMatch` and `BedtoolsClosestSV` are the three comparison rounds driven by "
        "`AnnotateCallsetOverlap`; each consumes what the previous round left unmatched. `ScatterVcf` is a "
        "general sharding helper.",
        ["ExactMatch", "TruvariMatch", "BedtoolsClosestSV", "ScatterVcf"],
    ),
]

ARCHIVE_TITLE = "Long-Read Annotation"
ARCHIVE_INTRO = (
    "This document describes each retired WDL workflow, including its purpose, inputs and outputs. Archived "
    "workflows are retained as historical reference only. They are not active pipeline entry points, and they "
    "are excluded from Dockstore registration and from the validation and style checks that run over `wdl/`."
)
ARCHIVE_NOTE = (
    "This file is generated from the `meta` and `parameter_meta` blocks of each workflow by "
    "[`generate_workflows_doc.py`](../../.github/scripts/generate_workflows_doc.py). Edit those blocks rather "
    "than this document. Inputs described as `From references.` are the shared reference files listed in "
    "[references](../../docs/references.md)."
)
# Archived workflows have no .dockstore.yml entries, so each section is listed alphabetically
ARCHIVE_SECTIONS = [
    ("Annotations", "archive/wdl/annotation"),
    ("Annotation Utilities", "archive/wdl/annotation_utils"),
    ("Tools", "archive/wdl/tools"),
]


class Site:
    """One generated document: where it is written and how its sections are ordered."""

    def __init__(self, title, intro, note, output, directories):
        self.title = title
        self.intro = intro
        self.note = note
        self.output = REPO_ROOT / output
        self.directories = directories


ACTIVE_SITE = Site(TITLE, INTRO, GENERATED_NOTE, "docs/workflows.md", wdl_meta.WORKFLOW_DIRS)
ARCHIVE_SITE = Site(
    ARCHIVE_TITLE,
    ARCHIVE_INTRO,
    ARCHIVE_NOTE,
    "archive/docs/workflows.md",
    [directory for _, directory in ARCHIVE_SECTIONS],
)
SITES = {"active": ACTIVE_SITE, "archive": ARCHIVE_SITE}

DOCKSTORE_GROUP_RE = re.compile(r"^#\s+(\w[\w ]*?)\s*$")
DOCKSTORE_NAME_RE = re.compile(r"^\s*-?\s*name:\s*(\w+)\s*$")
QUOTED_RE = re.compile(r"""^(?P<quote>["'])(?P<value>.*)(?P=quote)$""")


def dockstore_order():
    """Map each .dockstore.yml comment block to the workflow names listed under it."""
    groups = {}
    current = None
    for line in DOCKSTORE_YML.read_text().split("\n"):
        stripped = line.strip()
        match = DOCKSTORE_GROUP_RE.match(stripped)
        if match:
            current = match.group(1)
            groups.setdefault(current, [])
            continue
        name = DOCKSTORE_NAME_RE.match(line)
        if name and current is not None:
            groups[current].append(name.group(1))
    return groups


def render_default(default):
    """Render an input default for a bullet, unwrapping a quoted string."""
    quoted = QUOTED_RE.match(default.strip())
    if not quoted:
        return f" (default `{default.strip()}`)"
    value = quoted.group("value")
    if not value:
        return " (default empty)"
    return f" (default `{value}`)"


def render_bullet(decl, description):
    default = render_default(decl.default) if decl.default is not None else ""
    return f"- `{decl.type} {decl.name}`: {description}{default}"


def exempt_bullets(workflow):
    """Render the standard bullets for the inputs authors never describe."""
    bullets = []
    prefix = [d for d in workflow.inputs if d.type == "String" and d.name == wdl_meta.PREFIX_INPUT]
    dockers = [d for d in workflow.inputs if d.type == "String" and wdl_meta.DOCKER_INPUT_RE.match(d.name)]
    runtime = [d for d in workflow.inputs if d.type == wdl_meta.RUNTIME_ATTR_TYPE]
    for decl in prefix:
        bullets.append(render_bullet(decl, "Prefix for output file names."))
    if dockers:
        names = ", ".join(f"`{d.type} {d.name}`" for d in dockers)
        label = "Container image." if len(dockers) == 1 else "Container images."
        bullets.append(f"- {names}: {label}")
    if runtime:
        count = f" ({len(runtime)})" if len(runtime) > 1 else ""
        bullets.append(f"- `{wdl_meta.RUNTIME_ATTR_TYPE} runtime_attr_*`: Optional per-task runtime overrides{count}.")
    return bullets


def render_workflow(workflow, path, site):
    """Render one `###` section for a workflow, linking to it from the document's directory."""
    relative = os.path.relpath(str(path), str(site.output.parent))
    lines = [f"### [{workflow.name}]({relative})"]
    for index, paragraph in enumerate(workflow.description):
        if index:
            lines.append("")
        lines.append(paragraph)

    described = dict((key, value) for key, value, _ in workflow.param_meta)
    inputs = []
    for decl in workflow.inputs:
        if wdl_meta.is_exempt(decl):
            continue
        inputs.append(render_bullet(decl, described[decl.name]))
    inputs.extend(exempt_bullets(workflow))
    lines.append("")
    lines.append("Inputs:")
    lines.extend(inputs)

    if workflow.outputs:
        lines.append("")
        lines.append("Outputs:")
        for decl in workflow.outputs:
            lines.append(f"- `{decl.type} {decl.name}`: {described[decl.name]}")
    return lines


def load_workflows(directories):
    """Parse every workflow, keyed by name, failing on anything the style checker would reject."""
    workflows = {}
    problems = []
    for path in wdl_meta.find_workflow_files(directories):
        workflow = wdl_meta.parse_workflow(path)
        if workflow is None:
            continue
        described = set(key for key, _, _ in workflow.param_meta)
        for decl in workflow.inputs:
            if not wdl_meta.is_exempt(decl) and decl.name not in described:
                problems.append(f"{path.name}: input '{decl.name}' has no parameter_meta entry")
        for decl in workflow.outputs:
            if decl.name not in described:
                problems.append(f"{path.name}: output '{decl.name}' has no parameter_meta entry")
        if not workflow.description:
            problems.append(f"{path.name}: workflow '{workflow.name}' has no meta description")
        for line_no, code, message in workflow.errors:
            problems.append(f"{path.name}:{line_no}: {code} {message}")
        workflows[workflow.name] = (workflow, path)
    return workflows, problems


def build(site):
    """Render the whole document for a site, returning its text."""
    workflows, problems = load_workflows(site.directories)
    if problems:
        print("generate_workflows_doc.py: cannot build the document:\n")
        for problem in sorted(problems):
            print(f"  - {problem}")
        print(f"\n{len(problems)} problem(s). Run check_wdl_style.py for the full report.")
        sys.exit(2)

    lines = [f"# {site.title}", site.intro, "", site.note]
    used = set()
    if site is ARCHIVE_SITE:
        for title, directory in ARCHIVE_SECTIONS:
            names = sorted(
                name for name, (_, path) in workflows.items()
                if Path(path).parent.as_posix() == (REPO_ROOT / directory).as_posix()
            )
            if not names:
                continue
            lines.extend(["", "", f"## {title}"])
            for index, name in enumerate(names):
                workflow, path = workflows[name]
                lines.append("")
                if index == 0:
                    lines.append("")
                lines.extend(render_workflow(workflow, path, site))
                used.add(name)
        missing = sorted(set(workflows) - used)
        if missing:
            sys.exit("generate_workflows_doc.py: not placed in any section: {}".format(", ".join(missing)))
        return "\n".join(lines) + "\n"

    groups = dockstore_order()
    for title, group, directory in SECTIONS:
        names = groups.get(group)
        if not names:
            sys.exit(f"generate_workflows_doc.py: no '{group}' block in .dockstore.yml")
        lines.extend(["", "", f"## {title}"])
        for index, name in enumerate(names):
            workflow, path = workflows[name]
            if Path(path).parent.as_posix() != (REPO_ROOT / directory).as_posix():
                sys.exit(f"generate_workflows_doc.py: '{name}' is registered under '{group}' but lives elsewhere")
            lines.append("")
            if index == 0:
                lines.append("")
            lines.extend(render_workflow(workflow, path, site))
            used.add(name)

    lines.extend(["", "", f"## {SUBWORKFLOWS_TITLE}", SUBWORKFLOWS_INTRO])
    for group_index, (heading, intro, names) in enumerate(SUBWORKFLOW_GROUPS):
        lines.append("")
        if group_index == 0:
            lines.append("")
        lines.extend([f"### {heading}", intro])
        for name in names:
            workflow, path = workflows[name]
            lines.append("")
            lines.extend(render_workflow(workflow, path, site))
            used.add(name)

    missing = sorted(set(workflows) - used)
    if missing:
        sys.exit(
            "generate_workflows_doc.py: not placed in any section: {}. Register it in .dockstore.yml "
            "or add it to SUBWORKFLOW_GROUPS.".format(", ".join(missing))
        )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--site", choices=sorted(SITES), default="active", help="which document to generate")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="fail if the document is not up to date")
    mode.add_argument("--dry-run", action="store_true", help="build the document without writing it")
    args = parser.parse_args()

    site = SITES[args.site]
    text = build(site)
    output = site.output
    relative = output.relative_to(REPO_ROOT).as_posix()
    if args.dry_run:
        print(f"generate_workflows_doc.py: {relative} builds cleanly ({len(text.splitlines())} lines).")
        return
    if args.check:
        current = output.read_text() if output.exists() else ""
        if current == text:
            print(f"generate_workflows_doc.py: {relative} is up to date.")
            return
        diff = difflib.unified_diff(
            current.split("\n"),
            text.split("\n"),
            fromfile=f"{relative} (committed)",
            tofile=f"{relative} (generated)",
            lineterm="",
        )
        print(f"{relative} is out of date:\n")
        print("\n".join(diff))
        print(f"\nRun 'python .github/scripts/generate_workflows_doc.py --site {args.site}' and commit the result.")
        sys.exit(1)

    output.write_text(text)
    print(f"generate_workflows_doc.py: wrote {relative} ({len(text.splitlines())} lines).")


if __name__ == "__main__":
    main()
