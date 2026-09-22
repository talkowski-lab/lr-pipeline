#!/usr/bin/env python3
"""Move the hand-written content of docs/workflows.md into each workflow's WDL file.

This is the one-time migration that made docs/workflows.md generated. It reads a
snapshot of the hand-written document, converts each section into a `meta` block and a
`parameter_meta` block, and writes those blocks into the matching WDL file. It is kept
in the repository as the record of where the generated content came from, and it is
safe to re-run: existing `meta` and `parameter_meta` blocks are replaced, not appended.

Descriptions are rewritten to satisfy the authoring rules in docs/conventions.md:
Markdown links become plain text or a bare tool URL, emphasis is dropped, block quotes
and numbered lists are folded into paragraphs, and double quotes become single quotes.

A parameter with no bullet in the snapshot gets a mechanical stub where its name
implies one - an index file, or a pass-through of a GATK gCNV command-line flag - and
`TODO.` otherwise, which the style checker then surfaces for an author to replace.

Usage: migrate_workflows_doc.py SNAPSHOT [--report]
"""
import argparse
import re
from pathlib import Path

import wdl_meta

REPO_ROOT = Path(__file__).resolve().parents[2]

SECTION_RE = re.compile(r"^### \[(\w+)\]\(([^)]+)\)\s*$")
GROUP_HEADING_RE = re.compile(r"^### (?!\[)")
HEADING_RE = re.compile(r"^#{1,2} ")
BULLET_RE = re.compile(r"^- (.*)$")
TOKEN_RE = re.compile(r"^`([^`]+)`(?:, `([^`]+)`)*")
DEFAULT_RE = re.compile(r"\s*\(default [^)]*\)")
LINK_RE = re.compile(r"\[([^\]]*)\]\(([^)]+)\)")
# Underscore emphasis needs a non-word boundary on both sides, so snake_case names survive
EMPHASIS_RE = re.compile(r"\*\*([^*]+)\*\*|\*([^*]+)\*|(?<![\w`])_(?![\s_])([^_]+?)_(?![\w`])")
NOTE_RE = re.compile(r"^>\s*(?:\*\*)?Note:?(?:\*\*)?:?\s*", re.IGNORECASE)
LIST_ITEM_RE = re.compile(r"^(\d+)\.\s+")
BRACE_RE = re.compile(r"^(\w*)\{([\w,]+)\}(\w*)$")
INDEX_SUFFIXES = ["_idx", "_fai", "_bai", "_tbi", "_index"]

# Sub-workflow sections are grouped under a plain `###` heading rather than one per workflow
SKIP_HEADINGS = {"Depth-based CNV pipeline", "Callset matching and sharding"}


class Section:
    """One `###` workflow section of the hand-written document."""

    def __init__(self, name, link):
        self.name = name
        self.link = link
        self.prose = []
        self.inputs = []
        self.outputs = []


def parse_snapshot(text):
    """Split the hand-written document into workflow sections."""
    sections = []
    current = None
    bucket = None
    for line in text.split("\n"):
        match = SECTION_RE.match(line)
        if match:
            current = Section(match.group(1), match.group(2))
            sections.append(current)
            bucket = "prose"
            continue
        if GROUP_HEADING_RE.match(line) or HEADING_RE.match(line):
            current = None
            continue
        if current is None:
            continue
        if line.strip() == "Inputs:":
            bucket = "inputs"
            continue
        if line.strip() == "Outputs:":
            bucket = "outputs"
            continue
        bullet = BULLET_RE.match(line)
        if bullet and bucket in ("inputs", "outputs"):
            getattr(current, bucket).append(bullet.group(1))
        elif bucket == "prose":
            current.prose.append(line)
        elif line.strip():
            # Prose after the bullet lists belongs to the description
            current.prose.append(line)
    return sections


def delink(text):
    """Replace Markdown links with plain text, keeping a tool URL as a bare citation."""

    def replace(match):
        label, target = match.group(1), match.group(2)
        label = label.strip("`")
        if target.startswith("http"):
            if "github.com/broadinstitute/gatk-sv" in target:
                return f"`{label}`"
            return f"{label} ({target})"
        if target.endswith("references.md"):
            return label
        if target.startswith("#") or target.endswith(".wdl"):
            return f"`{label}`"
        return f"`{target.lstrip('./')}`"

    return LINK_RE.sub(replace, text)


def strip_emphasis(text):
    return EMPHASIS_RE.sub(lambda m: m.group(1) or m.group(2) or m.group(3), text)


def clean(text):
    """Apply every description rewrite the authoring rules require."""
    text = delink(text)
    text = strip_emphasis(text)
    text = text.replace('"', "'")
    return re.sub(r"\s+", " ", text).strip()


def build_description(section):
    """Fold the section prose into one paragraph per blank-line-separated block."""
    paragraphs = []
    buffer = []
    for raw in section.prose:
        line = raw.strip()
        if not line:
            if buffer:
                paragraphs.append(" ".join(buffer))
                buffer = []
            continue
        note = NOTE_RE.match(line)
        if note:
            if buffer:
                paragraphs.append(" ".join(buffer))
                buffer = []
            buffer.append("Note: " + line[note.end():])
            paragraphs.append(" ".join(buffer))
            buffer = []
            continue
        item = LIST_ITEM_RE.match(line)
        if item:
            buffer.append(f"({item.group(1)})" + " " + line[item.end():])
            continue
        buffer.append(line)
    if buffer:
        paragraphs.append(" ".join(buffer))
    return [clean(p) for p in paragraphs if clean(p)]


def expand_tokens(token, names):
    """Expand a bullet head token into the declaration names it covers.

    Brace and glob forms are resolved case-insensitively, because the document writes
    the mobile-element classes uppercase while the declarations are snake case.
    """
    token = token.split(" ")[-1]
    brace = BRACE_RE.match(token)
    candidates = [token]
    if brace:
        head, body, tail = brace.groups()
        candidates = [head + part + tail for part in body.split(",")]
    resolved = []
    for candidate in candidates:
        if "*" in candidate:
            pattern = re.compile("^" + re.escape(candidate).replace(r"\*", ".*") + "$", re.IGNORECASE)
            resolved.extend(name for name in names if pattern.match(name))
            continue
        lowered = candidate.lower()
        resolved.extend(name for name in names if name.lower() == lowered)
    return resolved


def parse_bullet(bullet, names):
    """Split a bullet into the declaration names it documents and its description text."""
    head, _, tail = bullet.partition(": ")
    if not tail or "`" not in head:
        return [], ""
    tokens = re.findall(r"`([^`]+)`", head)
    if not tokens or re.sub(r"`[^`]+`|,|\band\b|\s", "", head):
        return [], ""
    description = clean(DEFAULT_RE.sub("", tail))
    covered = []
    for token in tokens:
        covered.extend(expand_tokens(token, names))
    return [name for name in covered if name in names], description


def split_pair(covered, description):
    """Give a paired `file, index` bullet a description per name."""
    if len(covered) != 2:
        return {}
    first, second = covered
    if not any(second.endswith(suffix) for suffix in INDEX_SUFFIXES):
        return {}
    text = re.sub(r",? and (?:its )?index(?:es)?\.?$", ".", description, flags=re.IGNORECASE)
    text = re.sub(r",? with (?:its|their) index(?:es)?\.?$", ".", text, flags=re.IGNORECASE)
    return {first: text, second: f"Index for {first}."}


def stub(name, described):
    """Derive a description for a parameter the snapshot never documented."""
    for suffix in INDEX_SUFFIXES:
        if name.endswith(suffix):
            base = name[: -len(suffix)]
            if base in described:
                return f"Index for {base}."
    for prefix, tool in (("gcnv_", "GermlineCNVCaller"), ("ploidy_", "DetermineGermlineContigPloidy")):
        if name.startswith(prefix):
            flag = name[len(prefix):].replace("_", "-")
            return f"{tool} --{flag}."
    return "TODO."


def section_descriptions(section, names):
    """Map the parameters a single section documents to their description text."""
    described = {}
    for bullet in section.inputs + section.outputs:
        covered, description = parse_bullet(bullet, names)
        if not covered or not description:
            continue
        pair = split_pair(covered, description)
        for name in covered:
            described.setdefault(name, pair.get(name, description))
    return described


def shared_descriptions(sections, workflows):
    """Collect descriptions for parameter names that every section agrees on.

    Workflows in this repository reuse parameter names for the same thing, and the
    hand-written document often described a shared input in one section only. A name
    whose sections disagree is left out, so it surfaces as a TODO for an author.
    """
    seen = {}
    for name, section in sections.items():
        entry = workflows.get(name)
        if entry is None:
            continue
        declared = [decl.name for decl in entry[0].declarations()]
        for parameter, description in section_descriptions(section, declared).items():
            seen.setdefault(parameter, set()).add(description)
    return {parameter: texts.pop() for parameter, texts in seen.items() if len(texts) == 1}


def describe(section, workflow, shared):
    """Map every documented parameter of a workflow to its description text."""
    names = [decl.name for decl in workflow.declarations()]
    described = section_descriptions(section, names)

    ordered = {}
    todo = []
    for decl in workflow.inputs:
        if wdl_meta.is_exempt(decl):
            continue
        ordered[decl.name] = described.get(decl.name) or shared.get(decl.name) or stub(decl.name, described)
    for decl in workflow.outputs:
        ordered[decl.name] = described.get(decl.name) or shared.get(decl.name) or stub(decl.name, described)
    for name, text in ordered.items():
        if text == "TODO.":
            todo.append(name)
    return ordered, todo


def render_blocks(description, described):
    """Render the meta and parameter_meta blocks to insert into a WDL file."""
    lines = ["    meta {", "        description: ["]
    for index, paragraph in enumerate(description):
        comma = "," if index < len(description) - 1 else ""
        lines.append(f'            "{paragraph}"{comma}')
    lines.extend(["        ]", "    }", "", "    parameter_meta {"])
    for name, text in described.items():
        lines.append(f'        {name}: "{text}"')
    lines.extend(["    }", ""])
    return lines


def rewrite(path, blocks):
    """Replace any existing meta and parameter_meta blocks with the rendered ones."""
    lines = Path(path).read_text().split("\n")
    output = []
    skipping = False
    inserted = False
    for line in lines:
        if skipping:
            if wdl_meta.BLOCK_CLOSE_RE.match(line):
                skipping = False
            continue
        opened = wdl_meta.BLOCK_OPEN_RE.match(line)
        if opened and opened.group(1) in ("meta", "parameter_meta"):
            skipping = True
            continue
        output.append(line)
        if not inserted and wdl_meta.WORKFLOW_RE.match(line):
            output.extend(blocks)
            inserted = True
    while len(output) > 1 and output[-1] == "" and output[-2] == "":
        output.pop()
    # A replaced block leaves the blank line that followed it, so collapse any run
    collapsed = []
    for line in output:
        if line == "" and collapsed and collapsed[-1] == "":
            continue
        collapsed.append(line)
    return "\n".join(collapsed)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("snapshot", help="hand-written docs/workflows.md to migrate from")
    parser.add_argument("--report", action="store_true", help="report what would change without writing")
    args = parser.parse_args()

    sections = {s.name: s for s in parse_snapshot(Path(args.snapshot).read_text())}
    workflows = {}
    for path in wdl_meta.find_workflow_files():
        workflow = wdl_meta.parse_workflow(path)
        if workflow is not None:
            workflows[workflow.name] = (workflow, path)

    unmatched = sorted(set(sections) - set(workflows))
    missing = sorted(set(workflows) - set(sections))
    shared = shared_descriptions(sections, workflows)
    todos = {}
    written = 0
    for name, (workflow, path) in sorted(workflows.items()):
        section = sections.get(name)
        if section is None:
            continue
        description = build_description(section)
        described, todo = describe(section, workflow, shared)
        if todo:
            todos[name] = todo
        if not args.report:
            Path(path).write_text(rewrite(path, render_blocks(description, described)))
            written += 1

    if unmatched:
        print("Sections with no workflow in this repository (dropped):")
        for name in unmatched:
            print(f"  - {name}")
    if missing:
        print("\nWorkflows with no section in the snapshot (need a description written by hand):")
        for name in missing:
            print(f"  - {name}")
    if todos:
        total = sum(len(names) for names in todos.values())
        print(f"\n{total} parameter(s) need a description written by hand:")
        for name in sorted(todos):
            print(f"  - {name}: {', '.join(sorted(todos[name]))}")
    print(f"\n{'Would write' if args.report else 'Wrote'} {len(workflows) - len(missing) if args.report else written} file(s).")


if __name__ == "__main__":
    main()
