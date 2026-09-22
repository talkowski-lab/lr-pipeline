#!/usr/bin/env python3
"""Check that WDL files follow the style conventions in docs/conventions.md.

Only the mechanically checkable conventions are enforced: whitespace, 4-space
indentation, section ordering within a task, and workflow/task naming. Rules that
need judgement - input grouping, comment wording, runtime sizing - are left to review.

The body of a `command <<< >>>` block holds Bash and Python rather than WDL, so it is
exempt from the indentation rules; a nested heredoc payload such as `python3 <<'CODE'`
deliberately starts at column 0. Tabs and trailing whitespace are still checked there.

Usage: check_wdl_style.py [PATH ...]   (defaults to wdl/; directories recurse for *.wdl)
"""
import re
import sys
from pathlib import Path

import wdl_meta

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PATHS = ["wdl"]
INDENT_UNIT = 4
TASK_SECTIONS = ["input", "command", "output", "runtime"]
REQUIRED_TASK_SECTIONS = ["command", "runtime"]
DOCUMENTED_BLOCKS = ["meta", "parameter_meta", "input"]

BLOCK_RE = re.compile(r"^(task|workflow|struct)\s+(\w+)")
COMMAND_OPEN_RE = re.compile(r"^(\s*)command\s*<<<\s*$")
COMMAND_BRACE_RE = re.compile(r"^\s*command\s*\{")
SECTION_RE = re.compile(r"^ {4}(input|output|runtime)\s*\{")
IMPORT_AS_RE = re.compile(r"^import\s+.* as ")
BLANK_COMMENT_RE = re.compile(r"^\s*#{2,}\s*$")
PASCAL_CASE_RE = re.compile(r"^[A-Z][A-Za-z0-9]*$")
CALL_ALIAS_RE = re.compile(r"^\s*call\s+[\w.]+\s+as\s+(\w+)")
RUNTIME_ATTR_INPUT_RE = re.compile(r"^ {8}RuntimeAttr\?\s+(\w+)\s*$")

# Markdown that docs/workflows.md is generated from, so descriptions must not carry it
MARKDOWN_LINK_RE = re.compile(r"\]\(")
MARKDOWN_EMPHASIS_RE = re.compile(r"\*")
MARKDOWN_PREFIX_RE = re.compile(r"^(>|[-*+] |#{1,6} |\d+[.)] )")
# A backticked span holds code, where `*` is multiplication rather than emphasis
CODE_SPAN_RE = re.compile(r"`[^`]*`")
URL_RE = re.compile(r"https?://")
# A tool URL is written as `ToolName (https://...)`, the one link form descriptions allow
TOOL_URL_RE = re.compile(r"\(https?://[^\s()]+\)")

# Tasks whose command block deliberately opens with something other than `set -euo pipefail`.
# Add an entry only when the standard flags would change runtime behavior, and say why.
SET_LINE_ALLOWLIST = set()


class Checker:
    """Collect style violations for a single WDL file."""

    def __init__(self, path, text):
        self.path = path
        self.text = text
        self.lines = text.split("\n")
        if self.lines and self.lines[-1] == "":
            self.lines.pop()
        self.errors = []
        self.block_kind = None
        self.block_name = None
        self.sections = []
        self.in_command = False
        self.command_indent = 0
        self.command_body_started = False

    def error(self, line_no, code, message):
        self.errors.append((self.path, line_no, code, message))

    def run(self):
        self.check_file_level()
        for index, line in enumerate(self.lines):
            self.check_line(index, line)
        self.flush_task(len(self.lines))
        return self.errors

    def check_file_level(self):
        if "\t" in self.text:
            for index, line in enumerate(self.lines):
                if "\t" in line:
                    self.error(index + 1, "W001", "tab character; indent with 4 spaces")
        if not self.text.endswith("\n") or self.text.endswith("\n\n"):
            self.error(len(self.lines), "W004", "file must end with exactly one newline")

        # Locate `version 1.0`, skipping an exempt provenance or license header
        for index, line in enumerate(self.lines):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if line.rstrip() != "version 1.0":
                self.error(index + 1, "W011", "first non-comment line must be 'version 1.0'")
            elif index + 1 < len(self.lines) and self.lines[index + 1].strip():
                self.error(index + 2, "W011", "'version 1.0' must be followed by a blank line")
            break

        workflows = [m.group(2) for m in (BLOCK_RE.match(ln) for ln in self.lines) if m and m.group(1) == "workflow"]
        stem = Path(self.path).stem
        for name in workflows:
            if name != stem:
                self.error(1, "W013", "workflow '{}' must match the file name '{}.wdl'".format(name, stem))

        self.check_unused_runtime_attrs()
        if workflows:
            self.check_workflow_meta()
        else:
            self.check_alphabetical_tasks()

    def check_alphabetical_tasks(self):
        """In a task library - a file with no workflow - tasks must be declared alphabetically."""
        declared = []
        in_command = False
        for index, line in enumerate(self.lines):
            if in_command:
                if line.strip() == ">>>":
                    in_command = False
                continue
            if COMMAND_OPEN_RE.match(line):
                in_command = True
                continue
            block_match = BLOCK_RE.match(line)
            if block_match and block_match.group(1) == "task":
                declared.append((block_match.group(2), index + 1))
        for (name, line_no), (previous, _) in zip(declared[1:], declared):
            if name < previous:
                self.error(line_no, "W017", "task '{}' must be declared before '{}'".format(name, previous))

    def check_unused_runtime_attrs(self):
        """Flag workflow-level RuntimeAttr? inputs that no call passes as runtime_attr_override."""
        in_workflow = False
        for index, line in enumerate(self.lines):
            block_match = BLOCK_RE.match(line)
            if block_match:
                in_workflow = block_match.group(1) == "workflow"
                continue
            if not in_workflow:
                continue
            attr_match = RUNTIME_ATTR_INPUT_RE.match(line)
            if not attr_match:
                continue
            name = attr_match.group(1)
            uses = len(re.findall(r"\b" + re.escape(name) + r"\b", self.text))
            if uses < 2:
                self.error(index + 1, "W016", "nothing passes RuntimeAttr? input '{}' to a call".format(name))

    def check_workflow_meta(self):
        """Enforce the meta and parameter_meta rules that docs/workflows.md is generated from."""
        if not any(self.path.startswith(directory + "/") for directory in wdl_meta.WORKFLOW_DIRS):
            return
        workflow = wdl_meta.parse_workflow(self.path, self.text)
        if workflow is None:
            return
        for line_no, code, message in workflow.errors:
            self.error(line_no, code, message)

        present = [name for name, _ in workflow.blocks if name in DOCUMENTED_BLOCKS]
        if present != DOCUMENTED_BLOCKS:
            self.error(
                workflow.line_no,
                "W019",
                "workflow '{}' must open with {}, found {}".format(
                    workflow.name, " -> ".join(DOCUMENTED_BLOCKS), " -> ".join(present) or "none"
                ),
            )
        self.check_description(workflow)
        self.check_parameter_meta(workflow)

    def block_line(self, workflow, name):
        return next((line_no for block, line_no in workflow.blocks if block == name), workflow.line_no)

    def check_description(self, workflow):
        line_no = self.block_line(workflow, "meta")
        if not workflow.description:
            self.error(line_no, "W018", "workflow '{}' has no meta description".format(workflow.name))
            return
        if not workflow.description_is_array:
            self.error(line_no, "W018", "meta description must be an array of strings, one per paragraph")
        for index, paragraph in enumerate(workflow.description):
            self.check_prose(line_no, "description paragraph {}".format(index + 1), paragraph)

    def check_prose(self, line_no, label, text):
        """Reject the Markdown and quoting that the generator cannot carry into the document."""
        if not text.strip():
            self.error(line_no, "W024", "{} is empty".format(label))
            return
        prose = CODE_SPAN_RE.sub("", text)
        if MARKDOWN_LINK_RE.search(prose):
            self.error(line_no, "W024", "{} must not contain a Markdown link".format(label))
        if MARKDOWN_EMPHASIS_RE.search(prose):
            self.error(line_no, "W024", "{} must not contain Markdown emphasis".format(label))
        if MARKDOWN_PREFIX_RE.match(text):
            self.error(line_no, "W024", "{} must be a plain paragraph".format(label))
        if '"' in text:
            self.error(line_no, "W024", "{} must use single quotes, not escaped double quotes".format(label))
        if len(URL_RE.findall(text)) != len(TOOL_URL_RE.findall(text)):
            self.error(line_no, "W024", "{} may only cite a URL as 'ToolName (https://...)'".format(label))

    def check_parameter_meta(self, workflow):
        """Every non-exempt input and output needs one entry, in declaration order."""
        exempt = {decl.name for decl in workflow.inputs if wdl_meta.is_exempt(decl)}
        expected = [decl.name for decl in workflow.inputs if not wdl_meta.is_exempt(decl)]
        expected += [decl.name for decl in workflow.outputs]

        documented = []
        for key, value, line_no in workflow.param_meta:
            self.check_prose(line_no, "parameter_meta '{}'".format(key), value)
            if key in exempt:
                self.error(line_no, "W023", "'{}' is documented automatically; remove this entry".format(key))
            elif key in documented:
                self.error(line_no, "W020", "duplicate parameter_meta entry '{}'".format(key))
            elif workflow.declaration(key) is None:
                self.error(line_no, "W020", "parameter_meta '{}' is not a declared input or output".format(key))
            else:
                documented.append(key)

        for name in expected:
            if name not in documented:
                self.error(workflow.declaration(name).line_no, "W022", "'{}' has no parameter_meta entry".format(name))

        wanted = [name for name in expected if name in documented]
        for found, want in zip(documented, wanted):
            if found != want:
                line_no = next(entry[2] for entry in workflow.param_meta if entry[0] == found)
                self.error(
                    line_no,
                    "W021",
                    "parameter_meta '{}' is out of order; expected '{}' here".format(found, want),
                )
                break

    def check_line(self, index, line):
        line_no = index + 1
        if line != line.rstrip():
            self.error(line_no, "W002", "trailing whitespace")

        if self.in_command:
            self.check_command_body(index, line)
            return

        self.check_blank_run(index, line)

        open_match = COMMAND_OPEN_RE.match(line)
        if open_match:
            self.enter_command(index, line, open_match)
            return
        if COMMAND_BRACE_RE.match(line):
            self.error(line_no, "W006", "command block must use the '<<< ... >>>' heredoc form")

        if line.strip():
            indent = len(line) - len(line.lstrip(" "))
            if indent % INDENT_UNIT:
                self.error(line_no, "W005", "indent of {} spaces is not a multiple of 4".format(indent))

        if BLANK_COMMENT_RE.match(line):
            self.error(line_no, "W014", "blank comment line")
        if IMPORT_AS_RE.match(line):
            self.error(line_no, "W010", "imports must not be renamed with 'as'")

        alias_match = CALL_ALIAS_RE.match(line)
        if alias_match and not PASCAL_CASE_RE.match(alias_match.group(1)):
            self.error(line_no, "W015", "call alias '{}' must be in Pascal case".format(alias_match.group(1)))

        block_match = BLOCK_RE.match(line)
        if block_match:
            self.flush_task(line_no)
            self.block_kind, self.block_name = block_match.group(1), block_match.group(2)
            if not PASCAL_CASE_RE.match(self.block_name):
                self.error(line_no, "W012", "{} '{}' must be in Pascal case".format(self.block_kind, self.block_name))
            return

        section_match = SECTION_RE.match(line)
        if section_match and self.block_kind == "task":
            self.sections.append((section_match.group(1), line_no))

    def enter_command(self, index, line, open_match):
        line_no = index + 1
        self.in_command = True
        self.command_indent = len(open_match.group(1))
        self.command_body_started = False
        if self.command_indent % INDENT_UNIT:
            self.error(line_no, "W005", "indent of {} spaces is not a multiple of 4".format(self.command_indent))
        if self.block_kind == "task":
            self.sections.append(("command", line_no))
        previous = self.lines[index - 1].strip() if index else ""
        if previous and not previous.startswith("#"):
            self.error(line_no, "W009", "'command <<<' must be preceded by a blank line or a comment")

    def check_command_body(self, index, line):
        line_no = index + 1
        if line.strip() == ">>>":
            self.in_command = False
            indent = len(line) - len(line.lstrip(" "))
            if indent != self.command_indent:
                self.error(line_no, "W005", "'>>>' must align with 'command <<<'")
            if index + 1 < len(self.lines) and self.lines[index + 1].strip():
                self.error(line_no + 1, "W009", "'>>>' must be followed by a blank line")
            return
        if self.command_body_started or not line.strip():
            return
        self.command_body_started = True
        indent = len(line) - len(line.lstrip(" "))
        if indent != self.command_indent + INDENT_UNIT:
            self.error(line_no, "W006", "command body must be indented 4 spaces past 'command <<<'")
        if (self.path, self.block_name) in SET_LINE_ALLOWLIST:
            return
        if line.strip() != "set -euo pipefail":
            self.error(line_no, "W007", "command block must begin with 'set -euo pipefail'")
        elif index + 1 < len(self.lines) and self.lines[index + 1].strip():
            self.error(line_no + 1, "W007", "'set -euo pipefail' must be followed by a blank line")

    def check_blank_run(self, index, line):
        if line.strip() or index == 0:
            return
        if not self.lines[index - 1].strip():
            self.error(index + 1, "W003", "consecutive blank lines")

    def flush_task(self, line_no):
        if self.block_kind != "task":
            self.reset_task()
            return
        seen = {}
        for name, section_line in self.sections:
            if name in seen:
                self.error(section_line, "W008", "task '{}' has more than one '{}' block".format(self.block_name, name))
            else:
                seen[name] = section_line
        for name in REQUIRED_TASK_SECTIONS:
            if name not in seen:
                self.error(line_no, "W008", "task '{}' is missing a '{}' block".format(self.block_name, name))
        ordered = [(seen[n], n) for n in TASK_SECTIONS if n in seen]
        if ordered != sorted(ordered):
            expected = " -> ".join(n for n in TASK_SECTIONS if n in seen)
            found = " -> ".join(n for _, n in sorted(ordered))
            self.error(
                ordered[0][0],
                "W008",
                "task '{}' orders its blocks {} but must use {}".format(self.block_name, found, expected),
            )
        self.reset_task()

    def reset_task(self):
        self.block_kind = None
        self.block_name = None
        self.sections = []


def collect_files(paths):
    files = []
    for raw in paths:
        path = Path(raw)
        if not path.is_absolute():
            path = REPO_ROOT / path
        if path.is_dir():
            files.extend(sorted(path.rglob("*.wdl")))
        elif path.is_file():
            files.append(path)
        else:
            sys.exit("check_wdl_style.py: no such file or directory: {}".format(raw))
    return files


def main():
    paths = sys.argv[1:] or DEFAULT_PATHS
    files = collect_files(paths)
    errors = []
    for path in files:
        rel = str(path.relative_to(REPO_ROOT)) if REPO_ROOT in path.parents else str(path)
        errors.extend(Checker(rel, path.read_text()).run())

    if errors:
        print("WDL style check failed:\n")
        for rel, line_no, code, message in sorted(errors):
            print("  {}:{}: {} {}".format(rel, line_no, code, message))
        affected = len({e[0] for e in errors})
        print("\n{} error(s) in {} file(s).".format(len(errors), affected))
        sys.exit(1)

    print("WDL style check passed: {} files.".format(len(files)))


if __name__ == "__main__":
    main()
