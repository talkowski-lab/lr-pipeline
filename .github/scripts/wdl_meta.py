#!/usr/bin/env python3
"""Parse the workflow-level meta, parameter_meta, input and output blocks of a WDL file.

Shared by check_wdl_style.py and generate_workflows_doc.py so the checker and the
generator always read the same structure. Task bodies are skipped: only the blocks
declared directly inside `workflow <Name> {` are parsed.

Declarations are matched strictly. A line inside an input or output block that does
not look like `<Type> <name>[ = <default>]` is reported as an error rather than
silently dropped, so a future syntax change fails loudly instead of quietly losing a
parameter from the generated documentation.
"""
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Directories holding workflow files; wdl/utils additionally holds task libraries
WORKFLOW_DIRS = ["wdl/annotation", "wdl/annotation_utils", "wdl/tools", "wdl/utils"]
META_BLOCKS = ["meta", "parameter_meta", "input", "output"]

WORKFLOW_RE = re.compile(r"^workflow (\w+)\s*\{")
WORKFLOW_SEARCH_RE = re.compile(r"^workflow (\w+)\s*\{", re.MULTILINE)
BLOCK_OPEN_RE = re.compile(r"^ {4}(\w+)\s*\{\s*$")
BLOCK_CLOSE_RE = re.compile(r"^ {4}\}\s*$")
DECL_RE = re.compile(r"^ {8}(\S+) (\w+)(?: = (.*))?$")
PARAM_META_RE = re.compile(r"^ {8}(\w+): (.*)$")
DESCRIPTION_OPEN_RE = re.compile(r"^ {8}description: (.*)$")
STRING_RE = re.compile(r'^"((?:[^"\\]|\\.)*)",?$')

# Inputs the generator documents with a standard bullet, so authors never describe them
DOCKER_INPUT_RE = re.compile(r"^(\w+_)?docker$")
RUNTIME_ATTR_TYPE = "RuntimeAttr?"
PREFIX_INPUT = "prefix"


class Decl:
    """A single workflow-level input or output declaration."""

    def __init__(self, type_name, name, default, line_no):
        self.type = type_name
        self.name = name
        self.default = default
        self.line_no = line_no

    @property
    def optional(self):
        return self.type.endswith("?")

    def __repr__(self):
        return f"Decl({self.type} {self.name})"


class Workflow:
    """The parsed workflow-level structure of a single WDL file."""

    def __init__(self, path, name, line_no):
        self.path = path
        self.name = name
        self.line_no = line_no
        self.description = []
        self.description_is_array = False
        self.param_meta = []
        self.inputs = []
        self.outputs = []
        self.blocks = []
        self.errors = []

    def error(self, line_no, code, message):
        self.errors.append((line_no, code, message))

    @property
    def param_meta_keys(self):
        return [key for key, _, _ in self.param_meta]

    def declarations(self):
        return self.inputs + self.outputs

    def declaration(self, name):
        for decl in self.declarations():
            if decl.name == name:
                return decl
        return None


def unescape(value):
    return value.replace('\\"', '"').replace("\\\\", "\\")


def is_exempt(decl):
    """Report whether the generator documents this input with a standard bullet."""
    if decl.type == RUNTIME_ATTR_TYPE:
        return True
    if decl.type != "String":
        return False
    return decl.name == PREFIX_INPUT or bool(DOCKER_INPUT_RE.match(decl.name))


def parse_workflow(path, text=None):
    """Parse the workflow in a WDL file, or return None when the file declares none."""
    if text is None:
        text = Path(path).read_text()
    lines = text.split("\n")

    workflow = None
    block = None
    in_description = False
    for index, line in enumerate(lines):
        line_no = index + 1

        match = WORKFLOW_RE.match(line)
        if match:
            if workflow is not None:
                workflow.error(line_no, "W026", f"second workflow '{match.group(1)}' in one file")
                return workflow
            workflow = Workflow(str(path), match.group(1), line_no)
            continue
        if workflow is None:
            continue
        if line == "}":
            break

        if in_description:
            in_description = _read_description_line(workflow, line, line_no)
            continue
        if block is None:
            opened = BLOCK_OPEN_RE.match(line)
            if opened and opened.group(1) in META_BLOCKS:
                block = opened.group(1)
                workflow.blocks.append((block, line_no))
            continue
        if BLOCK_CLOSE_RE.match(line):
            block = None
            continue
        if not line.strip() or line.strip().startswith("#"):
            continue

        if block == "meta":
            in_description = _read_meta_line(workflow, line, line_no)
        elif block == "parameter_meta":
            _read_param_meta_line(workflow, line, line_no)
        else:
            _read_decl_line(workflow, block, line, line_no)

    return workflow


def _read_meta_line(workflow, line, line_no):
    """Read one line of a meta block, returning True when a description array is open."""
    match = DESCRIPTION_OPEN_RE.match(line)
    if not match:
        workflow.error(line_no, "W018", "meta may only declare 'description'")
        return False
    value = match.group(1).strip()
    if value == "[":
        workflow.description_is_array = True
        return True
    string = STRING_RE.match(value)
    if string:
        workflow.description.append(unescape(string.group(1)))
    else:
        workflow.error(line_no, "W018", "meta description must be an array of strings")
    return False


def _read_description_line(workflow, line, line_no):
    """Read one line inside an open description array, returning False at its close."""
    value = line.strip()
    if value == "]":
        return False
    string = STRING_RE.match(value)
    if string:
        workflow.description.append(unescape(string.group(1)))
    else:
        workflow.error(line_no, "W018", "meta description entries must be double-quoted strings")
    return True


def _read_param_meta_line(workflow, line, line_no):
    match = PARAM_META_RE.match(line)
    if not match:
        workflow.error(line_no, "W025", "unparseable parameter_meta line")
        return
    string = STRING_RE.match(match.group(2).strip())
    if not string:
        workflow.error(line_no, "W025", f"parameter_meta '{match.group(1)}' must be a double-quoted string")
        return
    workflow.param_meta.append((match.group(1), unescape(string.group(1)), line_no))


def _read_decl_line(workflow, block, line, line_no):
    match = DECL_RE.match(line)
    if not match:
        workflow.error(line_no, "W025", f"unparseable {block} declaration")
        return
    decl = Decl(match.group(1), match.group(2), match.group(3), line_no)
    if block == "input":
        workflow.inputs.append(decl)
    else:
        workflow.outputs.append(decl)


def find_workflow_files():
    """Return every WDL file under WORKFLOW_DIRS that declares a workflow, sorted by name."""
    files = []
    for directory in WORKFLOW_DIRS:
        for path in sorted((REPO_ROOT / directory).glob("*.wdl")):
            if WORKFLOW_SEARCH_RE.search(path.read_text()):
                files.append(path)
    return files
