#!/usr/bin/env python3
"""Check that Markdown files follow the style conventions in docs/conventions.md.

Only the mechanically checkable conventions are enforced: whitespace, blank-line
spacing around headings, thematic dividers, and nested list indentation. Rules that
need judgement - wording, section ordering, link targets - are left to review.

Fenced code blocks hold code rather than Markdown, so they are exempt from the
heading and list rules; tabs and trailing whitespace are still checked there.

Usage: check_markdown_style.py [PATH ...]   (defaults to README.md, AGENTS.md and
docs/; directories recurse for *.md). Symlinks are skipped so the AGENTS.md mirrors
at .claude/CLAUDE.md and .github/copilot-instructions.md are reported once.
"""
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PATHS = ["README.md", "AGENTS.md", "docs"]

FENCE_RE = re.compile(r"^\s*(```|~~~)")
HEADING_RE = re.compile(r"^(#{1,6})\s+\S")
LIST_ITEM_RE = re.compile(r"^( *)([-*+]|\d+[.)])(\s+)\S")
DIVIDER_RE = re.compile(r"^\s{0,3}(-{3,}|\*{3,}|_{3,})\s*$")


class Checker:
    """Collect style violations for a single Markdown file."""

    def __init__(self, path, text):
        self.path = path
        self.text = text
        self.lines = text.split("\n")
        self.errors = []
        # Line numbers whose content sits inside a fenced code block
        self.fenced = set()

    def error(self, line_no, code, message):
        self.errors.append((self.path, line_no, code, message))

    def run(self):
        self.mark_fences()
        self.check_whitespace()
        self.check_eof()
        self.check_headings()
        self.check_dividers()
        self.check_list_indentation()
        return self.errors

    def mark_fences(self):
        in_fence = False
        for index, line in enumerate(self.lines):
            if FENCE_RE.match(line):
                self.fenced.add(index)
                in_fence = not in_fence
            elif in_fence:
                self.fenced.add(index)

    def is_blank(self, index):
        return index < 0 or not self.lines[index].strip()

    def next_content_is_heading(self, index):
        for offset in range(index + 1, len(self.lines)):
            if self.lines[offset].strip():
                return bool(HEADING_RE.match(self.lines[offset]))
        return False

    def check_whitespace(self):
        for index, line in enumerate(self.lines):
            if "\t" in line:
                self.error(index + 1, "M001", "tab character; indent with spaces")
            if line != line.rstrip():
                self.error(index + 1, "M002", "trailing whitespace")

    def check_eof(self):
        if not self.text:
            return
        if not self.text.endswith("\n"):
            self.error(len(self.lines), "M003", "file does not end with a newline")
        elif self.text.endswith("\n\n"):
            self.error(len(self.lines), "M003", "file ends with more than one newline")

    def check_headings(self):
        """Enforce the blank-line spacing around headings, tracking `###` position in its `##` section."""
        seen_sub_in_section = False
        for index, line in enumerate(self.lines):
            if index in self.fenced:
                continue
            match = HEADING_RE.match(line)
            if not match:
                continue
            level = len(match.group(1))
            # A heading whose next content is another heading is governed by the blank-line rule below
            if self.is_blank(index + 1) and not self.next_content_is_heading(index):
                self.error(index + 2, "M004", "blank line between a heading and its content")
            if level == 2:
                required = 2
                seen_sub_in_section = False
            elif level == 3:
                required = 2 if not seen_sub_in_section else 1
                seen_sub_in_section = True
            else:
                continue
            if index == 0:
                continue
            blanks = 0
            while self.is_blank(index - blanks - 1) and index - blanks - 1 >= 0:
                blanks += 1
            if blanks != required:
                self.error(
                    index + 1,
                    "M005",
                    "{} blank line(s) before this heading, expected {}".format(blanks, required),
                )

    def check_dividers(self):
        for index, line in enumerate(self.lines):
            if index in self.fenced:
                continue
            if DIVIDER_RE.match(line):
                self.error(index + 1, "M006", "thematic divider; use a heading instead")

    def check_list_indentation(self):
        """Nested list items start at their parent's content column - 2 under `-`, 3 under `1.`."""
        stack = []
        for index, line in enumerate(self.lines):
            if index in self.fenced:
                continue
            match = LIST_ITEM_RE.match(line)
            if not match:
                if not line.strip() or line[:1] == " ":
                    continue
                stack = []
                continue
            indent, marker = match.group(1), match.group(2)
            column = len(indent)
            while stack and stack[-1] > column:
                stack.pop()
            expected = stack[-1] if stack else 0
            if column != expected:
                self.error(
                    index + 1,
                    "M007",
                    "list item indented {}, expected {} to match its parent's content column".format(column, expected),
                )
            stack.append(expected + len(marker) + 1)


def collect_files(paths):
    files = []
    for raw in paths:
        path = Path(raw)
        if not path.is_absolute():
            path = REPO_ROOT / path
        if path.is_dir():
            files.extend(sorted(p for p in path.rglob("*.md") if not p.is_symlink()))
        elif path.is_file():
            files.append(path)
        else:
            sys.exit("check_markdown_style.py: no such file or directory: {}".format(raw))
    return [p for p in files if not p.is_symlink()]


def main():
    paths = sys.argv[1:] or DEFAULT_PATHS
    files = collect_files(paths)
    errors = []
    for path in files:
        rel = str(path.relative_to(REPO_ROOT)) if REPO_ROOT in path.parents else str(path)
        errors.extend(Checker(rel, path.read_text()).run())

    if errors:
        print("Markdown style check failed:\n")
        for rel, line_no, code, message in sorted(errors):
            print("  {}:{}: {} {}".format(rel, line_no, code, message))
        affected = len({e[0] for e in errors})
        print("\n{} error(s) in {} file(s).".format(len(errors), affected))
        sys.exit(1)

    print("Markdown style check passed: {} files.".format(len(files)))


if __name__ == "__main__":
    main()
