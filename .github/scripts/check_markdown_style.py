#!/usr/bin/env python3
"""Check that Markdown files follow the style conventions in docs/conventions.md.

Only the mechanically checkable conventions are enforced: whitespace, blank-line
spacing around headings, thematic dividers, nested list indentation, and that every
Markdown file is referenced from somewhere. Rules that need judgement - wording,
section ordering, link targets - are left to review.

Fenced code blocks hold code rather than Markdown, so they are exempt from the
heading and list rules; tabs and trailing whitespace are still checked there.

Usage: check_markdown_style.py [PATH ...]   (defaults to README.md, AGENTS.md,
CONTRIBUTING.md, CODE_OF_CONDUCT.md, SECURITY.md, .github/PULL_REQUEST_TEMPLATE.md,
docs/ and tests/; directories recurse for *.md). Symlinks are skipped so the AGENTS.md
mirrors at .claude/CLAUDE.md and .github/copilot-instructions.md are reported once.
With no PATH, the repository-wide reference rule also runs over every tracked *.md.
"""

import posixpath
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PATHS = [
    "README.md",
    "AGENTS.md",
    "CONTRIBUTING.md",
    "CODE_OF_CONDUCT.md",
    "SECURITY.md",
    ".github/PULL_REQUEST_TEMPLATE.md",
    "docs",
    "tests",
]

FENCE_RE = re.compile(r"^\s*(```|~~~)")
HEADING_RE = re.compile(r"^(#{1,6})\s+\S")
LIST_ITEM_RE = re.compile(r"^( *)([-*+]|\d+[.)])(\s+)\S")
DIVIDER_RE = re.compile(r"^\s{0,3}(-{3,}|\*{3,}|_{3,})\s*$")
MD_MENTION_RE = re.compile(r"[\w./-]*\.md(?!\w)")


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


def tracked_texts(root):
    """Read every tracked text file under root that is not a symlink, keyed by repository-relative path."""
    listing = subprocess.run(["git", "ls-files", "-z"], cwd=str(root), check=True, capture_output=True, text=True).stdout
    texts = {}
    for rel in filter(None, listing.split("\0")):
        path = root / rel
        if path.is_symlink() or not path.is_file():
            continue
        try:
            texts[rel] = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
    return texts


def find_unreferenced(texts):
    """Return the Markdown files in texts ({repository-relative path: text}) that no other file mentions.

    A mention is any path ending in `.md` that resolves to the file from the repository root or, as a
    relative link does, from the mentioning file's directory. A bare basename counts too when no other
    Markdown file shares it, so `README.md` alone never vouches for a nested README. The root README.md
    is the repository's landing page and needs no referrer.
    """
    by_name = {}
    for path in texts:
        if path.endswith(".md"):
            by_name.setdefault(posixpath.basename(path), set()).add(path)
    referenced = set()
    for referrer, text in texts.items():
        base = posixpath.dirname(referrer)
        for mention in set(MD_MENTION_RE.findall(text)):
            candidates = {posixpath.normpath(mention.lstrip("/")), posixpath.normpath(posixpath.join(base, mention))}
            same_name = by_name.get(posixpath.basename(mention), set())
            if len(same_name) == 1:
                candidates |= same_name
            referenced |= candidates - {referrer}
    targets = {path for path in texts if path.endswith(".md") and path != "README.md"}
    return sorted(targets - referenced)


def check_references(root):
    return [
        (rel, 1, "M008", "not referenced by any other tracked file; link or mention it from a relevant document")
        for rel in find_unreferenced(tracked_texts(root))
    ]


def main():
    paths = sys.argv[1:] or DEFAULT_PATHS
    files = collect_files(paths)
    errors = []
    for path in files:
        rel = str(path.relative_to(REPO_ROOT)) if REPO_ROOT in path.parents else str(path)
        errors.extend(Checker(rel, path.read_text()).run())
    if not sys.argv[1:]:
        errors.extend(check_references(REPO_ROOT))

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
