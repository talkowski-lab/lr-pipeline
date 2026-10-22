#!/usr/bin/env python3
"""Check that .claude/CLAUDE.md and .github/copilot-instructions.md stay in sync with AGENTS.md.

Both files are expected to be symlinks to ../AGENTS.md, but this checks byte-for-byte
content equality rather than symlink-ness, so it also catches the case where a tool
(e.g. GitHub's web editor, a symlink-unaware checkout) replaces the symlink with a
materialized file whose content has drifted from AGENTS.md.
"""

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CANONICAL = REPO_ROOT / "AGENTS.md"
MIRRORS = [REPO_ROOT / ".claude" / "CLAUDE.md", REPO_ROOT / ".github" / "copilot-instructions.md"]


def main():
    errors = []
    canonical_content = CANONICAL.read_text()

    for mirror in MIRRORS:
        rel = mirror.relative_to(REPO_ROOT)
        if not mirror.is_file():
            errors.append(f"{rel} is missing")
            continue
        if mirror.read_text() != canonical_content:
            errors.append(f"{rel} content has drifted from AGENTS.md")

    if errors:
        print("Agents sync check failed:\n")
        for e in errors:
            print(f"  - {e}")
        print(f"\n{len(errors)} error(s) found. Each mirror must be a symlink to ../AGENTS.md (or an exact copy of its content).")
        sys.exit(1)

    print("Agents sync check passed: both mirrors match AGENTS.md.")


if __name__ == "__main__":
    main()
