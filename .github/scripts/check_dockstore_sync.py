#!/usr/bin/env python3
"""Check that .dockstore.yml and wdl/ stay in sync.

Every top-level workflow file under wdl/annotation, wdl/annotation_utils, and
wdl/tools must have exactly one corresponding entry in .dockstore.yml, and every
entry must point to a file that exists. Every entry must also list main under
filters.branches; --main-only additionally fails on any feature branch left
listed there, and is run before merging a feature branch into main.
"""
import argparse
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCKSTORE_YML = REPO_ROOT / ".dockstore.yml"
WORKFLOW_DIRS = ["wdl/annotation", "wdl/annotation_utils", "wdl/tools"]


def load_dockstore_entries():
    with open(DOCKSTORE_YML) as f:
        data = yaml.safe_load(f)
    return [w for w in data["workflows"] if w.get("subclass") == "WDL"]


def find_active_workflow_files():
    files = set()
    for d in WORKFLOW_DIRS:
        files.update(p for p in (REPO_ROOT / d).glob("*.wdl"))
    return files


def check_branch_filters(name, entry, main_only):
    branches = entry.get("filters", {}).get("branches")
    if not isinstance(branches, list):
        return [f".dockstore.yml entry '{name}' has no filters.branches list"]
    errors = []
    if "main" not in branches:
        errors.append(f".dockstore.yml entry '{name}' does not list 'main' under filters.branches")
    if main_only:
        extra = [b for b in branches if b != "main"]
        if extra:
            errors.append(
                f".dockstore.yml entry '{name}' still lists feature branch(es) under filters.branches: "
                f"{', '.join(extra)}"
            )
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--main-only",
        action="store_true",
        help="also fail if any entry lists a branch other than main",
    )
    args = parser.parse_args()

    errors = []
    entries = load_dockstore_entries()
    entries_by_name = {}
    for entry in entries:
        name = entry["name"]
        if name in entries_by_name:
            errors.append(f"duplicate .dockstore.yml entry for '{name}'")
        entries_by_name[name] = entry

    for name, entry in entries_by_name.items():
        path = REPO_ROOT / entry["primaryDescriptorPath"].lstrip("/")
        if not path.is_file():
            errors.append(
                f".dockstore.yml entry '{name}' points to missing file: {entry['primaryDescriptorPath']}"
            )
        errors.extend(check_branch_filters(name, entry, args.main_only))

    active_files = find_active_workflow_files()
    for path in sorted(active_files):
        name = path.stem
        entry = entries_by_name.get(name)
        if entry is None:
            errors.append(f"no .dockstore.yml entry found for active workflow file: {path.relative_to(REPO_ROOT)}")
            continue
        expected_path = "/" + str(path.relative_to(REPO_ROOT))
        if entry["primaryDescriptorPath"] != expected_path:
            errors.append(
                f".dockstore.yml entry '{name}' has primaryDescriptorPath "
                f"'{entry['primaryDescriptorPath']}', expected '{expected_path}'"
            )

    active_names = {p.stem for p in active_files}
    for name in entries_by_name:
        if name not in active_names:
            errors.append(f".dockstore.yml entry '{name}' has no matching active workflow file under wdl/")

    if errors:
        print("Dockstore sync check failed:\n")
        for e in errors:
            print(f"  - {e}")
        print(f"\n{len(errors)} error(s) found.")
        sys.exit(1)

    print(
        f"Dockstore sync check passed: {len(active_files)} active workflow files, "
        f"{len(entries_by_name)} .dockstore.yml entries."
    )


if __name__ == "__main__":
    main()
