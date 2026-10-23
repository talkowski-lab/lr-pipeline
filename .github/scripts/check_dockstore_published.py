#!/usr/bin/env python3
"""Check that every workflow in .dockstore.yml is published on Dockstore.

Dockstore only serves workflows to Terra once they are published, and registering an entry in
.dockstore.yml does not publish it. Each name is looked up through the public Dockstore API: 200 means
published and 404 means not. Server errors and timeouts are retried. --added-since <git-ref> downgrades
entries added to .dockstore.yml since that ref to warnings, since Dockstore needs time to sync new entries.

Usage: check_dockstore_published.py [--added-since GIT_REF]
"""

import argparse
import re
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCKSTORE_YML = REPO_ROOT / ".dockstore.yml"
API_URL = "https://dockstore.org/api/workflows/path/workflow/{path}/published"
REPOSITORY = "github.com/talkowski-lab/lr-pipeline"
ATTEMPTS = 3
RETRY_DELAY_SECONDS = 5
TIMEOUT_SECONDS = 30
ADDED_NAME_RE = re.compile(r"^\+\s*-?\s*name:\s*(\w+)\s*$")


def load_workflow_names():
    with open(DOCKSTORE_YML) as f:
        data = yaml.safe_load(f)
    return [w["name"] for w in data["workflows"] if w.get("subclass") == "WDL"]


def find_added_names(ref):
    """Names of entries that .dockstore.yml gained since the given git ref."""
    diff = subprocess.run(
        ["git", "diff", "--unified=0", ref, "--", ".dockstore.yml"],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return {m.group(1) for m in map(ADDED_NAME_RE.match, diff.splitlines()) if m}


def is_published(name):
    """True if Dockstore serves the workflow as published, False on 404; raises after repeated failures."""
    url = API_URL.format(path=urllib.parse.quote(f"{REPOSITORY}/{name}", safe=""))
    for attempt in range(1, ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT_SECONDS):
                return True
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return False
            if exc.code < 500:
                raise
            failure = exc
        except (urllib.error.URLError, TimeoutError) as exc:
            failure = exc
        if attempt == ATTEMPTS:
            raise RuntimeError(f"could not query Dockstore for '{name}': {failure}")
        time.sleep(RETRY_DELAY_SECONDS)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--added-since",
        metavar="GIT_REF",
        help="report entries added to .dockstore.yml since this ref as warnings instead of errors",
    )
    args = parser.parse_args()

    names = load_workflow_names()
    added = find_added_names(args.added_since) if args.added_since else set()
    unpublished = [n for n in names if not is_published(n)]
    errors = [n for n in unpublished if n not in added]
    warnings = [n for n in unpublished if n in added]

    for name in warnings:
        print(f"warning: '{name}' was added recently and is not published yet")
    if errors:
        print("Dockstore published check failed:\n")
        for name in errors:
            print(f"  - '{name}' is not published on Dockstore")
        print(f"\n{len(errors)} unpublished workflow(s) of {len(names)}.")
        sys.exit(1)

    print(f"Dockstore published check passed: {len(names) - len(warnings)} of {len(names)} workflows published.")


if __name__ == "__main__":
    main()
