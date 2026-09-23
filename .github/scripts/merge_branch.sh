#!/usr/bin/env bash
set -euo pipefail

# Merges a validated feature branch into main and deletes it, following the sequence in docs/ci-cd.md.
# Stops without merging if the rebase conflicts, leaving the conflict in place for a human to resolve.

usage() {
    echo "Usage: $(basename "$0") [<branch>]" >&2
    echo "Merges <branch>, defaulting to the current branch, into main. Run from anywhere in the repository." >&2
    exit 2
}

if [[ $# -gt 1 || ${1:-} == -* ]]; then
    usage
fi

cd "$(git rev-parse --show-toplevel)"
BRANCH="${1:-$(git rev-parse --abbrev-ref HEAD)}"

if [[ "$BRANCH" == "main" ]]; then
    echo "Refusing to merge main into itself; check out the feature branch first." >&2
    exit 1
fi

# Reject characters that would change the meaning of the sed address used to edit .dockstore.yml
if [[ ! "$BRANCH" =~ ^[A-Za-z0-9._-]+$ ]]; then
    echo "Branch name '$BRANCH' is not of the form kj-<kebab-topic>." >&2
    exit 1
fi

if ! git show-ref --verify --quiet "refs/heads/$BRANCH"; then
    echo "No local branch named '$BRANCH'." >&2
    exit 1
fi

if [[ -n "$(git status --porcelain)" ]]; then
    echo "Working tree is not clean; commit or stash before merging." >&2
    exit 1
fi

git switch "$BRANCH"

# Drop the branch from every .dockstore.yml filter so main never carries a feature-branch version
if grep -qE "^[[:space:]]*-[[:space:]]*${BRANCH}[[:space:]]*$" .dockstore.yml; then
    TMP_DOCKSTORE="$(mktemp)"
    sed -E "/^[[:space:]]*-[[:space:]]*${BRANCH}[[:space:]]*$/d" .dockstore.yml > "$TMP_DOCKSTORE"
    mv "$TMP_DOCKSTORE" .dockstore.yml
    python3 .github/scripts/check_dockstore_sync.py --main-only
    git add .dockstore.yml
    git commit -m "Remove $BRANCH Dockstore filter before merge"
fi

# Rebase onto current main so that the merge can fast-forward
git fetch origin main
if ! git rebase origin/main; then
    echo >&2
    echo "Rebasing '$BRANCH' onto main hit a conflict and nothing has been merged." >&2
    echo "Resolve the conflicts, finish with 'git rebase --continue', then run this script again." >&2
    echo "Run 'git rebase --abort' to back out instead." >&2
    exit 1
fi

git switch main
git pull --ff-only origin main
git merge --ff-only "$BRANCH"
git push origin main

# Delete the branch on both sides, since deleting the remote branch is what removes its Dockstore version
if [[ -n "$(git ls-remote --heads origin "$BRANCH")" ]]; then
    git push origin --delete "$BRANCH"
fi
git branch -d "$BRANCH"

echo "Merged $BRANCH into main and deleted it."
