#!/usr/bin/env bash
set -euo pipefail

# Merges a validated feature branch into main, removes its worktree and deletes it, following the sequence in docs/ci-cd.md.
# Stops without merging if the rebase conflicts, leaving the conflict in place in the worktree for a human to resolve.

usage() {
    echo "Usage: $(basename "$0") [<branch>]" >&2
    echo "Merges <branch> into main. Run from the main checkout with the branch name, or from inside the branch's worktree." >&2
    exit 2
}

if [[ $# -gt 1 || ${1:-} == -* ]]; then
    usage
fi

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

MAIN_WORKTREE="$(git worktree list --porcelain | sed -n '1s/^worktree //p')"
BRANCH_WORKTREE="$(git worktree list --porcelain |
    awk -v ref="refs/heads/$BRANCH" '$1 == "worktree" {path = $2} $1 == "branch" && $2 == ref {print path}')"

if [[ -z "$BRANCH_WORKTREE" ]]; then
    echo "'$BRANCH' is not checked out in any worktree; run .github/scripts/new_worktree.sh $BRANCH first." >&2
    exit 1
fi

if [[ "$BRANCH_WORKTREE" == "$MAIN_WORKTREE" ]]; then
    echo "'$BRANCH' is checked out in the main checkout, which must stay on main." >&2
    exit 1
fi

if [[ "$(git -C "$MAIN_WORKTREE" symbolic-ref --short HEAD)" != "main" ]]; then
    echo "Main checkout $MAIN_WORKTREE must stay on main; run 'git switch main' there first." >&2
    exit 1
fi

cd "$BRANCH_WORKTREE"

if [[ -n "$(git status --porcelain)" ]]; then
    echo "Worktree $BRANCH_WORKTREE is not clean; commit or stash before merging." >&2
    exit 1
fi

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

git -C "$MAIN_WORKTREE" pull --ff-only origin main
git -C "$MAIN_WORKTREE" merge --ff-only "$BRANCH"
git -C "$MAIN_WORKTREE" push origin main

# Leave the worktree before removing it, then delete the branch on both sides,
# since deleting the remote branch is what removes its Dockstore version
cd "$MAIN_WORKTREE"
if [[ -n "$(git ls-remote --heads origin "$BRANCH")" ]]; then
    git push origin --delete "$BRANCH"
fi
git worktree remove "$BRANCH_WORKTREE"
git branch -d "$BRANCH"

# Drop the branch-tagged images build_docker.sh pushed, last so that a gcloud failure leaves git state final
REGISTRY="us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers"
for image in $(gcloud artifacts docker images list "$REGISTRY" --include-tags --format='value(package)' \
    --filter="tags:$BRANCH" 2>/dev/null | sort -u); do
    echo "Deleting $image:$BRANCH"
    gcloud artifacts docker tags delete "$image:$BRANCH" --quiet
done

echo "Merged $BRANCH into main, removed its worktree and deleted it."
