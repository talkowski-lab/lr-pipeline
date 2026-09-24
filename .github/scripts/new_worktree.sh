#!/usr/bin/env bash
set -euo pipefail

# Creates the linked worktree a branch is developed in, beside the main checkout, following docs/conventions.md.
# Prints the worktree path on stdout; reuses the existing worktree when the branch already has one.

usage() {
    echo "Usage: $(basename "$0") <branch>" >&2
    echo "Creates or reuses the worktree for <branch>, named kj-<kebab-topic>. Run from anywhere in the repository." >&2
    exit 2
}

if [[ $# -ne 1 || ${1:-} == -* ]]; then
    usage
fi

BRANCH="$1"

if [[ ! "$BRANCH" =~ ^kj-[a-z0-9]+(-[a-z0-9]+)*$ ]]; then
    echo "Branch name '$BRANCH' is not of the form kj-<kebab-topic>." >&2
    exit 1
fi

MAIN_WORKTREE="$(git worktree list --porcelain | sed -n '1s/^worktree //p')"

# The main checkout is the shared hub every worktree fetches and fast-forwards through, so it stays on main
if [[ "$(git -C "$MAIN_WORKTREE" symbolic-ref --short HEAD)" != "main" ]]; then
    echo "Main checkout $MAIN_WORKTREE must stay on main; run 'git switch main' there first." >&2
    exit 1
fi

EXISTING="$(git worktree list --porcelain |
    awk -v ref="refs/heads/$BRANCH" '$1 == "worktree" {path = $2} $1 == "branch" && $2 == ref {print path}')"

if [[ -n "$EXISTING" ]]; then
    WORKTREE="$EXISTING"
    echo "Reusing existing worktree for $BRANCH." >&2
else
    WORKTREE="$(dirname "$MAIN_WORKTREE")/$(basename "$MAIN_WORKTREE")-worktrees/$BRANCH"
    git fetch origin >&2

    if git show-ref --verify --quiet "refs/heads/$BRANCH"; then
        git worktree add "$WORKTREE" "$BRANCH" >&2
    elif [[ -n "$(git ls-remote --heads origin "$BRANCH")" ]]; then
        git worktree add --track -b "$BRANCH" "$WORKTREE" "origin/$BRANCH" >&2
    else
        git worktree add -b "$BRANCH" "$WORKTREE" origin/main >&2
        git -C "$WORKTREE" push -u origin "$BRANCH" >&2
    fi
fi

# Share the local-only paths listed in .worktreelinks rather than copying them or starting empty ones
if [[ -f "$WORKTREE/.worktreelinks" ]]; then
    while read -r path _; do
        if [[ -z "$path" || "$path" == \#* ]]; then
            continue
        fi
        if [[ ! -e "$MAIN_WORKTREE/$path" ]]; then
            echo "Skipping '$path' from .worktreelinks: not present in $MAIN_WORKTREE." >&2
            continue
        fi
        if ! git -C "$WORKTREE" check-ignore -q -- "$path"; then
            echo "'$path' is listed in .worktreelinks but is not gitignored; linking it would shadow tracked content." >&2
            exit 1
        fi
        if [[ ! -e "$WORKTREE/$path" ]]; then
            mkdir -p "$(dirname "$WORKTREE/$path")"
            ln -s "$MAIN_WORKTREE/$path" "$WORKTREE/$path"
        fi
    done < "$WORKTREE/.worktreelinks"
fi

echo "$WORKTREE"
