#!/usr/bin/env bash
# SessionStart/UserPromptSubmit hook: point sessions sitting in the shared main checkout at the worktree flow.
set -uo pipefail

cwd="$(jq -r '.cwd // empty')"
top="$(git -C "${cwd:-.}" rev-parse --show-toplevel 2>/dev/null)" || exit 0

if [[ "$top" == */.claude/worktrees/* ]]; then
    echo "This is a Claude-native worktree, which this repository does not support (branch naming, .dockstore.yml filters and merge_branch.sh all assume .github/scripts/new_worktree.sh). Stop and ask the user to restart from the main checkout."
    exit 0
fi

# A linked worktree has a .git file rather than a directory, and needs no reminder
[[ -d "$top/.git" ]] || exit 0

cat <<'EOF'
This session is in the shared main checkout of lr-pipeline, which stays on `main`, is never switched, and is read-only: a hook denies Edit and Write here. Before changing anything, run `.github/scripts/new_worktree.sh kj-<kebab-topic>` with a fresh topic name, then do every edit, command and check inside the path it prints. If it reports "Reusing existing worktree" and you are not resuming that task, pick another name. Read-only investigation and `.github/scripts/merge_branch.sh` are the only work that belongs in this checkout.
EOF
