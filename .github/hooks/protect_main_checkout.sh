#!/usr/bin/env bash
# PreToolUse hook: deny edits to tracked files in the shared main checkout, which stays on main.
# Handles both payload shapes: a named file path (Claude Code's Edit/Write) and a patch whose
# only location is the session's working directory (Codex's apply_patch).
set -uo pipefail

deny() {
    jq -n --arg reason "$1" '{hookSpecificOutput: {hookEventName: "PreToolUse",
        permissionDecision: "deny", permissionDecisionReason: $reason}}'
    exit 0
}

REASON="The main checkout stays on main and is read-only; run .github/scripts/new_worktree.sh kj-<kebab-topic> and edit inside the worktree it prints."

input="$(cat)"
path="$(jq -r '.tool_input.file_path // .tool_input.notebook_path // empty' <<<"$input")"

if [[ -z "$path" ]]; then
    # No path in the payload, so judge the session's working directory instead
    cwd="$(jq -r '.cwd // empty' <<<"$input")"
    top="$(git -C "${cwd:-.}" rev-parse --show-toplevel 2>/dev/null)" || exit 0
    [[ -d "$top/.git" ]] || exit 0
    deny "$REASON"
fi

# Walk up to the nearest existing directory, since the target file may not exist yet
dir="$(dirname "$path")"
while [[ ! -d "$dir" && "$dir" != "/" ]]; do
    dir="$(dirname "$dir")"
done

top="$(git -C "$dir" rev-parse --show-toplevel 2>/dev/null)" || exit 0
[[ -d "$top/.git" ]] || exit 0

# Gitignored local files, such as the shared data/ directory, stay editable. Exit 1 means tracked;
# anything else means git could not judge the path, in which case the edit is allowed through
git -C "$top" check-ignore -q -- "$path" 2>/dev/null
case $? in
    1) ;;
    *) exit 0 ;;
esac

deny "$REASON"
