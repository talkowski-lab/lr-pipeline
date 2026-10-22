#!/usr/bin/env bash
set -euo pipefail

# Wraps podman for miniwdl, which pulls the image before every task; skip the pull when the image is already local
if [[ "${1:-}" == "pull" ]] && podman image exists "${2}"; then
    exit 0
fi
exec podman "$@"
