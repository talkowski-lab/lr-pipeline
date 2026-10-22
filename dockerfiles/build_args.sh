#!/usr/bin/env bash
set -euo pipefail

# Prints one NAME=value line for every ARG declared without an inline default
# in dockerfiles/Dockerfile.<image-name>, resolved from dockerfiles/versions.env.
# Shared by build_docker.sh and the CI image build.
#
# Usage: dockerfiles/build_args.sh <image-name>

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

image_name="$1"
dockerfile="${REPO_ROOT}/dockerfiles/Dockerfile.${image_name}"

source "${REPO_ROOT}/dockerfiles/versions.env"

{ grep -oE '^ARG [A-Z0-9_]+$' "${dockerfile}" || true; } | awk '{print $2}' | while read -r arg_name; do
    version_key="${image_name}__${arg_name}"
    if [[ -z "${!version_key:-}" ]]; then
        echo "${version_key} is not set in dockerfiles/versions.env" >&2
        exit 1
    fi
    echo "${arg_name}=${!version_key}"
done
