#!/usr/bin/env bash
set -euo pipefail

# Builds, tags and pushes a gnomad-lr Docker image to Artifact Registry.
# Auto-increments the image's kj_V<N> tag for the versioned history, and
# also (re)tags/pushes it as :latest, so any Dockerfile that inherits from
# it via "FROM .../<image>:latest" always picks up this build without edits.
#
# Tool/library versions are centralized in dockerfiles/versions.env, keyed
# as <image-name>__<ARG_NAME>. Every ARG declared in the target Dockerfile
# (with no inline default) is resolved from there and passed as --build-arg.
#
# Dockerfiles are named Dockerfile.<image-name>, so the image
# name is the only argument needed.
#
# Usage: dockerfiles/build_docker.sh <image-name>
# Example: dockerfiles/build_docker.sh utils

REGISTRY="us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

image_name="$1"
dockerfile_name="Dockerfile.${image_name}"
dockerfile="${REPO_ROOT}/dockerfiles/${dockerfile_name}"

source "${REPO_ROOT}/dockerfiles/versions.env"

build_args=()
for arg_name in $(grep -oE '^ARG [A-Z_]+' "${dockerfile}" | awk '{print $2}'); do
    version_key="${image_name}__${arg_name}"
    build_args+=(--build-arg "${arg_name}=${!version_key}")
done

current_version=$(gcloud artifacts docker tags list "${REGISTRY}/${image_name}" --format='value(tag)' 2>/dev/null \
    | grep -oE '^kj_V[0-9]+$' \
    | sed -E 's/kj_V//' \
    | sort -n \
    | tail -1) || true
new_tag="kj_V$(( ${current_version:-0} + 1 ))"

echo "Building ${image_name}:${new_tag} from ${dockerfile_name}"
podman build --platform linux/amd64 --network=host "${build_args[@]+"${build_args[@]}"}" -f "${dockerfile}" -t "${image_name}:${new_tag}" "${REPO_ROOT}"

for tag in "${new_tag}" latest; do
    podman tag "${image_name}:${new_tag}" "${REGISTRY}/${image_name}:${tag}"
    podman push "${REGISTRY}/${image_name}:${tag}"
done

echo "Pushed ${REGISTRY}/${image_name}:${new_tag} and :latest"
