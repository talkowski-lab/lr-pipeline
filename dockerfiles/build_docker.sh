#!/usr/bin/env bash
set -euo pipefail

# Builds, tags and pushes a gnomad-lr Docker image to Artifact Registry.
# On main, auto-increments the image's kj_V<N> tag for the versioned history,
# and also (re)tags/pushes it as :latest, so any Dockerfile that inherits from
# it via "FROM .../<image>:latest" always picks up this build without edits.
#
# On a feature branch, the image is pushed as :<branch> only, so that branches
# built in parallel worktrees never overwrite each other's images or :latest.
# A branch build also inherits from the branch's own base image when one exists.
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

branch="$(git -C "${REPO_ROOT}" symbolic-ref --short HEAD)"

from_args=()
if [[ "${branch}" == "main" ]]; then
    current_version=$(gcloud artifacts docker tags list "${REGISTRY}/${image_name}" --format='value(tag)' 2>/dev/null \
        | grep -oE '^kj_V[0-9]+$' \
        | sed -E 's/kj_V//' \
        | sort -n \
        | tail -1) || true
    new_tag="kj_V$(( ${current_version:-0} + 1 ))"
    push_tags=("${new_tag}" latest)
else
    new_tag="${branch}"
    push_tags=("${new_tag}")

    # Inherit the branch's own base image, so a branch that changes what the base bakes in is testable
    base_image=$(grep -m1 '^FROM ' "${dockerfile}" | sed -E "s#^FROM ${REGISTRY}/([^:]+):latest\$#\1#") || true
    if [[ -n "${base_image}" && "${base_image}" != FROM* ]] \
        && gcloud artifacts docker tags list "${REGISTRY}/${base_image}" --format='value(tag)' 2>/dev/null \
            | grep -qxF "${branch}"; then
        from_args=(--from "${REGISTRY}/${base_image}:${branch}")
        echo "Inheriting ${REGISTRY}/${base_image}:${branch} instead of :latest"
    fi
fi

echo "Building ${image_name}:${new_tag} from ${dockerfile_name}"
podman build --platform linux/amd64 --network=host "${from_args[@]+"${from_args[@]}"}" "${build_args[@]+"${build_args[@]}"}" -f "${dockerfile}" -t "${image_name}:${new_tag}" "${REPO_ROOT}"

for tag in "${push_tags[@]}"; do
    podman tag "${image_name}:${new_tag}" "${REGISTRY}/${image_name}:${tag}"
    podman push "${REGISTRY}/${image_name}:${tag}"
done

echo "Pushed ${REGISTRY}/${image_name} with tags: ${push_tags[*]}"
