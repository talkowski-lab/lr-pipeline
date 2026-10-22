#!/usr/bin/env bash
set -euo pipefail

# Builds, tags and pushes a gnomad-lr Docker image to Artifact Registry.
# On main, tags the build <YYYYMMDD>-<short-sha> from the HEAD commit, so every
# tag maps to exactly one commit, and also (re)tags/pushes it as :latest, so any
# Dockerfile that inherits from it via "FROM .../<image>:latest" always picks up
# this build without edits. A main build refuses a dirty tree for that reason.
#
# On a feature branch, the image is pushed as :<branch> only, so that branches
# built in parallel worktrees never overwrite each other's images or :latest.
# A branch build also inherits from the branch's own base image when one exists.
#
# Tool versions and upstream commits are centralized in dockerfiles/versions.env
# and resolved into --build-arg values by dockerfiles/build_args.sh.
#
# Dockerfiles are named Dockerfile.<image-name>, so the image
# name is the only argument needed. --no-push builds and tags locally only.
#
# Usage: dockerfiles/build_docker.sh [--no-push] <image-name>
# Example: dockerfiles/build_docker.sh utils

REGISTRY="us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

push=true
if [[ "${1:-}" == "--no-push" ]]; then
    push=false
    shift
fi

image_name="$1"
dockerfile_name="Dockerfile.${image_name}"
dockerfile="${REPO_ROOT}/dockerfiles/${dockerfile_name}"

resolved_args="$("${REPO_ROOT}/dockerfiles/build_args.sh" "${image_name}")"
build_args=()
while read -r build_arg; do
    [[ -n "${build_arg}" ]] && build_args+=(--build-arg "${build_arg}")
done <<< "${resolved_args}"

branch="$(git -C "${REPO_ROOT}" symbolic-ref --short HEAD)"

from_args=()
if [[ "${branch}" == "main" ]]; then
    if [[ "${push}" == true && -n "$(git -C "${REPO_ROOT}" status --porcelain)" ]]; then
        echo "Refusing to push a main build from a dirty tree, since its tag names the HEAD commit." >&2
        exit 1
    fi
    new_tag="$(git -C "${REPO_ROOT}" log -1 --format=%cd --date=format:%Y%m%d)-$(git -C "${REPO_ROOT}" rev-parse --short=7 HEAD)"
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

if [[ "${push}" == false ]]; then
    echo "Built localhost/${image_name}:${new_tag} without pushing"
    exit 0
fi

for tag in "${push_tags[@]}"; do
    podman tag "${image_name}:${new_tag}" "${REGISTRY}/${image_name}:${tag}"
    podman push "${REGISTRY}/${image_name}:${tag}"
done

echo "Pushed ${REGISTRY}/${image_name} with tags: ${push_tags[*]}"
