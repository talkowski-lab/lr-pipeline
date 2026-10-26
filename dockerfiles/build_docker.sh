#!/usr/bin/env bash
set -euo pipefail

# Builds, tags and pushes a lr-pipeline Docker image to Artifact Registry.
# On main, tags the build <YYYYMMDD>-<short-sha> from the HEAD commit, so every
# tag maps to exactly one commit, and also (re)tags/pushes it as :latest, so any
# Dockerfile that inherits from it (ARG BASE_IMAGE=.../<image>:latest) always picks up
# this build without edits. A main build refuses a dirty tree for that reason.
#
# On a feature branch, the image is pushed as :<branch> only, so that branches
# built in parallel worktrees never overwrite each other's images or :latest.
# A branch build also inherits from the branch's own base image when one exists.
#
# Tool versions and upstream commits are centralized in dockerfiles/versions.env
# and resolved into --build-arg values by dockerfiles/build_args.sh.
#
# Every build carries OCI labels naming the source repository, the commit, the tag
# it was built as and the build time, and a push prints the pushed manifest digest.
#
# Dockerfiles are named Dockerfile.<image-name>, so the image
# name is the only argument needed. --no-push, in any position, builds and tags locally only.
#
# Usage: dockerfiles/build_docker.sh [--no-push] <image-name>
# Example: dockerfiles/build_docker.sh utils
# Set DOCKER_REGISTRY to push somewhere other than the lab registry

REGISTRY="${DOCKER_REGISTRY:-us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

push=true
positional=()
for arg in "$@"; do
    case "${arg}" in
        --no-push) push=false ;;
        *) positional+=("${arg}") ;;
    esac
done
if [[ ${#positional[@]} -ne 1 ]]; then
    echo "Usage: dockerfiles/build_docker.sh [--no-push] <image-name>" >&2
    exit 2
fi

image_name="${positional[0]}"
dockerfile_name="Dockerfile.${image_name}"
dockerfile="${REPO_ROOT}/dockerfiles/${dockerfile_name}"

resolved_args="$("${REPO_ROOT}/dockerfiles/build_args.sh" "${image_name}")"
build_args=()
while read -r build_arg; do
    [[ -n "${build_arg}" ]] && build_args+=(--build-arg "${build_arg}")
done <<< "${resolved_args}"

branch="$(git -C "${REPO_ROOT}" symbolic-ref --short HEAD)"

base_args=()
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
    base_image=$(grep -m1 '^ARG BASE_IMAGE=' "${dockerfile}" | sed -E "s#^ARG BASE_IMAGE=.*/([^/:]+):latest\$#\1#") || true
    if [[ -n "${base_image}" && "${base_image}" != ARG* ]] \
        && gcloud artifacts docker tags list "${REGISTRY}/${base_image}" --format='value(tag)' 2>/dev/null \
            | grep -qxF "${branch}"; then
        base_args=(--build-arg "BASE_IMAGE=${REGISTRY}/${base_image}:${branch}")
        echo "Inheriting ${REGISTRY}/${base_image}:${branch} instead of :latest"
    fi
fi

label_args=(
    --label "org.opencontainers.image.source=https://github.com/talkowski-lab/lr-pipeline"
    --label "org.opencontainers.image.revision=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
    --label "org.opencontainers.image.version=${new_tag}"
    --label "org.opencontainers.image.created=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
)

echo "Building ${image_name}:${new_tag} from ${dockerfile_name}"
podman build --platform linux/amd64 --network=host "${base_args[@]+"${base_args[@]}"}" "${build_args[@]+"${build_args[@]}"}" "${label_args[@]}" -f "${dockerfile}" -t "${image_name}:${new_tag}" "${REPO_ROOT}"

if [[ "${push}" == false ]]; then
    echo "Built localhost/${image_name}:${new_tag} without pushing"
    exit 0
fi

# Every tag points at the same manifest, so the digest of any one push names them all
digest_file="$(mktemp)"
trap 'rm -f "${digest_file}"' EXIT
for tag in "${push_tags[@]}"; do
    podman tag "${image_name}:${new_tag}" "${REGISTRY}/${image_name}:${tag}"
    podman push --digestfile "${digest_file}" "${REGISTRY}/${image_name}:${tag}"
done

echo "Pushed ${REGISTRY}/${image_name}@$(cat "${digest_file}") with tags: ${push_tags[*]}"
