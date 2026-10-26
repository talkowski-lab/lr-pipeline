#!/usr/bin/env bash
set -euo pipefail

# Lists, for every image built from dockerfiles/Dockerfile.<image-name>, the
# manifest digest behind :latest and the <YYYYMMDD>-<short-sha> tag on the same
# manifest, which names the main commit it was built from. Exits non-zero if any
# image has no :latest, or a :latest with no dated tag, since then no commit can
# be recovered for it. Read-only; needs an authenticated gcloud.
#
# Usage: dockerfiles/audit_images.sh
# Set DOCKER_REGISTRY to audit somewhere other than the lab registry

REGISTRY="${DOCKER_REGISTRY:-us-central1-docker.pkg.dev/talkowski-sv-gnomad/kj-dockers}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

flagged=0
printf '%s\t%s\t%s\n' image latest_digest dated_tag
for dockerfile in "${REPO_ROOT}"/dockerfiles/Dockerfile.*; do
    image="${dockerfile##*/Dockerfile.}"
    row="$(gcloud artifacts docker images list "${REGISTRY}/${image}" --include-tags \
        --format='value(version,tags)' 2>/dev/null | awk -F'\t' '("," $2 ",") ~ /,latest,/' || true)"
    if [[ -z "${row}" ]]; then
        printf '%s\t-\t-\tno :latest\n' "${image}"
        flagged=1
        continue
    fi
    digest="${row%%$'\t'*}"
    dated="$(tr ',' '\n' <<< "${row#*$'\t'}" | grep -E '^[0-9]{8}-[0-9a-f]{7}$' | sort | tail -1 || true)"
    if [[ -z "${dated}" ]]; then
        printf '%s\t%s\t-\tno dated tag\n' "${image}" "${digest}"
        flagged=1
    else
        printf '%s\t%s\t%s\n' "${image}" "${digest}" "${dated}"
    fi
done

exit "${flagged}"
