#!/usr/bin/env bash
set -euo pipefail

sct_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
sct_base_image='lasermpnn-sdaa321-loong:torch2.12.0-sdaa3.2.1-py312'
sct_expected_base_id='sha256:3e425da47a74a03b589f8e1c18a67dcb90c8bcbde977cc1df09fa01209b7a524'
sct_output_image='sctenifoldknk-loong64:r2.1.0-cpu-adapt-sdaa'
sct_actual_base_id="$(docker image inspect "${sct_base_image}" --format '{{.Id}}')"

if [[ "${sct_actual_base_id}" != "${sct_expected_base_id}" ]]; then
  printf 'Base image mismatch: expected %s, got %s\n' \
    "${sct_expected_base_id}" "${sct_actual_base_id}" >&2
  exit 1
fi

test "$(uname -m)" = 'loongarch64'
sct_source_commit='08118337f6fc22cca5a22c443ccb78ef920dfa3b'
git -C "${sct_root}" cat-file -e "${sct_source_commit}^{commit}"
git -C "${sct_root}" diff --quiet "${sct_source_commit}" HEAD -- \
  . ':(exclude)sdaa' ':(exclude).gitignore'
test -z "$(GIT_OPTIONAL_LOCKS=0 git -C "${sct_root}" status --porcelain --untracked-files=all | awk '$2 !~ /^sdaa\// && $2 != ".gitignore"')"
(
  cd "${sct_root}/sdaa/vendor"
  sha256sum --check checksums.sha256
)

docker build --pull=false \
  --file "${sct_root}/sdaa/Dockerfile" \
  --tag "${sct_output_image}" \
  --build-arg "BASE_IMAGE=${sct_base_image}" \
  --build-arg "BASE_IMAGE_ID=${sct_expected_base_id}" \
  --build-arg "HOST_UID=$(id -u)" \
  --build-arg "HOST_GID=$(id -g)" \
  "${sct_root}"

docker image inspect "${sct_output_image}" --format \
  'image={{.Id}} arch={{.Architecture}} user={{json .Config.User}} base_id={{index .Config.Labels "org.sctenifoldknk.base.image-id"}}'
