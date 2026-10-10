#!/usr/bin/env bash
set -euo pipefail

sct_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
sct_project_root="$(dirname "${sct_root}")"
sct_container='sctenifoldknk-loong64-sdaa-adapt-smoke'
sct_image='sctenifoldknk-loong64:r2.1.0-cpu-adapt-sdaa'
sct_expected_base_id='sha256:3e425da47a74a03b589f8e1c18a67dcb90c8bcbde977cc1df09fa01209b7a524'
sct_output_dir="${sct_project_root}/work/loong64-smoke-test-adapt-sdaa"

export HOST_UID="$(id -u)"
export HOST_GID="$(id -g)"
export BASE_IMAGE='lasermpnn-sdaa321-loong:torch2.12.0-sdaa3.2.1-py312'
export BASE_IMAGE_ID="${sct_expected_base_id}"

test "$(docker image inspect "${BASE_IMAGE}" --format '{{.Id}}')" = \
  "${sct_expected_base_id}"
test "$(docker image inspect "${sct_image}" --format '{{.Architecture}}')" = 'loong64'
sct_source_commit='08118337f6fc22cca5a22c443ccb78ef920dfa3b'
git -C "${sct_root}" cat-file -e "${sct_source_commit}^{commit}"
git -C "${sct_root}" diff --quiet "${sct_source_commit}" HEAD -- \
  . ':(exclude)sdaa' ':(exclude).gitignore'
test -z "$(GIT_OPTIONAL_LOCKS=0 git -C "${sct_root}" status --porcelain --untracked-files=all | awk '$2 !~ /^sdaa\// && $2 != ".gitignore"')"

mkdir -p "${sct_output_dir}"
if docker container inspect "${sct_container}" >/dev/null 2>&1; then
  sct_expected_image_id="$(docker image inspect "${sct_image}" --format '{{.Id}}')"
  sct_container_image_id="$(docker inspect "${sct_container}" --format '{{.Image}}')"
  if [[ "${sct_container_image_id}" != "${sct_expected_image_id}" ]]; then
    printf 'Existing container uses %s, expected %s; refusing to replace it.\n' \
      "${sct_container_image_id}" "${sct_expected_image_id}" >&2
    exit 1
  fi
  docker start "${sct_container}" >/dev/null
else
  docker run --detach \
    --name "${sct_container}" \
    --init \
    --user "${HOST_UID}:${HOST_GID}" \
    --workdir /workspace/repro \
    --network none \
    --read-only \
    --tmpfs /tmp:rw,nosuid,nodev,size=2g \
    --cap-drop ALL \
    --security-opt no-new-privileges:true \
    --pids-limit 2048 \
    --env HOME=/tmp \
    --env SCTENIFOLDKNK_OUTPUT_DIR=/workspace/output \
    --env OMP_NUM_THREADS=1 \
    --env OPENBLAS_NUM_THREADS=1 \
    --env MKL_NUM_THREADS=1 \
    --env TZ=UTC \
    --env LANG=C.UTF-8 \
    --env LC_ALL=C.UTF-8 \
    --mount "type=bind,source=${sct_project_root},target=/workspace/repro,readonly" \
    --mount "type=bind,source=${sct_output_dir},target=/workspace/output" \
    "${sct_image}" sleep infinity >/dev/null
fi

sct_devices="$(docker inspect "${sct_container}" --format '{{json .HostConfig.Devices}}')"
sct_device_requests="$(docker inspect "${sct_container}" --format '{{json .HostConfig.DeviceRequests}}')"
test "${sct_devices}" = '[]'
if [[ "${sct_device_requests}" != 'null' && "${sct_device_requests}" != '[]' ]]; then
  printf 'Unexpected device requests: %s\n' "${sct_device_requests}" >&2
  exit 1
fi

test "$(docker inspect "${sct_container}" --format '{{.HostConfig.NetworkMode}}')" = 'none'
test "$(docker inspect "${sct_container}" --format '{{.HostConfig.ReadonlyRootfs}}')" = 'true'

docker exec -i "${sct_container}" bash -lc '
  set -euo pipefail
  test "${VIRTUAL_ENV}" = "/opt/venvs/scTenifoldKnk"
  test "$(python3 -c '\''import sys; print(sys.prefix)'\'')" = "/opt/venvs/scTenifoldKnk"
  python3 -c '\''import importlib.util; assert importlib.util.find_spec("torch") is None'\''
  if compgen -G "/dev/tcaicard*" >/dev/null; then
    echo "Unexpected SDAA device in CPU baseline container" >&2
    exit 1
  fi
  Rscript /workspace/repro/scTenifoldKnk/sdaa/scripts/smoke-test.R
'

stat -c '%u:%g %a %n' "${sct_output_dir}"/*
test "$(stat -c '%u:%g' "${sct_output_dir}/result.rds")" = "${HOST_UID}:${HOST_GID}"
sha256sum \
  "${sct_output_dir}/result.rds" \
  "${sct_output_dir}/diff_regulation.csv" \
  > "${sct_output_dir}/output-checksums.sha256"
printf 'devices=%s\ndevice_requests=%s\ngpu_required=false\n' \
  "${sct_devices}" "${sct_device_requests}" \
  > "${sct_output_dir}/device-audit.txt"

printf 'Container smoke test passed.\nOutput: %s\n' "${sct_output_dir}"
