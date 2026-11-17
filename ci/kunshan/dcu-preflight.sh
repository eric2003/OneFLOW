#!/bin/bash
# Validate the minimum device/toolchain contract before a standalone HIP build.

set -u

if [ "$#" -ne 1 ]; then
    echo "usage: $0 ARTIFACT_FILE" >&2
    exit 2
fi

artifact_file="$1"
artifact_dir=$(dirname "$artifact_file")
mkdir -p "$artifact_dir"

expected_arch="${ONEFLOW_HIP_ARCHITECTURES:-gfx906}"
amd_comgr_dir="${AMD_COMGR_DIR:-}"
rc=0

{
    echo "expected_arch=$expected_arch"
    echo "amd_comgr_dir=${amd_comgr_dir:-unset}"
    echo "rocminfo=$(command -v rocminfo 2>/dev/null || echo missing)"
} > "$artifact_file"

if ! command -v rocminfo >/dev/null 2>&1; then
    echo "DCU_PREFLIGHT_FAIL reason=rocminfo-missing" | tee -a "$artifact_file" >&2
    rc=1
else
    rocminfo_output="${artifact_file}.rocminfo"
    if ! rocminfo > "$rocminfo_output" 2>&1; then
        echo "DCU_PREFLIGHT_FAIL reason=rocminfo-failed" | tee -a "$artifact_file" >&2
        rc=1
    elif ! grep -Eq \
        "^[[:space:]]*Name:[[:space:]]*${expected_arch}[[:space:]]*$" \
        "$rocminfo_output"
    then
        echo "DCU_PREFLIGHT_FAIL reason=target-architecture-not-found" \
            | tee -a "$artifact_file" >&2
        rc=1
    else
        echo "target_architecture=$expected_arch" >> "$artifact_file"
    fi
fi

if [ -z "$amd_comgr_dir" ] || [ ! -d "$amd_comgr_dir" ]; then
    echo "DCU_PREFLIGHT_FAIL reason=amd_comgr-directory-missing" \
        | tee -a "$artifact_file" >&2
    rc=1
elif ! find "$amd_comgr_dir" -maxdepth 1 -type f \
    \( -name 'amd_comgr-config.cmake' -o -name 'amd_comgrConfig.cmake' \) \
    -print -quit | grep -q .
then
    echo "DCU_PREFLIGHT_FAIL reason=amd_comgr-config-missing" \
        | tee -a "$artifact_file" >&2
    rc=1
else
    echo "amd_comgr_config=present" >> "$artifact_file"
fi

if [ "$rc" -eq 0 ]; then
    echo "DCU_PREFLIGHT_PASS expected_arch=$expected_arch" \
        | tee -a "$artifact_file"
fi
exit "$rc"
