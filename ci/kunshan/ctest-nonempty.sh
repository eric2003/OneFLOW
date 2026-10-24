#!/usr/bin/env bash

# Run a CTest subset with an explicit non-empty discovery gate.
# This prevents CTest's zero-exit "No tests were found!" result from being
# mistaken for validation success after PRE_TEST GoogleTest discovery.

set -euo pipefail

usage()
{
    echo "usage: $0 BUILD_DIR [ctest options ...]" >&2
    echo "example: $0 build -L hardware -R HIP --output-on-failure" >&2
    exit 2
}

build_dir="${1:-}"
if [[ -z "$build_dir" || "$build_dir" == -* ]]; then
    usage
fi
shift

ctest_bin="${CTEST_BIN:-ctest}"
if [[ ! -d "$build_dir" ]]; then
    echo "CTest build directory does not exist: $build_dir" >&2
    exit 2
fi

# PRE_TEST discovery is regenerated when the generated file is older than the
# executable, but stale files can survive metadata changes. Remove only the
# generated discovery files inside this build tree.
find "$build_dir" -type f -name '*_tests.cmake' -delete

list_log="$(mktemp "${TMPDIR:-/tmp}/oneflow-ctest-list.XXXXXX")"
run_log="$(mktemp "${TMPDIR:-/tmp}/oneflow-ctest-run.XXXXXX")"
trap 'rm -f "$list_log" "$run_log"' EXIT

"$ctest_bin" --test-dir "$build_dir" -N "$@" >"$list_log" 2>&1

if grep -Eq 'No tests were found|No tests were found!' "$list_log"; then
    cat "$list_log" >&2
    echo "CTest discovery produced an empty test set" >&2
    exit 1
fi

listed_count="$(
    sed -nE 's/^[[:space:]]*Test[[:space:]]*#[0-9]+:[[:space:]]+.*$/&/p' \
        "$list_log" | wc -l | tr -d ' '
)"
if [[ "$listed_count" -eq 0 ]]; then
    cat "$list_log" >&2
    echo "CTest discovery produced no numbered tests" >&2
    exit 1
fi

if grep -Eq '(^|[[:space:]])[^[:space:]]+_NOT_BUILT([[:space:]]|$)' "$list_log"; then
    cat "$list_log" >&2
    echo "CTest discovery contains NOT_BUILT entries" >&2
    exit 1
fi

"$ctest_bin" --test-dir "$build_dir" "$@" >"$run_log" 2>&1
cat "$run_log"

if ! grep -Eq '100% tests passed' "$run_log"; then
    echo "CTest did not report a fully passing suite" >&2
    exit 1
fi

passed_count="$(
    sed -nE 's/.*100% tests passed, 0 tests failed out of ([0-9]+).*/\1/p;
             s/.*100% tests passed, 0 tests failed out of ([0-9]+) tests.*/\1/p;
             s/.*100% tests passed.*out of ([0-9]+).*/\1/p' \
        "$run_log" | tail -1
)"
if [[ -z "$passed_count" || "$passed_count" -eq 0 ]]; then
    echo "CTest summary did not contain a positive test count" >&2
    exit 1
fi
