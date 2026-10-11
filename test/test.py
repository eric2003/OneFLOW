#-*- coding:utf-8 -*-
'''
---------------------------------------------------------------------------
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------------
'''
import argparse
import json
import os
import sys
import logging
import math
import time
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import List, Optional

from residual_db import compare_case

# NOTE: default level is INFO so routine debug noise (variable dumps, argv
# echoes, etc.) is hidden by default. Set ONEFLOW_TEST_LOG_LEVEL=DEBUG to see
# everything that used to be printed unconditionally in the old script.
logging.basicConfig(
    level=os.environ.get("ONEFLOW_TEST_LOG_LEVEL", "INFO"),
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger("oneflow.test")


def list_directory_contents(file_dir: str) -> None:
    # Was `my_dir_cmd`: pure debug helper, kept but demoted to DEBUG level
    # instead of always printing to stdout.
    for file in os.listdir(file_dir):
        logger.debug(file)


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


FLOAT_ABS_TOLERANCE = 1.0e-8


def cmp_float_str(str1: str, str2: str) -> bool:
    """Compare two whitespace-separated lines token by token.

    Numeric tokens are compared as floats within FLOAT_ABS_TOLERANCE;
    non-numeric tokens are assumed equal already (caller only calls this
    after a plain str1 != str2 check).

    BUGFIX (step 2): the original version indexed tokens2[i] using only
    len(tokens1), so a row with a different number of columns raised an
    uncaught IndexError and crashed the whole test suite instead of being
    reported as a normal comparison failure.
    """
    tokens1 = str1.split()
    tokens2 = str2.split()
    if len(tokens1) != len(tokens2):
        return False
    for token1, token2 in zip(tokens1, tokens2):
        # Labels and integer/topology identifiers are semantic data, not
        # approximate numbers. They must match exactly.
        if not is_number(token1) or not is_number(token2):
            if token1 != token2:
                return False
            continue

        value1 = float(token1)
        value2 = float(token2)
        # NaN/Inf must never make a malformed result look equivalent because
        # comparisons involving NaN are false.
        if not math.isfinite(value1) or not math.isfinite(value2):
            if token1 != token2:
                return False
            # Even identical non-finite values are invalid numerical output.
            return False
        if abs(value1 - value2) > FLOAT_ABS_TOLERANCE:
            return False
    return True


MAX_REPORTED_LINE_DIFFS = 20


def cmp_file(
    file_name1: str,
    file_name2: str,
    max_reported_diffs: int = MAX_REPORTED_LINE_DIFFS,
    mode: str = "text",
) -> bool:
    """Compare text results numerically by line or binary artifacts byte-for-byte.

    Binary paths must be listed in the case's optional autotest/binary.txt
    manifest. Byte comparison is suitable only for deterministic binary
    formats; grid tests should additionally validate geometry and topology.
    """
    if mode == "binary":
        try:
            with open(file_name1, "rb") as f1, open(file_name2, "rb") as f2:
                expected = f1.read()
                actual = f2.read()
        except OSError as error:
            logger.info("Could not compare binary files %s vs %s: %s", file_name1, file_name2, error)
            return False
        if expected == actual:
            return True
        mismatch = next(
            (index for index, (left, right) in enumerate(zip(expected, actual)) if left != right),
            min(len(expected), len(actual)),
        )
        logger.info(
            "Binary mismatch: %s vs %s (expected_bytes=%d actual_bytes=%d first_difference_offset=%d)",
            file_name1, file_name2, len(expected), len(actual), mismatch,
        )
        return False

    differences: List[str] = []
    try:
        with open(file_name1, "r", encoding="utf-8-sig") as f1, \
                open(file_name2, "r", encoding="utf-8-sig") as f2:
            line_id = 0
            while True:
                str1 = f1.readline()
                str2 = f2.readline()
                if not str1 and not str2:
                    break
                if not str1 or not str2:
                    differences.append(
                        f"line {line_id}: different end-of-file position"
                    )
                    break
                if str1 != str2 and not cmp_float_str(str1, str2):
                    differences.append(
                        f"line {line_id}: line1={str1.rstrip()!r} "
                        f"line2={str2.rstrip()!r}"
                    )
                line_id += 1
    except OSError as error:
        logger.info("Could not compare %s vs %s: %s", file_name1, file_name2, error)
        return False

    if differences:
        logger.info("Differences found: %s vs %s", file_name1, file_name2)
        for diff in differences[:max_reported_diffs]:
            logger.info("  %s", diff)
        remaining = len(differences) - max_reported_diffs
        if remaining > 0:
            logger.info("  ... and %s more differing line(s)", remaining)
    return not differences


def get_filename_list(filename: str) -> List[str]:
    """Read one relative result-file path per line from a manifest."""
    with open(filename, "r", encoding="utf-8-sig") as f:
        return [line.strip() for line in f.readlines() if line.strip()]


def get_optional_filename_list(filename: str) -> List[str]:
    """Read an optional path manifest, returning an empty list when absent."""
    try:
        return get_filename_list(filename)
    except FileNotFoundError:
        return []


def _validate_case_manifests(
    test_project_dir: str,
    text_names: List[str],
    binary_names: List[str],
    required_names: List[str],
) -> List[str]:
    """Validate output manifests before running a solver.

    Manifest entries are case-relative paths. Comparison baselines are stored
    by output basename, so two compared paths with the same basename would
    otherwise silently share one reference file.
    """
    errors: List[str] = []
    case_root = os.path.abspath(test_project_dir)
    compared_basenames = {}

    for manifest_name, names, compared in (
        ("compare.txt", text_names, True),
        ("binary.txt", binary_names, True),
        ("required.txt", required_names, False),
    ):
        for relative_name in names:
            normalized = relative_name.replace("\\", "/")
            parts = normalized.split("/")
            if (
                not normalized
                or os.path.isabs(normalized)
                or any(part in {"", ".", ".."} for part in parts)
            ):
                errors.append(
                    f"{manifest_name}: invalid case-relative path {relative_name!r}"
                )
                continue

            resolved = os.path.abspath(os.path.join(case_root, normalized))
            if os.path.commonpath([case_root, resolved]) != case_root:
                errors.append(
                    f"{manifest_name}: path escapes case directory {relative_name!r}"
                )
                continue

            if compared:
                basename = os.path.basename(normalized).casefold()
                previous = compared_basenames.get(basename)
                if previous is not None:
                    errors.append(
                        f"{manifest_name}: {relative_name!r} and {previous!r} "
                        "share a baseline basename"
                    )
                else:
                    compared_basenames[basename] = relative_name

    if not text_names and not binary_names and not required_names:
        errors.append(
            "no output assertions configured; a successful process exit alone "
            "cannot validate a test case"
        )
    return errors


def _required_output_path(test_project_dir: str, relative_name: str) -> str:
    """Resolve a required output path without allowing paths outside the case."""
    normalized = relative_name.replace("\\", "/")
    parts = normalized.split("/")
    if not normalized or os.path.isabs(normalized) or any(part == ".." for part in parts):
        raise ValueError(f"Required output must be a case-relative path: {relative_name!r}")
    case_root = os.path.abspath(test_project_dir)
    output_path = os.path.abspath(os.path.join(case_root, normalized))
    if os.path.commonpath([case_root, output_path]) != case_root:
        raise ValueError(f"Required output escapes case directory: {relative_name!r}")
    return output_path


def _remove_previous_required_outputs(test_project_dir: str, relative_names: List[str]) -> None:
    """Remove old outputs before execution so stale files cannot satisfy a test."""
    for relative_name in relative_names:
        output_path = _required_output_path(test_project_dir, relative_name)
        if os.path.isfile(output_path) or os.path.islink(output_path):
            os.remove(output_path)


def _check_required_outputs(test_project_dir: str, relative_names: List[str]) -> bool:
    """Require each declared output to be a non-empty regular file."""
    all_ok = True
    for relative_name in relative_names:
        output_path = _required_output_path(test_project_dir, relative_name)
        if not os.path.isfile(output_path):
            logger.info("Required output is missing: %s", output_path)
            all_ok = False
        elif os.path.getsize(output_path) == 0:
            logger.info("Required output is empty: %s", output_path)
            all_ok = False
    return all_ok


@dataclass
class ComparisonPair:
    """One (baseline, generated) result-file pair to compare for a test case."""

    baseline_path: str
    generated_path: str
    mode: str = "text"


def build_comparison_pairs(
    test_project_dir: str,
    relative_names: List[str],
    binary_names: Optional[List[str]] = None,
) -> List[ComparisonPair]:
    """Build comparison pairs from text and optional binary file manifests."""
    pairs = []
    file_specs = [(name, "text") for name in relative_names]
    file_specs.extend((name, "binary") for name in (binary_names or []))
    for relative_name, mode in file_specs:
        base_name = relative_name.replace("\\", "/").split("/")[-1]
        baseline_file = os.path.join(test_project_dir, "autotest", base_name)
        generated_file = os.path.join(test_project_dir, relative_name)
        pairs.append(
            ComparisonPair(
                baseline_path=os.path.normpath(os.path.abspath(baseline_file)),
                generated_path=os.path.normpath(os.path.abspath(generated_file)),
                mode=mode,
            )
        )
    return pairs


@dataclass
class TestResult:
    """Outcome of one case, including its path for unambiguous reporting."""

    case_name: str
    case_dir: str
    passed: bool
    return_code: Optional[int] = None
    files_ok: Optional[bool] = None
    residual_ok: Optional[bool] = None
    required_outputs_ok: Optional[bool] = None
    failure_reasons: List[str] = field(default_factory=list)

    @property
    def failure_summary(self) -> str:
        return "; ".join(self.failure_reasons) if self.failure_reasons else ""


def _run_solver_process(cmd: str, process_timeout: Optional[float]) -> Optional[int]:
    """Run the solver command and return its exit code (None on timeout).

    BUGFIX (step 2): replaces the original Popen + busy-wait polling loop
    (`while poll() is None: sleep(0.5)`), which had no timeout and could
    hang the whole suite forever if a solver process got stuck.
    """
    try:
        completed = subprocess.run(cmd, shell=True, timeout=process_timeout)
        return completed.returncode
    except subprocess.TimeoutExpired:
        logger.info(
            "Test process timed out after %s second(s): %s", process_timeout, cmd
        )
        return None


def _compare_result_files(comparison_pairs: List[ComparisonPair]) -> bool:
    """Compare every (baseline, generated) pair for one case.

    REFACTOR (step 4): checks every pair (instead of stopping at the first
    failing file, as the original did) so a failing case's log shows *all*
    of its mismatching files at once, not just the first.
    """
    all_ok = True
    for pair in comparison_pairs:
        if not cmp_file(pair.baseline_path, pair.generated_path, mode=pair.mode):
            all_ok = False
    return all_ok


def _compare_residuals(
    test_project_dir: str,
    residual_db_path: Optional[str],
    residual_tolerance: Optional[float],
) -> Optional[bool]:
    """Compare residual files against the JSON baseline DB, if configured.

    Returns None when residual comparison isn't enabled for this run
    (no residual_db_path given), otherwise True/False for pass/fail.
    """
    if not residual_db_path:
        return None

    # Residual baselines are optional per case. Grid conversion cases validate
    # generated artifacts and required outputs but do not produce CFD residuals.
    try:
        with open(residual_db_path, "r", encoding="utf-8") as stream:
            database = json.load(stream)
    except (OSError, ValueError):
        residual_result = compare_case(
            test_project_dir, residual_db_path, absolute_tolerance=residual_tolerance
        )
        logger.info("residual_baseline=%s", residual_result)
        return residual_result["ok"]

    case_name = os.path.basename(os.path.normpath(test_project_dir))
    if case_name not in database.get("cases", {}):
        normalized_case_dir = os.path.normpath(test_project_dir)
        category = normalized_case_dir.split(os.sep, 1)[0]
        if category != "cfd":
            logger.info("residual_baseline=SKIPPED (no baseline applies to case=%s)", case_name)
            return None
        # CFD cases remain strict: a missing baseline entry is a test failure.

    residual_result = compare_case(
        test_project_dir, database, absolute_tolerance=residual_tolerance
    )
    logger.info("residual_baseline=%s", residual_result)
    return residual_result["ok"]


def run_test(
    test_project_dir: str,
    mpi_cmd: str,
    exe_cmd: str,
    residual_db_path: Optional[str] = None,
    residual_tolerance: Optional[float] = None,
    process_timeout: Optional[float] = None,
) -> TestResult:
    """Run one test case end to end and report exactly what passed/failed.

    REFACTOR (step 4): previously this one function did path building,
    process launching, file comparison, residual comparison, and result
    aggregation all inline, and returned a bare bool. It's now split into
    small single-purpose helpers (_run_solver_process /
    _compare_result_files / _compare_residuals) and returns a TestResult
    that records *why* a case failed, not just that it did.
    """
    case_name = os.path.basename(os.path.normpath(test_project_dir))
    start_time = time.perf_counter()
    logger.debug("test_project_dir=%s", test_project_dir)

    test_script = test_project_dir + "/autotest/compare.txt"
    abs_test_script = os.path.normpath(os.path.abspath(test_script))
    relative_names = get_filename_list(abs_test_script)
    binary_script = os.path.join(test_project_dir, "autotest", "binary.txt")
    binary_names = get_optional_filename_list(binary_script)
    comparison_pairs = build_comparison_pairs(test_project_dir, relative_names, binary_names)
    required_script = os.path.join(test_project_dir, "autotest", "required.txt")
    required_names = get_optional_filename_list(required_script)

    manifest_errors = _validate_case_manifests(
        test_project_dir, relative_names, binary_names, required_names
    )
    if manifest_errors:
        result = TestResult(
            case_name=case_name,
            case_dir=os.path.normpath(test_project_dir),
            passed=False,
            files_ok=False,
            failure_reasons=manifest_errors,
        )
        logger.error(
            "case=%s manifest validation failed: %s",
            case_name, "; ".join(manifest_errors),
        )
        return result

    _remove_previous_required_outputs(test_project_dir, required_names)

    cmd = mpi_cmd + " " + exe_cmd + " 0 " + test_project_dir
    logger.debug("cmd=%s", cmd)

    return_code = _run_solver_process(cmd, process_timeout)
    files_ok = _compare_result_files(comparison_pairs)
    residual_ok = _compare_residuals(test_project_dir, residual_db_path, residual_tolerance)
    required_outputs_ok = (
        _check_required_outputs(test_project_dir, required_names)
        if required_names else None
    )

    failure_reasons = []
    if return_code != 0:
        failure_reasons.append(f"solver exit code={return_code}")
    if not files_ok:
        failure_reasons.append("result file mismatch")
    if residual_ok is False:
        failure_reasons.append("residual baseline mismatch")
    if required_outputs_ok is False:
        failure_reasons.append("required output missing or empty")

    result = TestResult(
        case_name=case_name,
        case_dir=os.path.normpath(test_project_dir),
        passed=not failure_reasons,
        return_code=return_code,
        files_ok=files_ok,
        residual_ok=residual_ok,
        required_outputs_ok=required_outputs_ok,
        failure_reasons=failure_reasons,
    )
    elapsed_time = time.perf_counter() - start_time

    logger.info(
        "case=%s return_code=%s files_ok=%s residual_ok=%s "
        "required_outputs_ok=%s passed=%s elapsed=%.2f s",
        case_name,
        return_code,
        files_ok,
        residual_ok,
        required_outputs_ok,
        result.passed,
        elapsed_time,
    )    
    return result


def read_suite_file(filename: str) -> List[str]:
    """Expand a suite file into an ordered list of case directories.

    A line beginning with '@' includes another suite file, resolved relative
    to the including suite. This lets focused suites compose into CI suites
    without duplicating case lists.
    """
    def expand_suite(suite_filename: str, include_stack: List[str]) -> List[str]:
        suite_path = os.path.abspath(suite_filename)
        if suite_path in include_stack:
            cycle = include_stack + [suite_path]
            raise ValueError(
                "Suite include cycle: " + " -> ".join(cycle)
            )

        try:
            with open(suite_path, "r", encoding="utf-8-sig") as suite_file:
                lines = suite_file.readlines()
        except OSError as error:
            raise OSError(
                f"Cannot read suite file {suite_path}: {error}"
            ) from error

        case_dirs: List[str] = []
        for line_number, raw_line in enumerate(lines, start=1):
            entry = raw_line.strip()
            if not entry or entry.startswith("#"):
                continue

            if entry.startswith("@"):
                include_name = entry[1:].strip()
                if not include_name:
                    raise ValueError(
                        f"Empty suite include at {suite_path}:{line_number}"
                    )
                include_path = (
                    include_name
                    if os.path.isabs(include_name)
                    else os.path.join(os.path.dirname(suite_path), include_name)
                )
                case_dirs.extend(
                    expand_suite(include_path, include_stack + [suite_path])
                )
            else:
                case_dirs.append(entry)

        return case_dirs

    case_dirs = expand_suite(filename, [])
    seen = set()
    for case_dir in case_dirs:
        normalized = os.path.normcase(
            os.path.normpath(case_dir.replace("\\", os.sep).replace("/", os.sep))
        )
        if normalized in seen:
            raise ValueError(
                f"Duplicate case in expanded suite {filename}: {case_dir}"
            )
        seen.add(normalized)

    return case_dirs


def run_all_test(
    case_dirs: List[str],
    mpi_cmd: str,
    exe_cmd: str,
    residual_db_path: Optional[str] = None,
    residual_tolerance: Optional[float] = None,
    process_timeout: Optional[float] = None,
) -> List[TestResult]:
    """Run every case directory in `case_dirs` and collect their results.

    REFACTOR (step 4): now takes a plain list of case directories instead
    of a suite *filename*, so both the real suite-file path and the
    --dry-run fake cases can share this same function. Also takes
    mpi_cmd/exe_cmd explicitly instead of run_test() reaching into
    sys.argv itself, which made run_test impossible to call from anywhere
    other than main().
    """
    results = []
    for case_dir in case_dirs:
        logger.info("project_dir=%s", case_dir)
        results.append(
            run_test(
                case_dir, mpi_cmd, exe_cmd, residual_db_path, residual_tolerance, process_timeout
            )
        )
    return results


def _normalized_case_path(path: str) -> str:
    return os.path.normcase(os.path.normpath(path.replace("\\", os.sep).replace("/", os.sep)))


def read_suite_groups(filename: str):
    """Return named suite memberships for reporting, including included suites."""
    groups = []
    visited = set()

    def visit(suite_filename: str) -> None:
        suite_path = os.path.abspath(suite_filename)
        normalized_suite = os.path.normcase(os.path.normpath(suite_path))
        if normalized_suite in visited:
            return
        visited.add(normalized_suite)
        with open(suite_path, "r", encoding="utf-8-sig") as stream:
            entries = [line.strip() for line in stream if line.strip() and not line.strip().startswith("#")]
        case_dirs = read_suite_file(suite_path)
        groups.append((os.path.basename(suite_path), case_dirs))
        for entry in entries:
            if entry.startswith("@") and entry[1:].strip():
                include_name = entry[1:].strip()
                include_path = include_name if os.path.isabs(include_name) else os.path.join(os.path.dirname(suite_path), include_name)
                visit(include_path)

    visit(filename)
    return groups


def _summarize_results(results: List[TestResult]) -> str:
    passed = sum(1 for result in results if result.passed)
    failed = len(results) - passed
    return f"{passed}/{len(results)} passed, {failed} failed"


def report_results(
    results: List[TestResult],
    suite_groups=None,
    selected_suite: Optional[str] = None,
    total_elapsed: Optional[float] = None,
) -> None:
    """Print a layered report by case, suite membership, and case directory."""
    passed_results = [result for result in results if result.passed]
    failed_results = [result for result in results if not result.passed]

    print()
    print("=" * 78)
    print("OneFLOW Regression Test Report")
    print("=" * 78)
    if selected_suite:
        print(f"Selected suite : {selected_suite}")
    print(f"Overall        : {_summarize_results(results)}")
    if total_elapsed is not None:
        print(f"Elapsed        : {total_elapsed:.2f} s")
    print(f"Suite groups   : {len(suite_groups or [])}")
    print()

    if suite_groups:
        print("By suite manifest (groups may overlap; totals are not additive):")
        for group_name, group_cases in suite_groups:
            normalized = {_normalized_case_path(path) for path in group_cases}
            group_results = [result for result in results if _normalized_case_path(result.case_dir) in normalized]
            print(f"  {group_name:<28} {_summarize_results(group_results)}")
        print()

    directory_groups = {}
    for result in results:
        parts = result.case_dir.replace("\\", "/").split("/")
        category = parts[0] if len(parts) > 1 else "(root)"
        directory_groups.setdefault(category, []).append(result)
    print("By case directory:")
    for category in sorted(directory_groups):
        print(f"  {category + '/':<28} {_summarize_results(directory_groups[category])}")
    print()

    print(f"Passed ({len(passed_results)}):")
    if passed_results:
        for result in passed_results:
            print(f"  [PASS] {result.case_dir}")
    else:
        print("  (none)")

    print(f"Failed ({len(failed_results)}):")
    if failed_results:
        for result in failed_results:
            print(f"  [FAIL] {result.case_dir}")
            print(f"         {result.failure_summary or 'unknown reason'}")
    else:
        print("  (none)")

    print("-" * 78)
    status = "PASS" if not failed_results else "FAIL"
    print(f"FINAL STATUS: {status} | {len(passed_results)} passed / {len(failed_results)} failed / {len(results)} total")
    print("=" * 78)


# ---------------------------------------------------------------------------
# --dry-run: self-contained fake cases used to exercise manifest validation,
# process execution, output comparison, and final reporting without a real
# OneFLOW binary or suite file. Failure scenarios target distinct assertions.
# ---------------------------------------------------------------------------

_DUMMY_RESULT_RELATIVE_PATH = "results/residual.dat"
_DUMMY_BINARY_RELATIVE_PATH = "grid/generated.ofl"

_DUMMY_SCENARIOS = [
    {
        "name": "dummy_case_binary_grid_ok",
        "baseline_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n1 1.1\n',
        "result_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n1 1.1\n',
        "binary_hex": "0001027f80feff004f464c",
        "exit_code": 0,
    },
    {
        "name": "dummy_case_ok",
        "baseline_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n1 1.1\n',
        "result_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n1 1.1\n',
        "exit_code": 0,
    },
    {
        "name": "dummy_case_content_mismatch",
        "baseline_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n1 1.1\n',
        "result_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n1 9.9\n',
        "exit_code": 0,
    },
    {
        "name": "dummy_case_solver_crash",
        "baseline_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n',
        "result_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n',
        "exit_code": 1,
    },
    {
        "name": "dummy_case_missing_output",
        "baseline_content": 'Variables=\n"iter"\n"rho"\n0 1.0\n',
        "result_content": None,  # the fake solver won't write anything
        "exit_code": 0,
    },
    {
        "name": "dummy_case_no_output_assertions",
        "baseline_content": 'Variables=\\n"iter"\\n"rho"\\n0 1.0\\n',
        "result_content": 'Variables=\\n"iter"\\n"rho"\\n0 1.0\\n',
        "exit_code": 0,
        "no_assertions": True,
    },
]


def _write_fake_solver(root: str) -> str:
    """Create a tiny stand-in "solver" script for --dry-run.

    It reads a scenario.json placed in the test-case directory (see
    build_dummy_suite) and, based on that, writes a result file and/or
    exits with a chosen code - simulating a real solver's pass/fail
    behavior without needing OneFLOW itself.
    """
    stub_path = os.path.join(root, "fake_solver.py")
    stub_source = (
        "import json, os, sys\n"
        "project_dir = sys.argv[2]\n"
        "with open(os.path.join(project_dir, 'scenario.json'), encoding='utf-8') as f:\n"
        "    scenario = json.load(f)\n"
        "if scenario.get('write_result'):\n"
        "    result_path = os.path.join(project_dir, scenario['result_relative_path'])\n"
        "    os.makedirs(os.path.dirname(result_path), exist_ok=True)\n"
        "    with open(result_path, 'w', encoding='utf-8') as f:\n"
        "        f.write(scenario['result_content'])\n"
        "if scenario.get('binary_hex') is not None:\n"
        "    binary_path = os.path.join(project_dir, scenario['binary_relative_path'])\n"
        "    os.makedirs(os.path.dirname(binary_path), exist_ok=True)\n"
        "    with open(binary_path, 'wb') as f:\n"
        "        f.write(bytes.fromhex(scenario['binary_hex']))\n"
        "sys.exit(scenario.get('exit_code', 0))\n"
    )
    with open(stub_path, "w", encoding="utf-8") as f:
        f.write(stub_source)
    return stub_path


def build_dummy_suite(root: str):
    """Create fake cases for --dry-run, including a binary-output comparison."""
    fake_solver_path = _write_fake_solver(root)
    case_dirs = []
    for scenario in _DUMMY_SCENARIOS:
        case_dir = os.path.join(root, scenario["name"])
        os.makedirs(os.path.join(case_dir, "autotest"), exist_ok=True)
        if not scenario.get("no_assertions"):
            with open(os.path.join(case_dir, "autotest", "compare.txt"), "w", encoding="utf-8") as f:
                f.write(_DUMMY_RESULT_RELATIVE_PATH + "\n")
        with open(os.path.join(case_dir, "autotest", "residual.dat"), "w", encoding="utf-8") as f:
            f.write(scenario["baseline_content"])
        scenario_data = {
            "result_relative_path": _DUMMY_RESULT_RELATIVE_PATH,
            "exit_code": scenario["exit_code"],
            "write_result": scenario["result_content"] is not None,
            "result_content": scenario["result_content"] or "",
        }
        if scenario.get("binary_hex") is not None:
            scenario_data["binary_relative_path"] = _DUMMY_BINARY_RELATIVE_PATH
            scenario_data["binary_hex"] = scenario["binary_hex"]
            with open(os.path.join(case_dir, "autotest", "binary.txt"), "w", encoding="utf-8") as f:
                f.write(_DUMMY_BINARY_RELATIVE_PATH + "\n")
            baseline_binary = os.path.join(case_dir, "autotest", os.path.basename(_DUMMY_BINARY_RELATIVE_PATH))
            os.makedirs(os.path.dirname(baseline_binary), exist_ok=True)
            with open(baseline_binary, "wb") as f:
                f.write(bytes.fromhex(scenario["binary_hex"]))
        with open(os.path.join(case_dir, "scenario.json"), "w", encoding="utf-8") as f:
            json.dump(scenario_data, f)
        case_dirs.append(case_dir)

    mpi_cmd = sys.executable
    exe_cmd = fake_solver_path
    return case_dirs, mpi_cmd, exe_cmd


def _env_float(name: str) -> Optional[float]:
    value = os.environ.get(name)
    return float(value) if value is not None else None


def build_arg_parser() -> argparse.ArgumentParser:
    # REFACTOR (step 4): replaces manual sys.argv[1]/[2]/... indexing (which
    # raised an unhelpful IndexError if an argument was missing) with
    # argparse, and adds --dry-run / --timeout.
    parser = argparse.ArgumentParser(description="OneFLOW autotest runner")
    parser.add_argument(
        "mpi_cmd", nargs="?", help="MPI launcher, e.g. 'mpirun -np 1' (not needed with --dry-run)"
    )
    parser.add_argument(
        "exe_cmd", nargs="?", help="Path to the OneFLOW executable (not needed with --dry-run)"
    )
    parser.add_argument(
        "suite_file", nargs="?", default="suites/cpu-serial.txt",
        help="Path to the suite list file (default: suites/cpu-serial.txt)",
    )
    parser.add_argument(
        "residual_db", nargs="?", default=os.environ.get("ONEFLOW_RESIDUAL_DB"),
        help="Path to the residual baseline JSON database",
    )
    parser.add_argument(
        "residual_tolerance", nargs="?", type=float,
        default=_env_float("ONEFLOW_RESIDUAL_TOLERANCE"),
        help="Absolute tolerance override for residual comparison",
    )
    parser.add_argument(
        "--timeout", type=float, default=_env_float("ONEFLOW_TEST_TIMEOUT_SECONDS"),
        help="Per-case timeout in seconds (default: no timeout, i.e. wait forever)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help=(
            "Run a handful of self-contained fake test cases (one guaranteed to "
            "pass, three each failing for a different reason) instead of a real "
            "suite. No OneFLOW binary or suite file is needed - useful for "
            "smoke-testing this script itself and demonstrating the report format."
        ),
    )
    return parser


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    location = os.getcwd()
    logger.debug("location=%s", location)
    list_directory_contents(location)

    dry_run_dir = None
    suite_groups = []
    selected_suite = args.suite_file if not args.dry_run else "--dry-run"
    if args.dry_run:
        dry_run_dir = tempfile.mkdtemp(prefix="oneflow_test_dryrun_")
        case_dirs, mpi_cmd, exe_cmd = build_dummy_suite(dry_run_dir)
        residual_db, residual_tolerance = None, None
        logger.info(
            "Running --dry-run with %d fake case(s) in %s "
            "(no real suite file or OneFLOW binary is used)",
            len(case_dirs), dry_run_dir,
        )
    else:
        if not args.mpi_cmd or not args.exe_cmd:
            parser.error("mpi_cmd and exe_cmd are required unless --dry-run is given")
        try:
            case_dirs = read_suite_file(args.suite_file)
            suite_groups = read_suite_groups(args.suite_file)
        except (OSError, ValueError) as error:
            parser.error(str(error))
        mpi_cmd, exe_cmd = args.mpi_cmd, args.exe_cmd
        residual_db, residual_tolerance = args.residual_db, args.residual_tolerance

    logger.info("suite size=%s", len(case_dirs))
    logger.info("residual_db=%s", residual_db)
    logger.info("residual_tolerance=%s", residual_tolerance)
    logger.info("backend=%s", os.environ.get("ONEFLOW_ACCEL_BACKEND", "CPU"))
    logger.info(
        "residual_profile=%s",
        "strict"
        if os.environ.get("ONEFLOW_RESIDUAL_TEST_OUTPUT", "").lower() in {"1", "true", "on"}
        else "normal",
    )
    logger.info("process_timeout=%s", args.timeout)

    suite_start = time.perf_counter()
    try:
        results = run_all_test(
            case_dirs, mpi_cmd, exe_cmd, residual_db, residual_tolerance, args.timeout
        )
    finally:
        if dry_run_dir:
            shutil.rmtree(dry_run_dir, ignore_errors=True)

    report_results(
        results,
        suite_groups=suite_groups,
        selected_suite=selected_suite,
        total_elapsed=time.perf_counter() - suite_start,
    )
    n_fail = sum(1 for r in results if not r.passed)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
