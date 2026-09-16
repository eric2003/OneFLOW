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
        if not is_number(token1):
            continue
        if abs(float(token1) - float(token2)) > FLOAT_ABS_TOLERANCE:
            return False
    return True


MAX_REPORTED_LINE_DIFFS = 20


def cmp_file(
    file_name1: str,
    file_name2: str,
    max_reported_diffs: int = MAX_REPORTED_LINE_DIFFS,
) -> bool:
    """Line-by-line comparison of two result files.

    BUGFIX (step 2): collects every differing line instead of stopping at
    the first one, and uses `with` so a mid-read error can't leak a file
    handle.
    BUGFIX (step 4): missing/unreadable files (e.g. the solver never wrote
    its output) are now reported as a normal failure instead of raising an
    uncaught FileNotFoundError - this matters for --dry-run's
    "missing output file" scenario, and can happen for real if a case
    crashes before writing anything.
    """
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
    """Read one relative result-file path per line from `filename`."""
    with open(filename, "r", encoding="utf-8-sig") as f:
        return [line.strip() for line in f.readlines()]


@dataclass
class ComparisonPair:
    """One (baseline, generated) result-file pair to compare for a test case."""

    baseline_path: str
    generated_path: str


def build_comparison_pairs(
    test_project_dir: str, relative_names: List[str]
) -> List[ComparisonPair]:
    """Build the list of (baseline, generated) file pairs for one test case.

    REFACTOR (step 3): replaces three index-aligned parallel lists from the
    original script with a single list of paired values, removing the risk
    of the lists silently drifting out of alignment.
    """
    pairs = []
    for relative_name in relative_names:
        base_name = relative_name.split("/")[-1]
        baseline_file = test_project_dir + "/autotest/" + base_name
        generated_file = test_project_dir + "/" + relative_name
        pairs.append(
            ComparisonPair(
                baseline_path=os.path.normpath(os.path.abspath(baseline_file)),
                generated_path=os.path.normpath(os.path.abspath(generated_file)),
            )
        )
    return pairs


@dataclass
class TestResult:
    """Outcome of running one test case, with enough detail to say not just
    "did it pass" but "what exactly failed", for reporting purposes.
    """

    case_name: str
    passed: bool
    return_code: Optional[int] = None
    files_ok: Optional[bool] = None
    residual_ok: Optional[bool] = None
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
        if not cmp_file(pair.baseline_path, pair.generated_path):
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
    residual_result = compare_case(
        test_project_dir, residual_db_path, absolute_tolerance=residual_tolerance
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
    logger.debug("test_project_dir=%s", test_project_dir)

    test_script = test_project_dir + "/autotest/test.txt"
    abs_test_script = os.path.normpath(os.path.abspath(test_script))
    relative_names = get_filename_list(abs_test_script)
    comparison_pairs = build_comparison_pairs(test_project_dir, relative_names)

    cmd = mpi_cmd + " " + exe_cmd + " 0 " + test_project_dir
    logger.debug("cmd=%s", cmd)

    return_code = _run_solver_process(cmd, process_timeout)
    files_ok = _compare_result_files(comparison_pairs)
    residual_ok = _compare_residuals(test_project_dir, residual_db_path, residual_tolerance)

    failure_reasons = []
    if return_code != 0:
        failure_reasons.append(f"solver exit code={return_code}")
    if not files_ok:
        failure_reasons.append("result file mismatch")
    if residual_ok is False:
        failure_reasons.append("residual baseline mismatch")

    result = TestResult(
        case_name=case_name,
        passed=not failure_reasons,
        return_code=return_code,
        files_ok=files_ok,
        residual_ok=residual_ok,
        failure_reasons=failure_reasons,
    )
    logger.info(
        "case=%s return_code=%s files_ok=%s residual_ok=%s passed=%s",
        case_name, return_code, files_ok, residual_ok, result.passed,
    )
    return result


def read_suite_file(filename: str) -> List[str]:
    """Read one test-case directory per line from a suite list file."""
    with open(filename, "r", encoding="utf-8-sig") as f:
        return [line.strip() for line in f.readlines() if line.strip()]


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


def report_results(results: List[TestResult]) -> None:
    """Print a summary of which cases passed and which failed, by name -
    not just a pass/fail count.
    """
    passed_names = [r.case_name for r in results if r.passed]
    failed_results = [r for r in results if not r.passed]

    logger.info(
        "Passed (%d): %s", len(passed_names), ", ".join(passed_names) or "(none)"
    )
    if failed_results:
        logger.info("Failed (%d):", len(failed_results))
        for r in failed_results:
            logger.info("  - %s -> %s", r.case_name, r.failure_summary or "unknown reason")
    else:
        logger.info("Failed (0): (none)")

    if not failed_results:
        logger.info("Total tests passed!")
    else:
        logger.info("ERROR: Some tests failed")
    logger.info("%s tests passed! %s tests failed!", len(passed_names), len(failed_results))


# ---------------------------------------------------------------------------
# --dry-run: a handful of self-contained fake test cases, used to exercise
# run_test/run_all_test/report_results without needing a real OneFLOW
# binary or a real suite file. Each scenario is designed to fail for a
# DIFFERENT reason, so the summary output demonstrates that the pipeline
# correctly tells different kinds of failures apart.
# ---------------------------------------------------------------------------

_DUMMY_RESULT_RELATIVE_PATH = "results/residual.dat"

_DUMMY_SCENARIOS = [
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
        "sys.exit(scenario.get('exit_code', 0))\n"
    )
    with open(stub_path, "w", encoding="utf-8") as f:
        f.write(stub_source)
    return stub_path


def build_dummy_suite(root: str):
    """Create a handful of fake test cases under `root` for --dry-run.

    Returns (case_dirs, mpi_cmd, exe_cmd) that can be passed straight to
    run_all_test() - no suite file or real OneFLOW binary needed.
    """
    fake_solver_path = _write_fake_solver(root)
    case_dirs = []
    for scenario in _DUMMY_SCENARIOS:
        case_dir = os.path.join(root, scenario["name"])
        os.makedirs(os.path.join(case_dir, "autotest"), exist_ok=True)
        with open(os.path.join(case_dir, "autotest", "test.txt"), "w", encoding="utf-8") as f:
            f.write(_DUMMY_RESULT_RELATIVE_PATH + "\n")
        with open(os.path.join(case_dir, "autotest", "residual.dat"), "w", encoding="utf-8") as f:
            f.write(scenario["baseline_content"])
        with open(os.path.join(case_dir, "scenario.json"), "w", encoding="utf-8") as f:
            json.dump(
                {
                    "result_relative_path": _DUMMY_RESULT_RELATIVE_PATH,
                    "exit_code": scenario["exit_code"],
                    "write_result": scenario["result_content"] is not None,
                    "result_content": scenario["result_content"] or "",
                },
                f,
            )
        case_dirs.append(case_dir)

    mpi_cmd = sys.executable  # reuse the current Python interpreter as the "launcher"
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
        case_dirs = read_suite_file(args.suite_file)
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

    try:
        results = run_all_test(
            case_dirs, mpi_cmd, exe_cmd, residual_db, residual_tolerance, args.timeout
        )
    finally:
        if dry_run_dir:
            shutil.rmtree(dry_run_dir, ignore_errors=True)

    report_results(results)
    n_fail = sum(1 for r in results if not r.passed)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
