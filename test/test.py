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
import os
import sys
import logging
import subprocess
import time

from residual_db import compare_case

# NOTE: default level is INFO so routine debug noise (variable dumps, argv
# echoes, etc.) is hidden by default. Set ONEFLOW_TEST_LOG_LEVEL=DEBUG to see
# everything that used to be printed unconditionally in the old script.
logging.basicConfig(
    level=os.environ.get("ONEFLOW_TEST_LOG_LEVEL", "INFO"),
    format="%(levelname)s: %(message)s",
)
logger = logging.getLogger("oneflow.test")


def list_directory_contents(file_dir):
    # Was `my_dir_cmd`: pure debug helper, kept but demoted to DEBUG level
    # instead of always printing to stdout.
    for file in os.listdir(file_dir):
        logger.debug(file)


def is_number(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def cmp_float_str(str1, str2):
    # Was `CmpFloatStr`. Compares two whitespace-separated lines token by
    # token, treating numeric tokens as floats with a fixed tolerance and
    # everything else as exact string match (implicitly, via caller's
    # str1 != str2 check before this function is even invoked).
    tokens1 = str1.split()
    tokens2 = str2.split()
    for i in range(0, len(tokens1)):
        token1 = tokens1[i]
        token2 = tokens2[i]
        if not is_number(token1):
            continue
        value1 = float(token1)
        value2 = float(token2)
        if abs(value1 - value2) > 1.0e-8:
            return False
    return True


def cmp_file(file_name1, file_name2):
    # Was `CmpFile`. Line-by-line comparison of two result files.
    f1 = open(file_name1, "r", encoding="utf-8-sig")
    f2 = open(file_name2, "r", encoding="utf-8-sig")
    cmp_ok = True
    line_id = 0
    while True:
        str1 = f1.readline()
        str2 = f2.readline()
        if not str1 and not str2:
            break
        if not str1 or not str2:
            logger.info(file_name1)
            logger.info(file_name2)
            logger.info("Line %s has different end-of-file position", line_id)
            cmp_ok = False
            break
        if str1 != str2:
            if not cmp_float_str(str1, str2):
                logger.info(file_name1)
                logger.info(file_name2)
                logger.info("Line %s", line_id)
                logger.info("line1=%s", str1)
                logger.info("line2=%s", str2)
                cmp_ok = False
                break
        line_id += 1
    f1.close()
    f2.close()
    return cmp_ok


def get_filename_list(filename, filename_list):
    # Was `GetFileNameList`. NOTE: still mutates the list passed in
    # (out-parameter style) - this will be replaced with a return value
    # in a later refactor step.
    with open(filename, "r", encoding="utf-8-sig") as f:
        for line in f.readlines():
            filename_list.append(line.strip())


def run_test(test_project_dir, residual_db_path=None, residual_tolerance=None):
    # Was `RunTest`.
    logger.debug("test_project_dir=%s", test_project_dir)
    test_script = test_project_dir + "/autotest/test.txt"
    abs_test_script = os.path.normpath(os.path.abspath(test_script))
    logger.debug("test_script=%s", test_script)
    logger.debug("abs_test_script=%s", abs_test_script)

    relative_names = []
    get_filename_list(abs_test_script, relative_names)

    base_names = [name.split("/")[-1] for name in relative_names]

    baseline_files = []
    baseline_file_paths = []
    for base_name in base_names:
        baseline_file = test_project_dir + "/autotest/" + base_name
        baseline_files.append(baseline_file)
        baseline_file_paths.append(os.path.normpath(os.path.abspath(baseline_file)))

    generated_file_paths = []
    for relative_name in relative_names:
        generated_file = test_project_dir + "/" + relative_name
        generated_file_paths.append(os.path.normpath(os.path.abspath(generated_file)))

    logger.debug("argv (%s): %s", len(sys.argv), sys.argv)
    mpi_cmd = sys.argv[1]
    exe_cmd = sys.argv[2]
    cmd = mpi_cmd + " " + exe_cmd + " 0 " + test_project_dir
    logger.debug("cmd=%s", cmd)

    process = subprocess.Popen(cmd, shell=True)
    while process.poll() is None:
        time.sleep(0.5)
    return_code = process.poll()
    total_pass = return_code == 0

    for baseline_path, generated_path in zip(baseline_file_paths, generated_file_paths):
        cmp_ok = cmp_file(baseline_path, generated_path)
        if not cmp_ok:
            total_pass = False
            break

    if residual_db_path:
        residual_result = compare_case(
            test_project_dir, residual_db_path, absolute_tolerance=residual_tolerance
        )
        logger.info("residual_baseline=%s", residual_result)
        if not residual_result["ok"]:
            total_pass = False

    logger.info("return_code=%s", return_code)
    logger.info("total_pass=%s", total_pass)
    return total_pass


def run_all_test(filename, residual_db_path=None, residual_tolerance=None):
    # Was `RunAllTest`.
    pass_flags = []
    with open(filename, "r", encoding="utf-8-sig") as f:
        for line in f.readlines():
            project_dir = line.strip()
            logger.info("project_dir=%s", project_dir)
            pass_flags.append(run_test(project_dir, residual_db_path, residual_tolerance))
    return pass_flags


def main():
    location = os.getcwd()
    logger.debug("location=%s", location)
    list_directory_contents(location)

    suite_file = "suites/cpu-serial.txt"
    if len(sys.argv) >= 4:
        suite_file = sys.argv[3]

    residual_db = os.environ.get("ONEFLOW_RESIDUAL_DB")
    if len(sys.argv) >= 5:
        residual_db = sys.argv[4]

    residual_tolerance = os.environ.get("ONEFLOW_RESIDUAL_TOLERANCE")
    if len(sys.argv) >= 6:
        residual_tolerance = sys.argv[5]
    if residual_tolerance is not None:
        residual_tolerance = float(residual_tolerance)

    logger.info("suite_file=%s", suite_file)
    logger.info("residual_db=%s", residual_db)
    logger.info("residual_tolerance=%s", residual_tolerance)
    logger.info("backend=%s", os.environ.get("ONEFLOW_ACCEL_BACKEND", "CPU"))
    logger.info(
        "residual_profile=%s",
        "strict"
        if os.environ.get("ONEFLOW_RESIDUAL_TEST_OUTPUT", "").lower() in {"1", "true", "on"}
        else "normal",
    )

    pass_flags = run_all_test(suite_file, residual_db, residual_tolerance)
    n_test = len(pass_flags)
    n_pass = sum(1 for flag in pass_flags if flag)
    n_fail = n_test - n_pass

    if n_pass == n_test:
        logger.info("Total tests passed!")
    else:
        logger.info("ERROR: Some tests failed")
    logger.info("%s tests passed! %s tests failed!", n_pass, n_fail)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
