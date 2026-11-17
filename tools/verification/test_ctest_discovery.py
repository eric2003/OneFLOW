from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from .ctest_discovery import (
    CTestDiscoveryError, discover_configured_ctest, parse_ctest_discovery_json,
)


def payload():
    return {
        "kind": "ctestInfo",
        "version": {"major": 3, "minor": 29, "patch": 0},
        "backtraceGraph": {
            "files": ["tests/runtime/CMakeLists.txt"],
            "nodes": [
                {"file": 0, "line": 12, "command": 0, "parent": 1},
                {"file": 0, "line": 20, "command": 0},
            ],
            "commands": ["add_test"],
        },
        "tests": [{
            "name": "lifecycle_test/FullSimulationEntry", "config": "Debug",
            "command": ["OneFLOWTests.exe", "--gtest_filter=FullSimulationEntry"],
            "backtrace": 0,
            "properties": [{"name": "WORKING_DIRECTORY", "value": "build"}],
        }],
    }


class CTestDiscoveryTest(unittest.TestCase):
    def test_parses_name_configuration_command_properties_and_backtrace(self):
        snapshot = parse_ctest_discovery_json(json.dumps(payload()))
        self.assertEqual(snapshot.version.text, "3.29.0")
        self.assertEqual(len(snapshot.tests), 1)
        test = snapshot.tests[0]
        self.assertEqual(test.name, "lifecycle_test/FullSimulationEntry")
        self.assertEqual(test.configuration, "Debug")
        self.assertEqual(test.command[0], "OneFLOWTests.exe")
        self.assertEqual(test.properties[0]["name"], "WORKING_DIRECTORY")
        self.assertEqual(test.backtrace_file, "tests/runtime/CMakeLists.txt")
        self.assertEqual(test.backtrace_line, 12)

    def test_uses_parent_line_when_child_node_has_no_line(self):
        data = payload()
        data["backtraceGraph"]["nodes"][0].pop("line")
        test = parse_ctest_discovery_json(json.dumps(data)).tests[0]
        self.assertEqual(test.backtrace_line, 20)

    def test_optional_command_and_backtrace_can_be_absent(self):
        data = payload()
        data.pop("backtraceGraph")
        data["tests"][0].pop("backtrace")
        data["tests"][0].pop("command")
        test = parse_ctest_discovery_json(json.dumps(data)).tests[0]
        self.assertIsNone(test.command)
        self.assertIsNone(test.backtrace_file)

    def test_empty_test_list_is_valid(self):
        data = payload()
        data["tests"] = []
        self.assertEqual(parse_ctest_discovery_json(json.dumps(data)).tests, ())

    def test_rejects_invalid_json_kind_version_and_test_shapes(self):
        with self.assertRaises(CTestDiscoveryError):
            parse_ctest_discovery_json("{")
        cases = []
        data = payload(); data["kind"] = "other"; cases.append(data)
        data = payload(); data["version"]["major"] = -1; cases.append(data)
        data = payload(); data["tests"] = None; cases.append(data)
        data = payload(); data["tests"][0]["name"] = ""; cases.append(data)
        data = payload(); data["tests"][0]["command"] = "invalid"; cases.append(data)
        data = payload(); data["tests"][0]["properties"] = ["invalid"]; cases.append(data)
        for data in cases:
            with self.subTest(data=data), self.assertRaises(CTestDiscoveryError):
                parse_ctest_discovery_json(json.dumps(data))

    def test_rejects_out_of_range_backtrace_and_parent_cycle(self):
        data = payload(); data["tests"][0]["backtrace"] = 99
        with self.assertRaises(CTestDiscoveryError):
            parse_ctest_discovery_json(json.dumps(data))
        data = payload()
        data["backtraceGraph"]["nodes"][0].pop("line")
        data["backtraceGraph"]["nodes"][0]["parent"] = 0
        with self.assertRaises(CTestDiscoveryError):
            parse_ctest_discovery_json(json.dumps(data))

    def test_runner_uses_list_only_mode_and_never_executes_tests(self):
        with tempfile.TemporaryDirectory() as directory:
            completed = type("Completed", (), {
                "returncode": 0, "stdout": json.dumps(payload()), "stderr": "",
            })()
            with patch("tools.verification.ctest_discovery.shutil.which",
                       return_value="/usr/bin/ctest"), patch(
                "tools.verification.ctest_discovery.subprocess.run",
                return_value=completed,
            ) as run:
                result = discover_configured_ctest(Path(directory), "Debug")
            self.assertEqual(result.status, "discovered")
            args = run.call_args.args[0]
            self.assertIn("-N", args)
            self.assertIn("--show-only=json-v1", args)
            self.assertNotIn("--output-on-failure", args)

    def test_reports_missing_build_directory_and_ctest(self):
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(
                discover_configured_ctest(Path(directory) / "missing").status,
                "build_directory_missing",
            )
            with patch("tools.verification.ctest_discovery.shutil.which", return_value=None):
                self.assertEqual(
                    discover_configured_ctest(Path(directory)).status, "ctest_unavailable"
                )


if __name__ == "__main__":
    unittest.main()
