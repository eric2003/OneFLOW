from __future__ import annotations

import unittest

from .cmake_test_catalog import CMakeTestTarget
from .ctest_discovery import CTestTest, CTestDiscoveryResult, CTestDiscoverySnapshot, CTestVersion
from .ctest_evidence import match_configured_ctest_evidence, format_ctest_evidence_report


def target(discovery=True, mode="PRE_TEST"):
    return CMakeTestTarget(
        target="StlMeshIOTest",
        cmake_file="tests/geometry/CMakeLists.txt",
        source_file="tests/geometry/StlMeshIOTest.cpp",
        discovery_declared=discovery,
        discovery_mode=mode,
    )


def ctest(name="StlMeshIOTest.NativeStlToGmshGeneratesSingleHexahedron",
          command=("build/bin/StlMeshIOTest",), backtrace=None):
    return CTestTest(
        name=name, configuration="Debug", command=command,
        properties=(), backtrace_file=backtrace, backtrace_line=32,
    )


class ConfiguredCTestEvidenceTest(unittest.TestCase):
    def test_matches_exact_identity_and_executable_target(self):
        result = match_configured_ctest_evidence(
            "mesh/stl", "StlMeshIOTest.NativeStlToGmshGeneratesSingleHexahedron",
            "tests/geometry/StlMeshIOTest.cpp", (target(),), (ctest(),),
        )
        self.assertEqual(result.status, "discovered")
        self.assertEqual(result.targets, ("StlMeshIOTest",))

    def test_matches_exact_identity_using_cmake_backtrace_when_command_missing(self):
        test = ctest(command=(), backtrace="/checkout/tests/geometry/CMakeLists.txt")
        result = match_configured_ctest_evidence(
            "mesh/stl", test.name, "tests/geometry/StlMeshIOTest.cpp",
            (target(),), (test,),
        )
        self.assertEqual(result.status, "discovered")

    def test_does_not_accept_similar_test_names(self):
        test = ctest(name="StlMeshIOTest.NativeStlToGmshGeneratesSingleHexahedron/0")
        result = match_configured_ctest_evidence(
            "mesh/stl", "StlMeshIOTest.NativeStlToGmshGeneratesSingleHexahedron",
            "tests/geometry/StlMeshIOTest.cpp", (target(mode="POST_BUILD"),), (test,),
        )
        self.assertEqual(result.status, "not_discovered")

    def test_reports_unmapped_target_missing_discovery_and_provenance_mismatch(self):
        self.assertEqual(
            match_configured_ctest_evidence("c", "Suite.Test", "x.cpp", (), ()).status,
            "target_unmapped",
        )
        self.assertEqual(
            match_configured_ctest_evidence("c", "Suite.Test", "x.cpp",
                                            (target(False),), ()).status,
            "discovery_not_declared",
        )
        mismatch = ctest(command=("OtherTest",))
        self.assertEqual(
            match_configured_ctest_evidence("c", mismatch.name, "x.cpp",
                                            (target(),), (mismatch,)).status,
            "target_provenance_mismatch",
        )

    def test_pre_test_discovery_absence_is_reported_as_deferred(self):
        result = match_configured_ctest_evidence(
            "c", "Suite.Test", "x.cpp", (target(),), ()
        )
        self.assertEqual(result.status, "discovery_deferred")
        self.assertIn("run CTest once", result.diagnostic)

    def test_reports_absent_configured_test_without_calling_it_a_failure(self):
        result = match_configured_ctest_evidence(
            "c", "Suite.Test", "x.cpp", (target(mode="POST_BUILD"),), ()
        )
        self.assertEqual(result.status, "not_discovered")
        self.assertIn("configured CTest list", result.diagnostic)

    def test_report_states_that_discovery_does_not_execute_tests(self):
        discovery = CTestDiscoveryResult(
            "discovered", "build", "Debug",
            snapshot=CTestDiscoverySnapshot(CTestVersion(3, 30, 0), ()),
        )
        report = format_ctest_evidence_report(discovery, ())
        self.assertIn("No tests were executed", report)


if __name__ == "__main__":
    unittest.main()
