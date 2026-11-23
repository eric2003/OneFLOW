"""Tests for per-test CTest execution evidence capture."""
from __future__ import annotations
import hashlib,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from tools.oneflow_gui.ctest_execution_capture import capture_ctest_execution,ctest_junit_path,is_ctest_command
from tools.oneflow_gui.test_runner import TestCommand
from tools.verification.ctest_discovery import CTestDiscoveryResult,CTestDiscoverySnapshot,CTestTest,CTestVersion
from tools.verification.ctest_evidence import CTestEvidenceMatch

class CTestExecutionCaptureTest(unittest.TestCase):
    def _command(self,root):
        build=root/"build"
        return TestCommand("C++ unit tests (CTest)",("ctest","--output-on-failure","--output-junit",str(build/"ctest-results.xml"),"-C","Debug"),build)
    def test_detects_ctest_command_and_junit_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            command=self._command(Path(tmp))
            self.assertTrue(is_ctest_command(command))
            self.assertEqual(ctest_junit_path(command),Path(tmp)/"build"/"ctest-results.xml")
            self.assertFalse(is_ctest_command(TestCommand("Python",("python","-m","unittest"),Path(tmp))))
    def test_persists_cases_and_links_only_exact_discovered_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp).resolve(); build=root/"build"; build.mkdir()
            source=root/"tests"/"sample"/"SampleTest.cpp"; source.parent.mkdir(parents=True); source.write_text("// test\n",encoding="utf-8")
            cmake_source=source.parent/"CMakeLists.txt"; cmake_source.write_text("add_oneflow_gtest(SampleTests)\n",encoding="utf-8")
            command=self._command(root); report=ctest_junit_path(command)
            report.write_text('''<testsuite tests="3" failures="1" skipped="1">
<testcase name="Sample.Passes" status="run" time="0.01"><system-out>ok</system-out></testcase>
<testcase name="Sample.Fails" status="fail" time="0.02"><failure message="assertion failed">expected 1</failure><system-out>failure output</system-out></testcase>
<testcase name="Sample.Skips" status="notrun" time="0"><skipped message="disabled by configuration"/></testcase>
</testsuite>''',encoding="utf-8")
            tests=tuple(CTestTest(n,"Debug",("SampleTests.exe",),(),str(cmake_source),18) for n in ("Sample.Passes","Sample.Fails","Sample.Skips"))
            discovery=CTestDiscoveryResult("discovered",str(build),"Debug",CTestDiscoverySnapshot(CTestVersion(4,4,4),tests))
            matches=(CTestEvidenceMatch("runtime/process_lifecycle","Sample.Passes","tests/sample/SampleTest.cpp","discovered",("SampleTests",),"Sample.Passes","exact"),
                CTestEvidenceMatch("runtime/process_lifecycle","Sample.Fails","tests/sample/SampleTest.cpp","target_provenance_mismatch",("SampleTests",),diagnostic="wrong target"))
            records_path=root/".oneflow"/"records.jsonl"; logs=root/".oneflow"/"logs"
            def git(_root,*args):
                return {("remote","get-url","origin"):"https://github.com/eric2003/OneFLOW.git",("rev-parse","HEAD"):"abc123",("branch","--show-current"):"hexin"}.get(tuple(args),"")
            with patch("tools.oneflow_gui.ctest_execution_capture.discover_configured_ctest",return_value=discovery),patch("tools.oneflow_gui.ctest_execution_capture.audit_configured_ctest_evidence",return_value=matches),patch("tools.oneflow_gui.ctest_execution_capture._git_value",side_effect=git):
                records=capture_ctest_execution(command,repository_root=root,build_directory=build,junit_path=report,records_path=records_path,artifact_directory=logs,started_at="2026-11-20T10:00:00+00:00",finished_at="2026-11-20T10:00:03+00:00",exit_code=1,output_lines=("CTest failed",),configuration="Debug")
            self.assertEqual([r.outcome for r in records],["passed","failed","blocked"])
            self.assertEqual(records[0].contract_ids,("runtime/process_lifecycle",)); self.assertEqual(records[1].contract_ids,())
            self.assertEqual(records[0].source_path,"tests/sample/SampleTest.cpp"); self.assertIsNone(records[0].source_line)
            self.assertEqual(records[1].source_path,"tests/sample/CMakeLists.txt"); self.assertEqual(records[1].source_line,18)
            self.assertEqual(records[0].runtime_version,"CTest 4.4.4"); self.assertEqual(records[0].commit,"abc123")
            self.assertEqual(records[0].invocation_id,records[1].invocation_id)
            self.assertEqual(records[1].diagnostics[0],"assertion failed"); self.assertIn("disabled by configuration",records[2].diagnostics[0])
            self.assertEqual(records[0].artifacts[0].sha256,hashlib.sha256(report.read_bytes()).hexdigest())
            self.assertEqual(len(records_path.read_text(encoding="utf-8").splitlines()),3)
    def test_rejects_malformed_junit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); build=root/"build"; build.mkdir(); command=self._command(root); report=ctest_junit_path(command)
            report.write_text("<testsuite><testcase",encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"unable to read CTest JUnit report"):
                capture_ctest_execution(command,repository_root=root,build_directory=build,junit_path=report,records_path=root/"records.jsonl",artifact_directory=root/"logs",started_at="2026-11-20T10:00:00+00:00",finished_at="2026-11-20T10:00:01+00:00",exit_code=0,output_lines=(),configuration="Debug")

if __name__=="__main__": unittest.main()
