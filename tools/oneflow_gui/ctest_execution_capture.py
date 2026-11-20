"""Capture per-test execution evidence from CTest's JUnit report."""
from __future__ import annotations
import hashlib, math, platform, uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from ..verification.execution_records import ExecutionArtifact, ExecutionRecord, ExecutionRecordStore
from ..verification.ctest_discovery import discover_configured_ctest
from ..verification.ctest_evidence import audit_configured_ctest_evidence
from .execution_capture import _git_value, _repository_name
from .test_runner import TestCommand

def is_ctest_command(command: TestCommand) -> bool:
    if not command.arguments: return False
    exe=command.arguments[0].replace("\\","/").rsplit("/",1)[-1].casefold()
    return exe in {"ctest","ctest.exe"} and "--output-junit" in command.arguments

def ctest_junit_path(command: TestCommand) -> Path:
    try: return Path(command.arguments[command.arguments.index("--output-junit")+1]).expanduser()
    except (ValueError,IndexError) as exc: raise ValueError("CTest command has no --output-junit path") from exc

def _parse_cases(path: Path):
    try: root=ET.parse(path).getroot()
    except (OSError,ET.ParseError) as exc: raise ValueError(f"unable to read CTest JUnit report: {exc}") from exc
    cases=[]; seen=set()
    for node in root.iter("testcase"):
        name=node.attrib.get("name","").strip()
        if not name: raise ValueError("CTest JUnit testcase has no name")
        if name in seen: raise ValueError(f"duplicate CTest JUnit testcase: {name}")
        seen.add(name); status=node.attrib.get("status","").casefold()
        failure=node.find("failure"); error=node.find("error"); skipped=node.find("skipped")
        details=[]
        if failure is not None or error is not None or status=="fail":
            outcome="failed"; detail=failure if failure is not None else error
            if detail is not None:
                if detail.attrib.get("message"): details.append(detail.attrib["message"])
                if (detail.text or "").strip(): details.append(detail.text.strip())
            output=(node.findtext("system-out") or "").strip()
            if output: details.append(output[-4000:])
        elif skipped is not None or status in {"notrun","disabled"}:
            outcome="blocked"
            if skipped is not None:
                if skipped.attrib.get("message"): details.append(skipped.attrib["message"])
                if (skipped.text or "").strip(): details.append(skipped.text.strip())
            else: details.append("CTest marked this test as not run.")
        else: outcome="passed"
        try: duration=float(node.attrib.get("time","0"))
        except (TypeError,ValueError) as exc: raise ValueError(f"invalid duration for CTest testcase {name}") from exc
        if not math.isfinite(duration) or duration<0: raise ValueError(f"invalid duration for CTest testcase {name}")
        cases.append((name,outcome,tuple(details),duration))
    return tuple(cases)

def _artifact(path: Path, role: str, root: Path) -> ExecutionArtifact:
    digest=hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""): digest.update(chunk)
    resolved=path.resolve()
    try: value=resolved.relative_to(root).as_posix()
    except ValueError: value=str(resolved)
    return ExecutionArtifact(role,value,digest.hexdigest())

def _source(test, root: Path, build: Path):
    source=test.backtrace_file if test else None
    line=test.backtrace_line if test else None
    if source:
        path=Path(source)
        if not path.is_absolute():
            root_path=root/path; build_path=build/path
            path=root_path if root_path.exists() else build_path
        try: return path.resolve().relative_to(root).as_posix(),line
        except ValueError: return str(path),line
    return (build/"CTestTestfile.cmake").resolve().relative_to(root).as_posix(),None

def capture_ctest_execution(command: TestCommand, *, repository_root: Path, build_directory: Path,
    junit_path: Path, records_path: Path, artifact_directory: Path, started_at: str,
    finished_at: str, exit_code: int, output_lines, configuration: str=""):
    """Persist one schema-v1 record per JUnit testcase; link only audited exact identities."""
    root=Path(repository_root).expanduser().resolve(); build=Path(build_directory).expanduser().resolve()
    report=Path(junit_path).expanduser().resolve()
    if not report.is_file(): raise FileNotFoundError(f"CTest JUnit report was not produced: {report}")
    cases=_parse_cases(report)
    if not cases: return ()
    discovery=discover_configured_ctest(build,configuration or None,ctest_executable=command.arguments[0])
    discovered=(); links={}; source_paths={}; unresolved={}; runtime=""
    if discovery.status=="discovered" and discovery.snapshot is not None:
        discovered=discovery.snapshot.tests; runtime="CTest "+discovery.snapshot.version.text
        try: matches=audit_configured_ctest_evidence(root,discovered); audit_error=""
        except (OSError,RuntimeError,ValueError) as exc: matches=(); audit_error=str(exc)
        for item in matches:
            if item.status=="discovered" and item.discovered_name:
                links.setdefault(item.discovered_name,[]).append(item.contract_id)
                if item.source_path:
                    source_paths.setdefault(item.discovered_name,set()).add(item.source_path)
            elif item.reference:
                unresolved.setdefault(item.reference,[]).append(
                    f"Contract {item.contract_id} not linked ({item.status}): {item.diagnostic}")
        if audit_error: unresolved.setdefault("",[]).append("Contract evidence audit unavailable: "+audit_error)
    else:
        unresolved.setdefault("",[]).append("Configured CTest discovery unavailable: "+(discovery.diagnostic or discovery.status))
    test_map={item.name:item for item in discovered}
    invocation=str(uuid.uuid4()); artifact_root=Path(artifact_directory).expanduser(); artifact_root.mkdir(parents=True,exist_ok=True)
    log=artifact_root/(invocation+"-ctest.log")
    data=("\n".join(output_lines)+("\n" if output_lines else "")).encode("utf-8"); log.write_bytes(data)
    artifacts=(_artifact(report,"ctest-junit",root),_artifact(log,"stdout",root))
    repository=_repository_name(_git_value(root,"remote","get-url","origin"))
    commit=_git_value(root,"rev-parse","HEAD"); branch=_git_value(root,"branch","--show-current")
    records=[]
    for name,outcome,details,_duration in cases:
        candidates=source_paths.get(name,set())
        diagnostics=list(details)+unresolved.get(name,[])+unresolved.get("",[])
        if len(candidates)==1:
            source_path=next(iter(candidates)); source_line=None
        else:
            source_path,source_line=_source(test_map.get(name),root,build)
            if len(candidates)>1:
                diagnostics.append("Multiple Contract evidence source paths match this CTest identity; using CMake registration provenance.")
        records.append(ExecutionRecord(
            test_id="ctest:"+name,framework="CTest",test_name=name,source_path=source_path,
            source_line=source_line,command=tuple(command.arguments),outcome=outcome,
            started_at=started_at,finished_at=finished_at,
            exit_code=0 if outcome=="passed" else 1 if outcome=="failed" else None,
            repository=repository,commit=commit,branch=branch,configuration=configuration,
            platform=platform.platform(),runtime_version=runtime,working_directory=str(build),
            contract_ids=tuple(sorted(set(links.get(name,())))),artifacts=artifacts,
            diagnostics=tuple(item for item in diagnostics if item),invocation_id=invocation))
    store=ExecutionRecordStore(Path(records_path))
    for record in records: store.append(record)
    return tuple(records)
