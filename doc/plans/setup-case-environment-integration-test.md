# SetupCaseEnvironment integration-test investigation

Date: 2026-10-26
Branch: `hexin`
Inspected HEAD: [`82124625359b72f72910039807ffb209dcff3281`](https://github.com/eric2003/OneFLOW/commit/82124625359b72f72910039807ffb209dcff3281)

## Question

Can a test exercise the real `SimuContext::SetupCaseEnvironment()` configuration path with a temporary case and a real control file, without reimplementing solver startup?

## Verified call path

`SimuImp::Run() -> SetupProcessEnvironment() -> RunCase() -> PreProcess() -> InitSimu() -> SimuContext::SetupCaseEnvironment()`

The case setup implementation in `codes/main/src/SimuContextEnv.cpp` performs this sequence:

1. Marks case setup active so the existing exception path attempts cleanup after partial initialization.
2. Binds the explicit case directory through `Prj::SetPrjBaseDir(caseDir_)` and binds `logFile`.
3. Calls `ReadControlInfo(caseDir_)`.
4. Reads optional root settings from `RequireGlobalDataBase().RequireDataPara()`.
5. Applies project, generic, and mesh/flow-specific roots; conditionally reads `sourceGridCaseDir`.
6. Writes one or two `io_paths.txt` diagnostics, depending on output-root resolution.

## What makes the real entry point non-trivial

The case-specific overload `ReadControlInfo(caseDir)` in `codes/io/src/ParaFile.cpp` is not a pure file reader. It:

- Reads the script list and parameter files only on `Parallel::IsServer()`.
- Calls `Parallel::TestSayHelloFromEveryProcess()`.
- Broadcasts the control database with `HXBcast(CompressData, DecompressData, ...)`.
- Dumps the resulting parameter database to `log/database.log` under the case directory.

The broadcast callbacks resolve the process-global database through `RequireGlobalDataBase()`. Therefore a test that invokes the production overload needs a valid database singleton and a valid parallel runtime, or a deliberately narrow seam around the coordination boundary. Merely constructing `SimuContext` does not initialize either. The current `tests/main/simu_context_test.cpp` deliberately tests context/task plumbing without MPI or control-file I/O, so it is not an integration fixture today.

The control reader also has two overloads (implicit legacy project paths and explicit `caseDir` paths). The setup path correctly uses the explicit overload, which should remain the path under test rather than being replaced with a test-only parser.

## Smallest next experiment

Before changing `SimuContext` responsibilities, inspect the existing global database fixture/initialization API and the test target's linked libraries. If the database can be initialized independently and the test process is guaranteed serial, create a temporary case containing `script/control.txt` plus one parameter script, initialize the same database object used by `RequireGlobalDataBase()`, and call `SetupCaseEnvironment()` directly. Assert at least:

- a real control parameter is available from the global parameter database after setup;
- explicit `meshInputRoot` overrides `sourceGridCaseDir`;
- absent/empty mesh input selects the legacy source-grid fallback;
- a second case with omitted roots does not inherit the first case's roots;
- setup teardown clears case-bound state after success and after a deliberately missing control file.

Do not call `SetupProcessEnvironment()` unless the test is explicitly designed to exercise MPI/accelerator initialization. The intended first integration target is case configuration, not the whole solver or device runtime.

If database initialization cannot be performed without the process/parallel bootstrap, the next option is a narrow injectable control-loading boundary whose production implementation still calls `ReadControlInfo(caseDir)`. Do not copy the script parsing, broadcast, or database commit logic into the test. A change to `SimuContext` ownership is justified only if the dependency prevents this focused test and the same boundary is needed by production callers.

## Current outcome

**R1 remains Partial.** This inspection identifies the concrete integration dependencies, but does not yet establish a runnable integration fixture or change production behavior. No tests were run in this session. The reported 700/700 local pass predates the latest fallback-value test commit and is not evidence for the current HEAD.
