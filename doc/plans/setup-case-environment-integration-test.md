# SetupCaseEnvironment integration-test investigation

Date: 2026-10-26
Branch: `hexin`
Inspected source: [`SimuContextEnv.cpp`](https://github.com/eric2003/OneFLOW/blob/hexin/codes/main/src/SimuContextEnv.cpp), [`ParaFile.cpp`](https://github.com/eric2003/OneFLOW/blob/hexin/codes/io/src/ParaFile.cpp), [`DataBase.cpp`](https://github.com/eric2003/OneFLOW/blob/hexin/codes/database/src/DataBase.cpp)

## Question

Can a test exercise the real `SimuContext::SetupCaseEnvironment()` configuration path with a temporary case and a real control file, without reimplementing solver startup?

## Verified call path

`SimuImp::Run() -> SetupProcessEnvironment() -> RunCase() -> PreProcess() -> InitSimu() -> SimuContext::SetupCaseEnvironment()`

The case setup implementation in `codes/main/src/SimuContextEnv.cpp`:

1. Marks case setup active so the existing exception path attempts cleanup after partial initialization.
2. Binds the explicit case directory through `Prj::SetPrjBaseDir(caseDir_)` and binds `logFile`.
3. Calls `ReadControlInfo(caseDir_)`.
4. Reads optional root settings from `RequireGlobalDataBase().RequireDataPara()`.
5. Applies project, generic, and mesh/flow-specific roots; conditionally reads `sourceGridCaseDir`.
6. Writes one or two `io_paths.txt` diagnostics, depending on output-root resolution.

## Database initialization finding

`codes/database/src/DataBase.cpp` defines the global `std::unique_ptr<DataBase> globalDataBase` and a static `HXInitGlobalDataBase` object that constructs a `DataBase` during static initialization. Therefore, the database is not obviously dependent on `SetupProcessEnvironment()` merely to exist. This is a promising path toward a focused test, but it does not prove that all static initialization order, parallel calls, or linked test-target dependencies are safe.

The database exposes `RequireGlobalDataBase()`; its `DataPara` store can be cleared between cases. `SimuContext::TeardownCase()` also clears parameter and field stores, resets solver-related global state, clears the log case binding, and calls `Prj::ClearPrjBaseDir()`. The test fixture must respect that this cleanup is broader than a path-only reset and must not run concurrently with other tests that use these process-global objects.

## Why the real control reader needs care

The case-specific overload `ReadControlInfo(caseDir)` in `codes/io/src/ParaFile.cpp` is not a pure file reader. It:

- Reads the script list and parameter files only on `Parallel::IsServer()`.
- Calls `Parallel::TestSayHelloFromEveryProcess()`.
- Broadcasts the control database with `HXBcast(CompressData, DecompressData, ...)`.
- Dumps the resulting parameter database to `log/database.log` under the case directory.

The broadcast callbacks resolve the process-global database through `RequireGlobalDataBase()`. Thus, the test must determine whether the existing serial test runtime leaves `Parallel` in a valid single-process state. It must not call the production reader with an uninitialized MPI state just because the database singleton exists.

The current `tests/main/simu_context_test.cpp` tests context/task plumbing without control-file I/O. The production setup uses the explicit `ReadControlInfo(caseDir)` overload; keep that as the path under test rather than replacing it with a test-only parser.

## Concrete control-file fixture from rae2822_roe_sa

The supplied `rae2822_roe_sa` database dump is evidence that ten script files are merged into one parameter database. Its `script/control.txt` lists `simu.txt`, `dimension.txt`, `tolerence.txt`, `cfd.txt`, `grid.txt`, `vismodel.txt`, `multigrid.txt`, `inflow.txt`, `turb.txt`, and `walldist.txt`.

Use this as a source of representative contracts, not as a 119-entry golden dump. A small synthetic fixture is preferable for the first integration test because it isolates behavior and avoids requiring actual RAE2822 mesh and gas-model assets. Keep the fixture's real multi-file loading shape:

- `script/control.txt` includes at least a simulation task file, dimension file, and path/control file.
- One file defines `string simutask = "Solve";`.
- One file defines `int dimension = 2;` and `Real global_dt = 0.001;`.
- The path file defines `string sourceGridCaseDir = "../producer";` or an explicit `meshInputRoot`.

After setup, query the real database and assert representative values with their actual types: `simutask == "Solve"`, `dimension == 2`, `global_dt == 0.001`, and the selected mesh input root. Also verify the generated `log/database.log` exists. Do not run the Solve task; the test is for configuration only.

## Smallest next experiment

1. Inspect the test target's CMake sources/link dependencies and how `Parallel::IsServer()`, `TestSayHelloFromEveryProcess()`, and `HXBcast()` behave before process initialization.
2. If a serial invocation is valid, add a non-parallel temporary-case integration test that calls `SetupCaseEnvironment()` directly and reads the values from `RequireGlobalDataBase()`. Mark it non-parallel / serialize it if the test framework supports that.
3. Run two contexts sequentially. The first case defines explicit roots and a legacy source-grid fallback; the second omits them. Assert the second case does not inherit the first case's roots.
4. Ensure teardown after both success and an intentionally missing control file, while preserving the original setup exception.
5. Only if the process-global parallel boundary makes this impossible, introduce a narrow injectable control-loading boundary whose production implementation still calls `ReadControlInfo(caseDir)`. Do not copy parsing, broadcast, or database commit logic into the test.

Do not call `SetupProcessEnvironment()` unless the test is explicitly designed to exercise MPI/accelerator initialization. The intended first integration target is case configuration, not the whole solver or device runtime.

## Current outcome

**R1 remains Partial.** The global database is statically constructed, making a serial integration fixture plausible, but parallel-runtime validity and test-target linkage still need verification. No tests were run in this session. The reported 700/700 local pass predates the latest fallback-value test commit and is not evidence for the current HEAD.
