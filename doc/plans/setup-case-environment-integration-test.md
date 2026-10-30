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

## Test-target and serial-runtime evidence

The `tests/main/CMakeLists.txt` target split is deliberate:

- `simu_context_test` compiles `SimuContext.cpp`, `SimuTask.cpp`, and accelerator-state registry/lifecycle files. It does not compile `SimuContextEnv.cpp`, so it cannot call the production case setup method.
- `simu_context_parameter_test` compiles the database implementation, `BasicParallel.cpp`, `Parallel.cpp`, and zone state, but does not compile `ParaFile.cpp` or `SimuContextEnv.cpp`. It tests the parameter helper contracts only.
- `tests/io/CMakeLists.txt` currently has a narrow `TextFileParser_test`; there is no existing test target that links the complete control-reader path.

This means a direct test cannot be added safely by simply placing a new assertion in either existing executable. The production environment translation unit references many solver/runtime components (`SolverMap`, `SolverNameClass`, `Zone`, field managers, solver reset APIs, accelerator runtime, logging, project path handling), while `ReadControlInfo` itself pulls in script parsing, configuration commit, database serialization, logging, and parallel facilities. A new target would need a coherent production link set, not an improvised subset of source files.

In `Parallel.cpp`, `Parallel::IsServer()` compares `GetPid()` and `GetServerid()`. The serial implementation in `BasicParallel.cpp` makes `HXRank()` return 0 and `HXSize()` return 1, but those functions are only applied to the `Parallel` fields by `SetUpParallelEnvironment()`. `HXBcast(callbacks, root)` returns early when `Parallel::GetNProc() <= 1`; `TestSayHelloFromEveryProcess()` still calls `CollectString()`, which consults the stored process count. Static-zero defaults make a serial path appear likely, but relying on implicit initialization is not a sound test contract. The test should establish a valid serial parallel state through the existing setup function or a narrow, explicit serial test fixture before invoking the production reader.

## Smallest next experiment

1. Build a focused integration target with the real `SimuContextEnv.cpp` and `ParaFile.cpp` production sources plus their actual target dependencies, or identify an existing library target that already supplies that complete path. Do not duplicate the parser, broadcast, or database commit logic.
2. Establish serial runtime state explicitly before calling `SetupCaseEnvironment()`. Avoid accelerator initialization unless the selected link/test path genuinely requires it; if the environment class's full translation unit makes that unavoidable, consider extracting the case-configuration operation behind a narrow production dependency boundary rather than linking the whole solver stack into a unit test.
3. Use a temporary multi-file case fixture and assert real database values and `log/database.log` creation.
4. Run two contexts sequentially. The first case defines explicit roots and a legacy source-grid fallback; the second omits them. Assert the second case does not inherit the first case's roots.
5. Ensure teardown after success and an intentionally missing control file, preserving the original setup exception.

Do not introduce a test-only parser. If an injectable seam becomes necessary, it should represent the control-loading operation, and the production implementation must continue to call `ReadControlInfo(caseDir)`. The seam should not spread across all of `SimuContext` or become a generic dependency-injection framework.

## Current outcome

**R1 remains Partial.** The global database is statically constructed, making a serial integration fixture plausible, but parallel-runtime validity and test-target linkage still need verification. No tests were run in this session. The reported 700/700 local pass predates the latest fallback-value test commit and is not evidence for the current HEAD.

## Follow-up: executable source graph and cleanup boundary

The root `codes/CMakeLists.txt` constructs `PRJ_SRC_LIST` by calling `ConstructSolutionDirTree(...)` and compiles that source set directly into the `OneFLOW` executable. It does not expose the production runtime as a reusable library target. A test that links the full source list would therefore need to recreate the executable's compile definitions, include paths, third-party libraries, and platform options, while excluding the production entry point to avoid competing with GoogleTest's main. That is a poor first R1 test target and would duplicate build configuration.

The real exception boundary is in `SimuImp::RunCase()`: it calls `PreProcess()` and `MainProcess()` inside a `try`; on any exception it calls `PostProcess()` in a nested `try/catch`, then rethrows the original exception. `PostProcess()` calls `SimuContext::TeardownCase()`. Therefore a correct failure-cleanup assertion should exercise this boundary (or the smallest extracted production boundary preserving it), not only call `SetupCaseEnvironment()` directly and inspect its exception.

### Updated implementation decision

Do not create a test executable by copying the complete solver source list into `tests/main/CMakeLists.txt`. First extract a narrow, production-owned case-setup operation with explicit dependencies for the control-load step and path binding; keep the default production implementation wired to `ReadControlInfo(caseDir)`. Then test the real control reader separately with a small multi-file fixture, and test the case orchestration's success/failure cleanup without initializing the entire solver stack. This is a deliberate boundary extraction, not a parser replacement or a generic dependency-injection framework.

This follow-up is still investigation, not completion of R1. No source code was changed and no tests were run in this step.

## Follow-up: use the repository's real CFD cases without turning them into unit fixtures

The repository's `test/` directory and `tests/` directory serve different purposes and must remain distinct:

- `tests/` contains automated unit and module tests. Keep small, deterministic contract tests here.
- `test/` contains real CFD cases and their control inputs. Use these to validate production initialization and, at a later stage, numerical regression.
- Do not copy or edit a real CFD case merely to make a test pass. If a case needs a small integration fixture, create a temporary case under the test build/runtime directory and leave the reference case unchanged.

Use three levels of evidence, in this order:

1. **Production control-reader integration:** a temporary minimal case with the same multi-file `script/control.txt` shape as real cases. Call the actual `ReadControlInfo(caseDir)`; assert typed database values and the generated `log/database.log`. This isolates file loading and database commit from the mesh and solver.
2. **Case-environment initialization smoke test:** after the production-owned setup boundary exists, select one existing case from `test/` based on its actual required files and initialization cost. Exercise setup and teardown only; do not execute the CFD task. Record the exact case path and prerequisites rather than guessing the directory name.
3. **Numerical regression:** run the established case workflow and compare its defined reference quantities. This is a separate and more expensive check; successful initialization alone is not evidence of numerical correctness.

For sequential-case isolation, the strongest small fixture is two cases run in one process: case A sets explicit project and mesh/flow roots, while case B omits them. Assert that B resolves according to its own documented defaults and never inherits A's settings. Also test a missing or malformed control input and verify cleanup preserves the original exception.

## Sandbox capability and validation boundary

The current execution sandbox provides Linux, CMake, GCC, Clang, Python, and Git. It does not have the OneFLOW checkout mounted, and outbound DNS/network access from the shell is unavailable. Repository source can be inspected and changed through the GitHub integration, but that alone does not allow a truthful local build or CTest run against the complete source tree. Consequently:

- Do not report the user-supplied 967/967 baseline as validation of new changes.
- Do not claim the real-case integration test passed until it has run against an actual checkout.
- Prefer a source change whose dependency boundary is clear from inspection, then validate it in the repository's normal Windows/CI environment. If a checkout/archive is made available in the sandbox, use it to build and run the narrow target before broadening the test.

The next implementation slice remains the narrow production-owned setup boundary described above. It should preserve the real `ReadControlInfo(caseDir)` implementation and the existing exception/cleanup semantics. Extracting a reusable production runtime library is a possible later architecture improvement, but should not be smuggled into this first integration-test change.
## Follow-up: exact checked-in case paths

Source inspection confirms two suitable existing cases under `test/`:

- `test/grid/rae2822convert/script/control.txt` lists `simu.txt`, `dimension.txt`, `tolerence.txt`, `grid.txt`, and `part.txt`. This is the smaller candidate for a real-reader/configuration smoke test, subject to checking every referenced file and any required relative assets.
- `test/cfd/rae2822_roe_sa/script/control.txt` lists `simu.txt`, `dimension.txt`, `tolerence.txt`, `cfd.txt`, `vismodel.txt`, `multigrid.txt`, `inflow.txt`, `turb.txt`, and `walldist.txt`. It is a more representative CFD control case but has a broader parameter set and likely more external assets.

The CFD case's `simu.txt` explicitly sets `string simutask = "Solve";`. That value should be asserted as a parsed database value, not executed by the setup-only smoke test. The grid conversion case is the better first candidate if the test needs to minimize dependencies; the CFD case should follow when the first reader integration is stable.

These paths are verified from the current `hexin` file contents. The sandbox still has no checkout, so the referenced directory trees have not been copied, built, or executed here.


## Follow-up: distinguish setup orchestration from cleanup ownership

Rechecked against the current production source on `hexin` (HEAD `4b41404e8325a667c4a4b39398b54a74071417b`).

The cleanup contract has two different production callers:

- `SimuImp::RunCase()` wraps `PreProcess()` and `MainProcess()` in `ExecuteWithCleanup`; its cleanup calls `PostProcess() -> SimuContext::TeardownCase()`. The outer `SimuImp::Run()` independently finalizes the process runtime.
- The compatibility entry point `SimuContext::SetupEnvironment()` wraps `SetupCaseEnvironment()` in `ExecuteWithRollback`; on failure its rollback calls `TeardownEnvironment()`, which tears down the case and finalizes the process runtime.

Do not combine these two policies into one new generic lifecycle abstraction. They have different ownership scopes: the normal multi-case runner must clean up one case while retaining process runtime for the next case; the compatibility single-case entry point rolls back both case and process setup.

The first production-owned extraction should therefore be limited to the case configuration stage inside `SetupCaseEnvironment()`:

1. Keep binding `caseDir_` and `logFile` before the real `ReadControlInfo(caseDir_)` call.
2. Keep the production reader call explicit and unchanged.
3. Extract only the post-read operation that snapshots `DataPara`, applies the configured roots, and writes the path diagnostics into a named operation with a narrow contract.
4. Test that operation using the real `Prj` path API and a controlled database; separately test the real reader with a minimal multi-file case.
5. Only after the boundary and link dependencies are understood, add the end-to-end setup/teardown test. Do not make a unit test impersonate the full solver runtime.

This is intentionally smaller than extracting a reusable runtime library. It improves observability of the existing production sequence without changing which layer owns cleanup.

## Current validation record

- User-reported baseline before this slice: 969/969 tests passed.
- Current inspected branch HEAD: `4b41404e8325a667c4a4b39398b54a74071417b`.
- This follow-up is a source-backed design refinement; it does not claim a new build, test, or CI result.


## Follow-up: extract the post-reader I/O application operation

The case-I/O application step is now a named production operation in
`codes/main/include/SimuContextIO.h`:

`SimuContextDetail::ApplyCaseIOSettings(const CaseIOSettings&)`

`SimuContext::SetupCaseEnvironment()` still owns the ordering contract:

1. Bind the selected case root and log directory.
2. Call the real `ReadControlInfo(caseDir)` implementation.
3. Read and validate the complete `CaseIOSettings` snapshot.
4. Call `ApplyCaseIOSettings()` to mutate the legacy `Prj` roots and write diagnostics.

The extracted operation deliberately retains the existing `Prj` compatibility APIs. It does not parse control files, initialize the parallel runtime, or own case/process cleanup. The focused parameter test target now links the same path implementation sources used by `prj_cmdline_test` and adds filesystem-backed checks for the default diagnostic location and explicit diagnostic names under distinct mesh/flow output roots.

This is a first production-boundary slice, not proof of the full reader integration. The new target and tests still require a real Windows/MSVC build and CTest run. Next, inspect that result before extracting anything else; then proceed to a minimal temporary multi-file fixture that calls the actual `ReadControlInfo(caseDir)` only after the serial parallel state is explicitly initialized.


## Follow-up: keep the production operation out of the header

The case-I/O application operation is now declared in
`codes/main/include/SimuContextIO.h` and implemented in
`codes/main/src/SimuContextIO.cpp`. The focused test target compiles that same
production source file. This preserves a narrow testable operation without putting
`Prj` and filesystem implementation dependencies into every translation unit that
includes the declaration. The production source collection is recursive, so the new
main-source file is included by the existing `codes` build configuration.

This is an architectural cleanup, not a new integration-evidence tier. The current
focused tests still do not call the real `ReadControlInfo(caseDir)`; R1 remains
partial. Before attempting to link the entire reader into a narrow test target, inspect
the serial parallel setup and its transitive dependencies. Avoid manually reproducing
the solver's full source list just to get a test binary.


## Follow-up: exact-HEAD CI and reader execution contract

Inspected branch HEAD: `d080aaa1b6a567ae3d0f2f5034f9332600c3e7ec`.

The Windows, Linux, and documentation workflows for this exact commit all completed
successfully:

- [Windows CI](https://github.com/eric2003/OneFLOW/actions/runs/40954345525)
- [Linux CI](https://github.com/eric2003/OneFLOW/actions/runs/40954345454)
- [Documentation check](https://github.com/eric2003/OneFLOW/actions/runs/40954345433)

This validates the current extracted `ApplyCaseIOSettings()` implementation and
its focused tests at that revision. It does not validate the full real-reader
integration described above.

A fresh source inspection confirms the precise serial contract. The production
`ReadControlInfo(caseDir)` performs the following sequence:

1. Calls `Parallel::IsServer()` and, on the server, reads
   `script/control.txt` followed by each listed script through
   `ReadMultiScriptFiles()`.
2. Calls `Parallel::TestSayHelloFromEveryProcess()` even in a one-process run.
   This path calls `HXGetProcessorName()` and then `CollectString()`; it is
   therefore not safe to assume that zero-initialized static parallel fields
   constitute a valid serial runtime.
3. Calls `BroadcastControlParameterToAllProcessors()`. Its callback-based
   `HXBcast()` returns without a collective when `Parallel::GetNProc() <= 1`.
4. Writes the merged parameter database to `log/database.log` under the supplied
   case directory.

Consequently, a reader-only integration test must explicitly call the existing
`SetUpParallelEnvironment()` before invoking the real reader. That operation calls
`HXInit()` and initializes the stored rank, process count, and server rank from
the selected parallel backend. The test must not replace the reader, broadcast,
configuration loader, or database commit with a test-only implementation.

The next experiment should be a dedicated reader integration target using the
production `ParaFile.cpp` and its actual parser/configuration/database/parallel
dependencies, plus a temporary multi-file case. It should verify typed database
values and `log/database.log`, then clear the database and temporary files. Keep
this separate from `SimuContextEnv.cpp`: that translation unit references solver,
field-manager, accelerator, and registry teardown APIs. The full case lifecycle
test should be a later step once its production dependency boundary is made
coherent. Do not copy the executable's complete source list into either test.

The user-reported local baseline remains 971/971 before these changes. The exact-HEAD
CI results above are separate evidence and must not be restated as a local test run.
R1 remains **Partial**: real reader invocation, sequential-case isolation through
the production case lifecycle, and setup-failure cleanup are still open.


## Follow-up: translation-unit boundary discovered on 2026-10-30

Inspected the exact current `hexin` source at
`69a03cc27dd7ed7ab78c7516d4265138d00b20e8`.

The serial initialization contract is confirmed in `codes/parallel/src/Parallel.cpp`:
`SetUpParallelEnvironment()` calls `Parallel::SetDefaultParallelParameter()`,
which calls `HXInit()` and initializes server rank, process rank, and process
count from the selected backend. In a serial build, `HXRank()` is 0 and
`HXSize()` is 1. This remains the correct setup for a real reader test; static
zero-initialized members are not a substitute for an explicit runtime contract.

A representative minimal real case exists at `grid/minimal_rectangle`.
Its `script/control.txt` lists five files (`simu.txt`, `dimension.txt`,
`tolerence.txt`, `grid.txt`, `part.txt`); the first files establish
`simutask = "Grid"` and `dimension = 2`. For an isolated test, create a
smaller temporary fixture with the same manifest structure rather than copying
the complete grid case and its geometry assets.

A key new dependency finding changes the implementation approach:
`ReadControlInfo(caseDir)`, `ReadPrjScript(caseDir)`,
`ReadScriptFileNameList(caseDir)`, `ReadMultiScriptFiles()`,
database compression/decompression, and database dumping currently live in
`codes/io/src/ParaFile.cpp` together with legacy parameter parsing and several
unrelated helper functions. Linking that translation unit pulls unresolved
references from its other functions too. The test cannot be made narrow merely
by adding `ParaFile.cpp` to `tests/main/CMakeLists.txt`; that would either
require an unnecessarily broad set of production sources or duplicate the
executable's dependency surface.

### Revised next implementation slice

Split by responsibility, not by test convenience:

1. Move the control-reader orchestration functions
   (`ReadControlInfo`, `ReadPrjScript`, `ReadScriptFileNameList`,
   `ReadMultiScriptFiles`) into a dedicated production source file.
2. Move the database transport/dump functions used by the reader
   (`BroadcastControlParameterToAllProcessors`, `CompressData`,
   `DecompressData`, and case-local database dump) into a narrowly named
   database/control-I/O source file if their ownership and includes support it.
   Keep the public declarations stable in `ParaFile.h` initially.
3. Leave legacy tokenization and parameter convenience helpers in
   `ParaFile.cpp`; avoid changing parser behavior in this slice.
4. Add the reader integration target using the exact extracted production source,
   the real `ConfigLoader`/adapter, database, parallel, and file/log dependencies.
   Do not mock parsing, database commit, or serial initialization.
5. Use a temporary multi-file manifest; assert representative typed values and
   `log/database.log`. Then clear the database and temporary directory.
6. Only after reader integration passes, add the two-case sequence through the
   case-owned setup/teardown boundary. Keep process runtime initialization and
   finalization outside each per-case iteration.

The architectural benefit is a reusable control-file ingestion boundary and an
explicit dependency graph. It does not change the numerical model or the meaning
of any control parameter. The current parameter pipeline is:

```text
legacy text files
    -> ConfigLoader
    -> ConfigDocument
    -> ConfigDatabaseAdapter
    -> global DataPara
    -> parallel broadcast
    -> case database diagnostic
```

The parser already has a useful separation (`ConfigLoader` creates a
`ConfigDocument`, and `ConfigDatabaseAdapter` commits typed entries). The
remaining coupling is the orchestration/transport code co-located with legacy
helpers. The refactor should preserve this pipeline and make its stages independently
linkable and testable.

Validation status at this point:
- User-reported local baseline: 971/971 before these changes.
- Documentation check for `69a03cc`: successful.
- Windows and Linux CI for `69a03cc`: in progress at last inspection.
- No source-code or CMake changes were made in this follow-up; R1 remains Partial.


## Implementation update: source responsibilities separated

Base inspected: `hexin` HEAD `f20725d1d6b5721393273ffc4ccaff9d78b30bed`.

The first production source split is now present:

- `codes/io/src/ControlFileReader.cpp` owns both `ReadControlInfo`
  overloads, `ReadPrjScript`, manifest enumeration, and ordered multi-script
  orchestration.
- `codes/io/src/ControlDatabaseIO.cpp` owns database broadcast, compression,
  decompression, and database dump operations.
- `codes/io/src/ParaFile.cpp` retains legacy syntax and parameter helpers.
- `codes/io/include/ParaFile.h` keeps the existing public declarations.

The extraction preserves function bodies and call order; it is intended to be
behavior-preserving. No parser semantics, parameter types, broadcast policy,
or lifecycle ownership were intentionally changed.

### Why this is more than moving functions

The original translation unit mixed three responsibilities: legacy syntax
processing, orchestration of a multi-file case manifest, and transport/diagnostic
operations on the runtime database. That made the dependency graph wider than
the conceptual job being tested. The new files establish two named boundaries:
the reader coordinates the ingestion stages; the database I/O module owns
cross-process transfer and diagnostic persistence. The existing
`ConfigLoader -> ConfigDocument -> ConfigDatabaseAdapter -> DataPara` chain
remains the semantic path.

The source collector `ConstructSolutionDirTree()` uses recursive source
discovery, so the new C++ files are expected to be included automatically.
This is a source-inspection finding, not a substitute for a clean build.

### Expected issues and test purpose

1. **Compile/source-discovery issue:** a new file might not enter the build,
   or may expose missing direct includes. A clean Windows and Linux build must
   establish the actual compilation boundary.
2. **Remaining parser linkage:** `ReadMultiScriptFiles()` still calls
   `ReadOneFLOWScriptFile()` in `ParaFile.cpp`. The reader is therefore not
   yet proven independently linkable. The next step is to inspect the parser
   adapter's dependencies and extract only that adapter if it yields a useful
   production boundary.
3. **Serial runtime issue:** `ReadControlInfo()` calls parallel coordination
   even for one process. A real test must initialize `SetUpParallelEnvironment()`;
   zero-initialized static fields are not a valid substitute.
4. **Behavior regression:** the ordered manifest must still load every file,
   commit declared types, broadcast the database, and write
   `log/database.log`. A temporary fixture should assert representative typed
   values and the diagnostic output.
5. **Cross-case contamination:** case A's database and roots must not leak into
   case B. A two-case test in one process is required.
6. **Failure cleanup:** a missing or malformed control input must preserve the
   original exception while leaving the next case able to initialize cleanly.

The tests therefore serve distinct purposes: source/build checks catch module
integration errors; the reader integration test pins the real ingestion
contract; the two-case test pins state isolation; and the lifecycle failure test
pins cleanup ownership. None of these tests establishes CFD numerical validity.

### Architectural evolution model

The target is a pipeline of small, explicit transformations rather than a
single monolithic reader:

```text
legacy files
  -> syntax adapter
  -> ConfigDocument
  -> typed parameter adapter
  -> DataPara snapshot
  -> distribution adapter
  -> diagnostic serializer
```

Each stage should declare its input/output representation, invariants, units
where applicable, ownership, side effects, lifecycle scope, and failure policy.
The proposed machine-readable starting point is
`doc/specs/computation-unit.schema.json`; the design discussion is in
`doc/source/evolutionary_computation_architecture.rst`. These are architectural
proposals, not a claim that a plugin runtime already exists.

Short term, stabilize the real control reader and its test seam. Medium term,
apply versioned intermediate representations and directional adapters to mesh,
topology, boundary, and result data. Long term, compose alternative solver
algorithms and dataflow stages against shared mathematical and data contracts.
Dynamic plugin loading is deferred until an actual extension requirement justifies
its cost.

### Validation status

- The source extraction is **implemented, unverified**.
- No local build or CTest run was performed in this environment.
- User-reported 971/971 is the pre-change baseline only.
- R1 remains **Partial**. Do not infer successful integration from source separation.


## Second extraction: parser adapter and array-size resolver

The next source-level step is implemented on `hexin`:

- `ControlScriptParser.cpp` now owns both `ReadOneFLOWScriptFile` overloads.
  It connects `TextFileParser` / file input to `ConfigLoader`, then commits
  the resulting `ConfigDocument` through `ConfigDatabaseAdapter`.
- `ControlParameterResolver.cpp` owns `GetParameterArraySize()`. Literal
  lengths are parsed locally; named lengths are looked up in the runtime
  `DataBase`.
- `ParaFile.cpp` now retains legacy parameter inspection and extraction
  helpers. Direct dependencies on `Word.h` and `FileUtils.h` are explicit.

This distinction is important: `ConfigLoader` already accepts an injected
array-size resolver. It does not intrinsically require the global database.
The production adapter supplies the legacy database-backed resolver, while
unit tests can supply a deterministic callback. That is a genuine seam for
testing parser behavior without pretending the runtime database is absent in
production.

The current dependency graph is:

```text
ControlFileReader
  -> ControlScriptParser
       -> ConfigLoader
            -> ConfigDocument
       -> ConfigDatabaseAdapter
            -> DataPara
       -> ControlParameterResolver
            -> DataBase (named array lengths)

ControlDatabaseIO
  -> DataBase + Parallel + diagnostic file output

ParaFile
  -> legacy inspection helpers
```

### What this enables and what it does not prove

- The source-level parser adapter is no longer defined in the broad legacy
  `ParaFile.cpp` translation unit.
- The resolver has an explicit runtime-state dependency.
- Existing `ConfigLoader_test` cases pin scalar parsing, literal arrays,
  injected named-array lengths, missing-resolver rejection, and document
  replacement. A new `control_parameter_resolver_test` target pins literal
  lengths, real named-length lookup through DataBase, and missing-name failure.
- These tests verify parser semantics, not the production commit to DataPara,
  the multi-file manifest, parallel broadcast, or database dump.
- A real reader integration target is still required to prove those runtime
  behaviors together.

The next test design should avoid two extremes: do not link the complete
solver executable just to parse a text fixture, but do not mock away the
database commit and then call that a production ingestion test. Keep a
parser-level test with an injected resolver, then a component/integration test
that links the real adapter, real DataBase implementation, and real serial
parallel setup only where the contract requires them.

### Expected failure modes

1. **Resolver semantics drift:** literal and named lengths could be interpreted
   differently after extraction. Keep explicit tests for both.
2. **Missing named length:** the runtime lookup should retain the existing
   database exception and identify the missing variable; do not silently
   substitute zero.
3. **Cross-file ordering:** a later script may use a length declared by an
   earlier script. The reader must preserve manifest order; do not parallelize
   script parsing unless dependencies are explicitly modeled.
4. **Partial commit:** if a later entry fails, earlier entries may already have
   been committed. This is an important transaction question, not something
   this extraction solves. Define and test transaction semantics before
   promising all-or-nothing control-file loading.
5. **Repeated declarations:** DataPara currently updates existing entries
   under its own type/size compatibility rules. The adapter must not introduce
   a second, conflicting policy.
6. **Diagnostic failure:** database dump failure occurs after the parameter
   commit and broadcast. Decide whether diagnostics are mandatory for a
   successful case setup, and keep that policy explicit.

### Architectural question opened by this stage: should configuration commit be transactional?

The current pipeline has a meaningful distinction between a parsed
`ConfigDocument` and runtime `DataPara`. This makes a future transactional
commit possible: validate all names, types, lengths, and duplicate declarations
against a snapshot before mutating live state, then apply as one controlled
operation or restore the previous snapshot on failure.

Do not implement that transaction speculatively in this extraction. First
characterize the current behavior with tests, especially malformed late entries
and duplicate declarations. Then decide whether atomic commit is a required
scientific reproducibility guarantee or whether documented incremental commit
is sufficient. This decision affects memory cost, error reporting, compatibility,
and how configuration reload works.

### Validation state

- Parser adapter and resolver extraction: **implemented, unverified**.
- Resolver tests added: literal lengths, database-backed named lengths, and
  missing-name exception. The target links the real database sources; execution
  awaits exact-revision CI.
- ConfigLoader unit tests remain the separate proof of injected-resolver parser
  semantics; no new test run is claimed.
- New source files rely on recursive CMake discovery with `CONFIGURE_DEPENDS`;
  exact-revision build evidence is still required.
- The 971/971 user baseline predates this work.


## 2026-10-30 update: typed adapter coverage added

The production configuration path is now split across these current source units:

- `ControlFileReader.cpp`: case manifest reading, ordered script-file orchestration, serial/process coordination, and reader entry points.
- `ControlScriptParser.cpp`: legacy script parsing into `ConfigDocument` followed by the production `ConfigDatabaseAdapter::Commit()`.
- `ControlParameterResolver.cpp`: literal or runtime-DataBase-backed array-size resolution.
- `ControlDatabaseIO.cpp`: database broadcast serialization and diagnostic dump.
- `ParaFile.cpp`: remaining legacy parameter-inspection and compatibility helpers.

The `control_parameter_resolver_test` target now compiles the production
`ConfigDatabaseAdapter.cpp` alongside the real database implementation. Its
new tests assert that integer, Real, and string entries are committed with
their runtime types and that repeated declarations are applied in document
order, with the later value replacing the earlier one.

This establishes a useful intermediate contract:

```text
ConfigDocument (ordered, format-neutral entries)
          |
          v
ConfigDatabaseAdapter
          |
          v
DataBaseType lookup + ProcessData conversion
          |
          v
DataPara replacement / ownership
          |
          v
runtime parameter lookup
```

The adapter is incremental rather than transactional. If a later entry fails,
earlier entries may already have changed the runtime database. This behavior
is now documented explicitly; whether configuration commits should become
transactional is an architectural decision that needs compatibility analysis
and a dedicated failure test before changing semantics.

**Evidence boundary:** these tests cover the production document-to-database
adapter and the existing array-size resolver. They do not prove physical
manifest order across multiple files, `ReadControlInfo(caseDir)` orchestration,
parallel broadcast, `log/database.log`, effective path application, or
sequential-case isolation. R1 remains Partial.

## Next integration step

Do not broaden the test target by copying the full executable source list.
First map the exact unresolved symbols pulled in by `ControlFileReader.cpp`,
`ControlDatabaseIO.cpp`, `ControlScriptParser.cpp`, `TextFileParser.cpp`,
and the project/logging helpers. Then decide whether the smallest sound
production boundary is:

1. a small configuration-runtime library target shared by the executable and
   tests; or
2. a narrower extracted production operation with explicit filesystem,
   database, and parallel dependencies.

Prefer the first only if the source graph is coherent and reusable by the
production executable. Avoid a test-only parser or duplicate implementation.
The integration fixture must initialize serial parallel state through
`SetUpParallelEnvironment()`, use at least two real script files, assert
typed runtime values, and verify the database diagnostic output.

After that, extend to `SetupCaseEnvironment()`, path diagnostics, two
sequential cases, and cleanup after a missing/invalid control file. Preserve
the separation between case teardown and process finalization.

## Documentation ownership decision

Keep implementation plans in `doc/plans/`, separate from the published
Sphinx knowledge base in `doc/source/`. A plan is allowed to contain
hypotheses, sequencing, risks, acceptance criteria, and unresolved questions.
The source manual should describe the current architecture and verified
behavior, mathematical/algorithmic explanations, and stable design contracts.
When a plan item is completed, promote its durable conclusions into the
relevant `doc/source` page and retain only useful decision history and
remaining actions in the plan. Link the two; do not duplicate them wholesale.

## 2026-10-30 follow-up: production dependency boundary inspection

Inspection of the current `hexin` source confirms that `ReadControlInfo(caseDir)`
executes the following sequence:

1. On the server process, read the ordered `script/control.txt` manifest and
   invoke `ReadOneFLOWScriptFile` for each listed file.
2. Run `Parallel::TestSayHelloFromEveryProcess()`, which calls
   `HXGetProcessorName()` and `CollectString()` even in a one-process run.
3. Call `BroadcastControlParameterToAllProcessors()`; its callback-based
   `HXBcast` returns early for one process, but only after parallel state has
   been initialized.
4. Dump the resulting runtime database to `log/database.log` under the case.

The serial setup contract is `SetUpParallelEnvironment()` ->
`Parallel::SetDefaultParallelParameter()` -> `HXInit()` / `HXRank()` /
`HXSize()`. The integration test must use this production entry point rather
than depending on static-zero values of `Parallel::pid` and `Parallel::nProc`.

The CMake structure also matters: `codes/CMakeLists.txt` recursively collects
production sources into the `OneFLOW` executable and applies project-wide
compile definitions, include paths, and external libraries there. The IO tests
are separately assembled from selected translation units. There is currently
no reusable configuration-runtime library target. This is a real architectural
boundary, not merely a test inconvenience: a shared target would need a coherent
public dependency surface and matching compile definitions, not a copied list
of arbitrary sources.

The smallest manifest fixture already present in the repository is
`test/grid/rae2822convert`; its control manifest lists five files, with
`simu.txt` declaring the Grid task and `dimension.txt` declaring dimension 2.
A setup-only integration test can use a temporary copy of this fixture and
must not run the grid conversion task itself.

Decision for the next implementation slice: first build a dependency map from
`ControlFileReader.cpp` through parser, database, project-file, logging, and
parallel support. Extract a shared production target only if that map shows a
stable cohesive subsystem. Otherwise prefer a narrow production operation with
explicit dependencies. Do not create a test-only imitation of the reader.

## 2026-10-30 follow-up: isolate and test manifest discovery

The two `ReadScriptFileNameList` overloads have been moved from
`ControlFileReader.cpp` to `ControlManifestReader.cpp`. This is a deliberately
narrow extraction: it separates filesystem-backed manifest discovery from the
process coordination and database broadcast performed by
`ReadControlInfo(caseDir)`. The public declarations remain in `ParaFile.h`.

A new test in `tests/io/control_manifest_reader_test.cpp`, registered with the
existing `TextFileParser_test` target, creates a temporary case manifest and
asserts that blank lines are ignored, file paths are resolved under the given
case root, and the listed order is preserved. The header now documents the
append-to-existing-vector and ordering contract. The extraction does not alter
manifest syntax or the subsequent script parsing semantics.

This is useful component evidence and a clean dependency seam, but it is not
the real multi-file integration test: it does not invoke `ReadMultiScriptFiles`,
`ReadControlInfo(caseDir)`, parallel setup, DataBase commit, or database dump.
The next step is to map the dependencies of those orchestration operations and
then select a production-shared target or a narrower seam based on that graph.


## 2026-10-30 follow-up: source-level dependency graph checkpoint

The current source split makes the configuration-ingestion path more explicit.
The following production edges were verified directly at `4b2e5d496ed359c1e3c3645413c4541dc81122f7`:

| Responsibility | Production entry / implementation | Direct dependency evidence |
|---|---|---|
| Manifest discovery | `ControlManifestReader.cpp` | `TextFileParser`, `Prj`, filesystem-backed case-file resolution |
| Ordered multi-file orchestration | `ControlFileReader.cpp` | `ReadScriptFileNameList`, `ReadMultiScriptFiles`, `Parallel` |
| Script parsing | `ControlScriptParser.cpp` | `ConfigLoader`, `GetParameterArraySize`, `ConfigDatabaseAdapter` |
| Array-size resolution | `ControlParameterResolver.cpp` | real global `DataBase` lookup for named array dimensions |
| Typed document commit | `ConfigDatabaseAdapter.cpp` | `DataBaseType` and `ProcessData` |
| Runtime coordination | `ControlFileReader.cpp` | `Parallel::TestSayHelloFromEveryProcess`, callback-based `HXBcast` |
| Database diagnostic output | `ControlDatabaseIO.cpp` | `DataBase`, `DataBook`, `Prj`, `PIO`, and global `logFile` |

This graph exposes two distinct testable boundaries:

1. **Configuration ingestion**: manifest discovery -> ordered script parsing -> typed database commit. This can validate the actual production parser and database adapter without linking the parallel broadcast and diagnostic-output path.
2. **Reader orchestration**: `ReadControlInfo(caseDir)` additionally requires initialized parallel state, broadcast support, logging, and database-file output. This is a broader integration contract and should remain a separate test layer.

The existing `control_parameter_resolver_test` already links the production resolver and adapter with the real database sources, but it does not load a manifest or parse multiple physical files. The existing `TextFileParser_test` links manifest discovery, but does not parse and commit the listed files. The gap is therefore specifically the composition between these existing component boundaries; it is not a lack of another isolated parser test.

CMake evidence: `tests/io/CMakeLists.txt` assembles the two current targets independently, while `codes/CMakeLists.txt` gathers production sources into the executable. No reusable configuration library target currently owns the complete ingestion dependency set. This argues for a first focused integration executable that composes the real ingestion translation units and the already-used database support sources. It does **not** yet justify extracting a public library or moving the executable's whole source graph.

A separate caution applies to the serial runtime fixture: `SetUpParallelEnvironment()` calls `HXInit()`. Whether this invokes MPI initialization depends on the compile-time parallel configuration. The integration target must inherit the same relevant compile definitions and link requirements as the intended production mode; it must not accidentally call MPI initialization in a process whose test harness or another fixture already owns MPI lifecycle. The serial and MPI-enabled contracts should not be conflated.

### Revised next implementation slice

Add one test target for the **configuration-ingestion** boundary first:

- create a temporary case directory with `script/control.txt` and at least two real script files;
- call production `ReadScriptFileNameList(caseDir, files)` and `ReadMultiScriptFiles(files)`;
- initialize/clear the real global database using its supported lifecycle, not a substitute store;
- assert typed scalar values and an ordered duplicate-key replacement across files;
- verify a missing listed script fails through the production parser path and leaves the database state behavior explicitly characterized;
- clean the temporary directory with RAII so failed assertions do not leak fixtures.

Do not call `ReadControlInfo(caseDir)` from this first target until the test target's parallel, logging, `PIO`, and `DataBook` dependency closure is explicit. Once ingestion is verified, add the orchestration test as a separate layer for serial runtime initialization, broadcast no-op semantics at one process, and `log/database.log` creation. Only then decide whether those dependencies form a reusable production library.

This narrows the next action from "link the whole reader" to a testable composition boundary while preserving the real production path. R1 remains **Partial**: no such multi-file ingestion integration target or post-change test evidence exists yet.


## 2026-10-30 implementation: production multi-file ingestion composition

Starting revision: `4b359903f4ffbf3bf6c595f066d920bf6b26e5b3`. The user reports 977/977 local tests passed before this work; that is a pre-change baseline only.

### Production boundary change

`ReadMultiScriptFiles` has moved from `ControlFileReader.cpp` into `ControlScriptParser.cpp`. The former now owns only runtime orchestration: server-side case script selection, process greeting/collection, database broadcast, and diagnostic dump. The parser translation unit owns the ordered sequence of file parsing and per-file database commit. The declaration now accepts a const vector reference and documents two important semantics:

- listed files are processed in their supplied order;
- each file is committed before the next file is read, so the whole list is not a transaction.

This is a narrow production extraction. It does not alter the `ReadControlInfo` orchestration sequence or create a parallel/runtime substitute.

### Integration test added

`tests/io/control_configuration_ingestion_test.cpp` creates a unique temporary case with a real `script/control.txt` manifest and two real script files. It invokes production `ReadScriptFileNameList(caseDir, files)` and `ReadMultiScriptFiles(files)`, then reads the actual global `DataBase` and asserts:

- integer, Real, and string values are committed with their intended types;
- a parameter declared in both files takes the value from the later manifest entry;
- values declared only in the first file remain present after the second file is processed.

The temporary case and global parameter store are cleaned up by scope-bound fixture cleanup. The test is included in the existing configuration/resolver executable, renamed from `control_parameter_resolver_test` to `control_configuration_test`; the existing resolver and adapter tests remain in that executable. CMake composes the real parser, manifest, adapter, resolver, database, and required low-level support translation units. It does not compile `ControlFileReader.cpp`, because that file intentionally remains coupled to parallel orchestration and database dumping.

### What this test does not establish

It does not call `ReadControlInfo(caseDir)`, initialize MPI/serial process state, exercise broadcast, write `log/database.log`, test `SetupCaseEnvironment()`, or prove sequential-case isolation and failure cleanup. It also has not yet been built or run in this environment. Treat the source change and test as **Implemented, unverified** until exact-HEAD Windows/Linux CI and the user's local MSVC/CTest run confirm it.

### Next validation and design step

1. Check exact-HEAD Windows and Linux CI for compile/link issues in the expanded test target.
2. Run the complete local Debug/Release CTest suite and report the post-change count separately from the 977/977 baseline.
3. Add explicit failure behavior coverage only after confirming how the production file-open/parser path reports a missing manifest entry; do not silently assume whole-list rollback. The current adapter and multi-file sequence are incremental.
4. Add a separate `ReadControlInfo(caseDir)` integration target only when its serial/MPI initialization, process coordination, broadcast, logging, PIO, and diagnostic dump dependencies are all represented faithfully.
5. Continue from reader integration to effective case path application, sequential-case isolation, and cleanup at the actual `RunCase` exception boundary.

R1 remains **Partial**. This is the first test that composes the real manifest reader, production parser, typed adapter, and global database; it is not yet a full runtime integration test.


## 2026-10-30 follow-up: database API contract documentation

The database API now has source-level documentation for type-erased storage, typed access preconditions, const/mutable read symmetry, ownership, and borrowed lookup lifetime. A Sphinx page `doc/source/database_architecture.rst` describes the `DataBase -> DataPara -> DataEntry -> DataObject -> TDataObject<T>` ownership chain, same-key update behavior, field pointer registry distinction, serialization boundary, and the tested configuration-ingestion path. The ingestion test now exercises both mutable and const scalar reads. User-reported 978/978 is the pre-change baseline for this documentation/API slice; revalidate the new exact HEAD before treating it as passing. Full `ReadControlInfo` and Case lifecycle coverage remains open.


## 2026-10-30 continuation: configuration-ingestion failure contracts

Starting revision: `3fadcbacafee01081b65bc889218045566dc6a46`. The user reports
978/978 local tests passed before this continuation; that is a pre-change
baseline, not validation of the new tests.

Source inspection established the failure boundaries:

- `ReadMultiScriptFiles` processes the ordered path list sequentially.
- Each `ReadOneFLOWScriptFile(file)` constructs a fresh `ConfigLoader`,
  parses the complete file into a `ConfigDocument`, and only then invokes
  `ConfigDatabaseAdapter::Commit`.
- A missing file fails in `Prj::OpenFile` via the existing
  `Fatal(...)` / `std::runtime_error` path.
- A parse-time failure (for example, a named array length that cannot be
  resolved) occurs before the current file's document is committed.
- Earlier files have already been committed and are not rolled back. The
  adapter itself also documents that a commit is not atomic if a later entry
  fails validation or conversion.

Added two production-composition tests to
`control_configuration_ingestion_test.cpp`:

1. A valid first file followed by a missing file must throw while retaining
   the first file's committed value.
2. A valid first file followed by a file whose later named-array length cannot
   be resolved must throw; the earlier file remains committed, while even a
   valid declaration parsed before the failing declaration in the current
   file is not committed.

These tests define the current incremental contract rather than imposing
whole-manifest rollback. They also distinguish parse-before-commit at the
file boundary from the adapter's non-transactional per-entry commit behavior.
The latter still merits a separate contract test if realistic conversion or
type/size mismatch failure can be induced without corrupting global test state.

No local build or CTest was run here. The 978/978 report predates these edits.
Exact starting-HEAD CI was still running on Windows and Linux; Sphinx had
passed. Recheck CI after the new test commit. R1 remains Partial: these tests
do not yet exercise `ReadControlInfo(caseDir)`, parallel coordination,
broadcast, diagnostic database output, Case path application, sequential-case
isolation, or setup-failure cleanup.


### Follow-up: adapter commit failure within one document

A third failure-path test was added after reviewing the adapter's documented non-atomic behavior. The fixture first commits an integer key, then parses a second file containing a new integer entry followed by a declaration that reuses the existing key with a different type. The second file parses completely, but its adapter commit throws on the type mismatch. The test asserts that the original key remains unchanged and that the earlier new entry from the second document remains committed. This distinguishes three boundaries: file-open failure, parse failure before commit, and adapter failure partway through commit. It confirms that neither the whole manifest nor an individual document commit is currently transactional. The test source newline literals were corrected in commit `971f14e40ad4933fba429350d71b485d8b4be3a8`.

The new test is source-level only until exact-HEAD CI and the user's local test run complete. Do not change the current behavior to rollback without an explicit compatibility decision; transactional commit would require staging/validation of all entries before mutating the global DataBase.


## 2026-10-30 continuation: separate script ingestion from runtime orchestration

The production source previously placed two different responsibilities in
`ControlFileReader.cpp`:

- `ReadPrjScript` composes ordered manifest discovery with parsing and
  database commit.
- `ReadControlInfo` adds process coordination, broadcast, and database
  diagnostics around that script reader.

This made the useful case-level reader composition unavailable to the focused
configuration test without also linking the runtime orchestration and its
parallel/PIO/logging dependencies.

The responsibilities are now separated:

- `codes/io/src/ControlScriptReader.cpp` owns the two `ReadPrjScript`
  overloads. Its dependency is the public control-reader interface and the
  manifest/multi-file reader functions.
- `codes/io/src/ControlFileReader.cpp` retains `ReadControlInfo` and its
  parallel/broadcast/diagnostic orchestration.
- The main executable's recursive source collection includes the new
  translation unit. The focused `control_configuration_test` explicitly
  compiles `ControlScriptReader.cpp`, and its end-to-end ingestion test now
  calls `ReadPrjScript(caseDir)` rather than manually composing
  `ReadScriptFileNameList` and `ReadMultiScriptFiles`.

This is a dependency-boundary improvement, not merely a test convenience:
manifest ingestion can now be linked and exercised independently of
process-level orchestration. It preserves the existing public function
signatures and ordered, non-transactional commit behavior.

Acceptance boundary: the test still does not call `ReadControlInfo(caseDir)`.
That next step must model initialized serial/MPI state and real database
diagnostics; do not replace it with mocks or copy the full executable source
list. The source split and new integration call are implemented but await
exact-HEAD CI and the user's local CTest result.


## 2026-10-30 continuation: reader orchestration and parallel file ownership

The user reports 981/981 tests passed before this continuation. That result is
the starting baseline, not evidence for the changes below. Starting revision:
`144cde3f472fc73d2a8633f84b5eaea874308de8` on `hexin`.

A source-level dependency review found that the case-specific database dump
did not use the parallel file-ownership boundary consistently:

- `DumpDataBase()` opens its output through `PIO::OpenPrjFile`, which returns
  immediately on non-owner ranks.
- `DumpDataBase(caseDir)` previously called `Prj::OpenCaseFile` directly,
  opening the same case-relative output on every rank, but closed it through
  `PIO::CloseFile`, which is owner-filtered.
- `Prj::OpenCaseFile` only resolves the case path, creates the parent
  directory, and opens the stream; it does not perform rank ownership checks.

The case overload now uses `PIO::OpenCaseFile` as well. This restores the
same ownership rule for both database dump entry points and avoids concurrent
multi-rank opens/truncation of `log/database.log`.

To make the explicit path operations independently linkable without dragging
the task-state-only PIO entry points into a focused runtime integration test,
the explicit file operations (`ParallelOpen`, explicit `OpenPrjFile`,
`OpenCaseFile`, and `CloseFile(fstream&)`) were moved into
`codes/io/src/PIOFileAccess.cpp`. `PIO.cpp` now retains the implicit
task-bound `OpenPrjFile()` and `CloseFile()` operations. The public PIO API
is unchanged; production recursive source collection discovers the new
translation unit.

The `control_configuration_test` target now links the production
`ControlFileReader.cpp`, `ControlDatabaseIO.cpp`, `PIOFileAccess.cpp`,
and `LogFile.cpp`. A new serial-runtime integration test calls
`SetUpParallelEnvironment()` and then the real
`ReadControlInfo(caseDir)` entry point. It verifies:
- two-file manifest ingestion and later-file override;
- typed values in the actual global DataBase;
- the case-relative `log/database.log` dump and its contents;
- the case-bound broadcast diagnostic log.

This test exercises the real serial process-coordination path and the
single-process broadcast short circuit. It does **not** validate MPI
collectives, multiple-rank file ownership at runtime, full
`SetupCaseEnvironment()`, sequential-case isolation, or exception cleanup.
Those remain separate R1 acceptance items. The test target is intentionally
serial and does not claim MPI coverage.

No local build or CTest was run in this environment. The user's 981/981
report predates these edits. Exact-HEAD Windows/Linux CI and a new local test
run are required before treating this slice as verified.


## 2026-10-30 continuation: map the real Case lifecycle boundary

The user reports 982/982 tests passed locally before this continuation. Treat
that as the baseline for the current work, not as evidence for any later
changes.

A fresh inspection of `Simulation.cpp`, `SimuImp.cpp`,
`SimuContextEnv.cpp`, `SimuContextIO.cpp`, and `Prj.cpp` confirms the
production multi-case lifecycle:

1. `Simulation::RunImpl` parses all case directories and initializes
   process-level project state once.
2. A dedicated process `SimuContext` owns
   `SetupProcessEnvironment` / `FinalizeEnvironment`.
3. Each directory creates a separate `SimuImp` and therefore a separate
   case `SimuContext`.
4. `SimuImp::RunCase` wraps `PreProcess` + `MainProcess` in
   `ExecuteWithCleanup`; its cleanup invokes `PostProcess`, which calls
   `TeardownCase`.
5. `PreProcess` calls `SetupCaseEnvironment`; a setup exception therefore
   still enters the same case cleanup path before the outer process lifecycle
   finalizes the shared runtime.
6. `SetupCaseEnvironment` binds `Prj::prjBaseDir`, binds the case log,
   calls real `ReadControlInfo(caseDir)`, reads the typed IO settings from
   the global DataBase, and applies them to the legacy `Prj` path state.
7. `TeardownCase` releases solver/zone/field/database state, clears the
   case log binding, and clears every case-specific `Prj` root. It does not
   finalize the process runtime; that remains owned by the process context.

The path lifecycle already has a focused contract test:
`PrjIOPaths.DomainRootsResetWhenBindingNextCase` configures project and
domain-specific roots for one case, binds a second case, and verifies that
mesh/flow input/output paths resolve under the second case. The production
`Prj::SetPrjBaseDir` resets project and domain roots, while
`Prj::ClearPrjBaseDir` clears all of them. This is useful component-level
evidence, but it is not a test of `SetupCaseEnvironment` or
`TeardownCase` as composed by `SimuImp`.

The key remaining integration gap is now narrower and explicit: the test
suite separately exercises real `ReadControlInfo(caseDir)`, IO parameter
snapshot/application, and `Prj` case-root reset, but does not invoke the
production composition in `SimuContext::SetupCaseEnvironment` or its
cleanup from `SimuImp::RunCase`.

Do not label the existing `PrjIOPaths` test as sequential-case runtime
isolation. Next, consider extracting the case bootstrap sequence into a small
production translation unit with a narrow dependency surface, then link that
same unit into a serial integration test. The extracted function must remain
the one called by `SimuContext::SetupCaseEnvironment`; do not write a parallel
test-only orchestration. The integration test should verify the actual
case-root/log binding, reader/database ingestion, effective generic and
domain IO paths, and diagnostic files. A separate higher-level test should
exercise cleanup after a failure during bootstrap. Keep process-level MPI and
accelerator initialization outside the serial case-bootstrap test, and retain
a separate MPI lifecycle validation plan.

No source change or post-baseline test result is claimed by this investigation
entry. Exact-HEAD Windows/Linux CI for the preceding PIO link fix is still
running at the time of this note.
