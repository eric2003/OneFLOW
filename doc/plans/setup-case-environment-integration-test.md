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
