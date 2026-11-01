# OneFLOW Testing Architecture

## Purpose

OneFLOW has several kinds of verification activity. They share a user-facing
entry point, but they do not share the same execution semantics. The design
therefore unifies discovery, classification, command construction, and result
presentation while retaining the native authority of each test framework.

## Test layers and authorities

| Layer | Source of truth | What it establishes | How it runs |
|---|---|---|---|
| C++ unit / integration tests | GoogleTest declarations and CMake targets | Local invariants and subsystem/production-path behavior within the test scope | CTest in the configured build tree |
| Python tests | Python test modules and unittest cases | Catalog validation, GUI services, workflow logic, tooling behavior | `python -m unittest` or purpose-built Python commands |
| Repository Test Inventory | Source scan of C++ and Python test declarations | What test declarations are present and where they live | Read-only discovery; does not execute tests |
| Contract / scenario fixtures | Numerical catalog entries and Geometry scenario files | Intended formula/fixture metadata and traceability links | Contract Workbench adapters; no implicit execution |
| CFD / grid regression | Existing `test/suites/*.txt` manifests and case data | Numerical behavior of selected representative cases | Existing `test/test.py` runner |
| Development Tests | Case catalog | Selected development cases run in a protected workspace | Existing OneFLOW case runner |

CTest and Python unittest are complementary, not competing frameworks.
CTest is the CMake-level test orchestration layer and discovers registered
native test cases. Python unittest executes Python test modules. A GUI checkbox
labelled "Unit Tests (CTest and Python tests)" means the selected CTest suite
plus the Python command set registered by the Workbench.

## The three Contract Workbench Python modules

These are ordinary executable Python unit-test modules, not special test
frameworks and not CTest tests.

### `tools.verification.test_curve_catalog`

Validates the numerical theory catalog:

- expected SA, SST, and NS-LUSGS catalog families are present;
- formula LaTeX and gated code sketches are populated;
- repository-relative path resolution is portable;
- declared pure API names resolve to source symbols;
- declared source, documentation, and data paths exist.

This protects the Numerical Contract Workbench's metadata and its links to
implementation/documentation. It does not run a CFD solver or prove a complete
turbulence model.

### `tools.verification.test_contract_catalog`

Validates the shared read-only contract adapter:

- adapts Numerical TheoryEntry and Geometry scenario data;
- assigns stable namespaced IDs;
- distinguishes declared APIs, available fixtures, and linked-but-unexecuted
  tests from actual execution evidence;
- rejects missing IDs, duplicate IDs, missing files, absolute paths, and
  repository path traversal;
- detects stale linked GoogleTest declarations;
- loads the combined catalog and checks its references;
- verifies that importing Geometry scenario metadata does not eagerly load
  the optional pipeline runner.

It does not execute GoogleTest, geometry production C++, or CFD cases.

### `tools.oneflow_gui.test_test_inventory`

Tests the source-based Repository Test Inventory:

- discovers C++ GoogleTest declarations and their source locations;
- attaches input/expected metadata only when a Geometry scenario explicitly
  links the test;
- leaves result unset until a real execution backend reports one;
- discovers Python unittest methods;
- produces deterministic discovery order;
- checks that Geometry fixture metadata is not presented as an executed result.

This is a test of test discovery. It is not itself a test runner.

## Inventory versus execution

The GUI's Repository Test Inventory and test execution actions serve different purposes:

- **Repository Test Inventory** scans source files and displays declarations,
  source locations, and explicit scenario links. It is read-only. Seeing a test
  in this list means it was discovered, not that it ran or passed.
- **C++ Unit Tests (CTest)** invokes the configured CTest suite.
- **Python Unit Tests (unittest)** invokes the explicit list of Python test
  commands through the existing command worker. These two selections are
  independent and can be run separately or together. A test appears in the
  inventory if its declaration can be found; it runs only if its executable
  target or Python command is registered with the corresponding runner.
- **CFD/Grid Regression** invokes the existing manifest-driven regression
  runner. It is not interchangeable with a unit test.
- **Development Tests** execute selected catalog cases from an isolated
  workspace. They are not source-discovered unit tests.

The source inventory intentionally has broader discovery semantics than the
execution allow-list. Some Python files need special arguments, dependencies,
working directories, or environment setup. Automatically running every
`test_*.py` found by an AST scan would be unsafe and could silently change the
meaning of the repository's test command. Add executable commands explicitly
and cover the command registry with tests.

## Current GUI integration

The three Contract Workbench modules are registered in the Workbench Python
unit-test command list. Running **Tests > Python Unit Tests** runs them
independently of CTest. The C++ Unit Tests selection runs only the configured
CTest suite. The command builder test asserts that all
three module paths remain registered.

The separate Repository Test Inventory remains the read-only source browser.
This preserves a clear distinction between discovery and execution while
making the three tests part of the normal GUI test action.

## Design principles and industry practices

The design borrows a small set of established practices:

1. **Native test authority:** CMake/CTest owns native test registration;
   unittest owns Python test execution; the regression harness owns case
   comparison. The GUI orchestrates rather than reimplements them.
2. **Discovery is not execution:** source scanning must never synthesize a
   pass/fail result. Results require a real process execution and exit status.
3. **Explicit execution registry:** commands with special setup are allow-listed
   and tested. Do not infer safe execution from a filename or AST declaration.
4. **Layered verification:** fast unit tests, subsystem integration tests,
   production-path tests, and representative CFD regressions establish
   different scopes of confidence.
5. **Traceability and claim scoping:** contract metadata links principle,
   implementation, documentation, and evidence, while each evidence record
   states its scope and limitations.
6. **Reproducibility:** execution records should ultimately retain revision,
   configuration, platform/backend, command, exit status, and artifact paths.
7. **Headless core:** discovery, catalog validation, and command construction
   remain independent of Qt so they can be tested and reused in CI.

## Validation commands

From the repository root:

```powershell
python -m unittest tools.verification.test_curve_catalog
python -m unittest tools.verification.test_contract_catalog
python -m unittest tools.oneflow_gui.test_test_inventory
python -m unittest tools.oneflow_gui.test_test_runner
ctest --test-dir build -C Debug --output-on-failure
```

The three focused Contract Workbench commands validate catalog and inventory
behavior. The command-builder test validates GUI registration. CTest validates
the configured native suite. CFD regression remains a separate validation
step.
