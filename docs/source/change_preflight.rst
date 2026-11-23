OneFLOW Change Preflight
========================

The Change Preflight is a read-only repository consistency check. It is
intended to run before a functional-slice commit and in Python Tooling CI.

Run from the repository root:

.. code-block:: console

   python -m tools.verification.preflight
   python -m tools.verification.preflight --changed

The optional --changed switch lists staged, unstaged, and untracked paths.
It does not narrow validation to those paths: the full contract catalog and
Python test registry are checked deliberately, because a changed implementation
can affect a contract or test declaration in another directory.

Checks performed
----------------

* Load all supported Numerical, Geometry, and versioned domain contracts.
* Reuse the Contract Catalog path and schema validation.
* Resolve supported evidence references, including Python unittest class/method
  declarations and linked GoogleTest declarations.
* Reject duplicate evidence identities within one contract.
* Discover Python unittest source files and compare them with the existing
  _PYTHON_UNIT_TESTS execution registry and documented exclusions.
* Require each resolved python_test contract evidence source to be covered by
  that execution registry.
* When a Python evidence source is missing, search repository Python files for
  the exact declared unittest class/method identity and report matching
  candidate paths.
* When a linked GoogleTest source is missing, search repository C++ test files
  for the exact GoogleTest macro identity and report matching candidate paths.
  Comments and string literals are masked before matching; build/cache
  directories are excluded. A unique match is suggested; multiple matches are
  listed for manual disambiguation.
* For a resolved GoogleTest source with a statically recognized CMake target,
  require a matching ``gtest_discover_tests`` declaration. Simple local
  CMake function/macro wrappers with literal arguments are expanded, including
  the repository's ``add_oneflow_gtest`` helper.
The command exits with status 0 when checks pass, 1 when consistency issues
are found, and 2 when the check cannot run (for example, Git changed-file
discovery fails).

Configured CTest evidence
-------------------------

After configuring and building a local tree, run CTest at least once, then
inspect the populated CTest registry without running tests:

.. code-block:: console

   python -m tools.verification.ctest_evidence --build-dir build --configuration Debug

The command uses ``ctest -N --show-only=json-v1``. The repository uses
``gtest_discover_tests(... DISCOVERY_MODE PRE_TEST)`` for its GoogleTest
executables, so individual cases are materialized when CTest first runs. The
native Linux and Windows CI workflows run this report after CTest execution.
It compares exact GoogleTest
contract identities with configured CTest test names and checks target
provenance using the discovered executable name or CMake backtrace file. It
reports ``discovered``, ``not_discovered``, ``target_unmapped``,
``discovery_not_declared``, ``discovery_deferred``, and ``stale_source`` separately. Optional build
features may legitimately be absent from a particular configured tree; use
``--require-all`` only when that build configuration is expected to contain
every declared GoogleTest.

A ``discovered`` result proves only that this configured build tree currently
lists the matching test. It does not prove that the test binary has been run
or that the test passed. The GUI Test Workbench writes CTest JUnit output into
the build directory and captures one durable ExecutionRecord per reported
testcase. Contract IDs are attached only after exact test-name and CMake target
provenance matching; other CTest cases remain unlinked records. Actual
execution outcomes belong in durable ExecutionRecords and CI test results.

Interpretation and boundaries
-----------------------------

A PASS means that the supported declarations and repository references are
internally consistent. It does not mean that C++ targets compile, CTest or
unittest cases execute successfully, CFD results are physically valid, or
contract evidence has passed. Runtime verification remains the responsibility
of the existing build and test workflows.

GoogleTest evidence is checked in distinct stages. The Contract Catalog checks
the source declaration. The preflight then checks CMake target mapping and
CTest discovery when the source-to-target relationship is expressible through
the supported static reader. A mapped target without ``gtest_discover_tests``
is a consistency failure. An unmapped source remains explicitly unresolved:
the reader does not evaluate arbitrary CMake variables, conditions, generated
lists, or general function semantics. This static check does not prove that
CMake configured, the target compiled, CTest discovered the test in a build
tree, or the test executed. The native CMake/CTest jobs remain the authorities
for those runtime stages.
When a new Python test module is added, either register it in
tools/oneflow_gui/test_runner.py or add a justified entry to the existing
_PYTHON_TEST_EXCLUSIONS map. Contract evidence should point to an existing
test declaration and its source file; the preflight then checks that the
module is included in the execution registry. If the path is stale but the same exact unittest or GoogleTest identity exists
elsewhere, the report suggests the candidate path(s). C++ matching ignores
comments and string literals to avoid treating examples as declarations.
Review the match and update the Contract explicitly; the preflight never
rewrites evidence or assumes that a matching identity proves the intended
test is the same after a refactor.

Recommended functional-slice workflow
-------------------------------------

#. Implement the related code and tests.
#. Add or update contract implementation, documentation, and evidence links.
#. Register executable Python tests (or document a deliberate exclusion).
#. Run the focused tests and python -m tools.verification.preflight.
#. Run the applicable CMake/CTest and Python test workflows.
#. Commit the implementation, tests, inventory/CI wiring, and documentation
   together.

The tool intentionally reuses existing catalogs and registries. Do not add a
second hand-maintained list of contract files, Python test modules, or
evidence counts here.


Contract execution status
-------------------------

The Contract Workbench computes execution status separately from the declared
Contract lifecycle status. Only explicit Contract-ID links and exact CTest test
identities (or exact declared Python test source modules) can contribute.
Current-commit results are distinguished from stale or unknown revisions.
Failed and blocked outcomes are not converted into passes, and Debug/Release
and platform contexts remain separate. A context with incomplete evidence is
reported as partially verified or not run. The evaluator is read-only and does
not infer an expected CI matrix or modify Contract declarations.

Validate the status evaluator with:

.. code-block:: console

   python -m unittest tools.verification.test_contract_verification
