OneFLOW Change Preflight and Change Gate
==========================================

The repository provides two complementary checks before a functional-slice
commit. The Change Preflight validates repository declarations and references.
The Change Gate adds change-aware Python test selection and execution. Neither
replaces the complete cross-platform CI workflows.

Change Preflight
----------------

Run from the repository root:

.. code-block:: console

   python -m tools.verification.preflight
   python -m tools.verification.preflight --changed

The read-only preflight reuses the Contract Catalog and Python test inventory.
It checks contract schema and paths, supported evidence declarations,
duplicate evidence identities, Python unittest declarations and registration,
and statically recognized GoogleTest/CMake discovery relationships.
`--changed` lists staged, unstaged, and untracked paths; it does not narrow
repository-wide consistency checks.

Change Gate
-----------

Run the affected Python tests for the current working tree:

.. code-block:: console

   python -m tools.verification.change_gate --changed

Preview the commands without running them:

.. code-block:: console

   python -m tools.verification.change_gate --changed --dry-run

For CI, compare the checked-out commit with the event base revision:

.. code-block:: console

   python -m tools.verification.change_gate --base <base-commit-sha>

The Python Tooling workflow fetches the required Git history and invokes this
mode on pushes and pull requests. Manual workflow dispatch runs `--all`.

Run the complete registered Python suite:

.. code-block:: console

   python -m tools.verification.change_gate --all

The gate always runs the repository-wide Change Preflight first. The Python
Tooling CI invokes the gate after Python compilation and before the detailed
service test steps; push and pull-request runs use the event base revision so
the clean CI checkout still has a meaningful changed-file set. If that fails,
test execution is stopped. It then selects commands from the existing
`_PYTHON_UNIT_TESTS` registry in `tools/oneflow_gui/test_runner.py`; it does
not maintain a second list of Python test modules.

Selection policy
----------------

* STL-to-Gmsh changes select the registered STL-to-Gmsh and STL-import tests,
  plus the directly related mesh validation, geometry asset, and mesh
  generation Workbench tests.
* STL-import changes select the registered STL-import and STL-to-Gmsh tests,
  plus the related cylinder-cut Workbench tests.
* Mesh-layout changes select registered mesh-layout tests and related
  Workbench tests.
* `tools/oneflow_gui/` changes select registered GUI/service tests.
* `tools/verification/` changes select registered verification tests.
* Documentation-only changes do not select Python tests.
* An unrecognized code area (including native solver, CMake, or other
  implementation paths) conservatively selects the complete registered
  Python suite. Mixed known and unknown code changes also select the full
  suite. The gate must never silently treat an unmapped code change as having
  no affected tests.

The selector only executes commands already present in the central registry.
Registry exclusions remain exclusions and are not silently converted into
tests; the full C++ build/CTest and platform CI workflows remain required for
native changes. A zero-test result is appropriate only when the changed paths
are documentation-only or the change set is empty.

Exit codes are 0 for a passing preflight and selected tests, 1 for consistency
or test failures, and 2 when Git change discovery or gate setup cannot run.
A dry run reports planned commands but is not test evidence.

Boundaries
----------

A passing Change Gate means the declared consistency checks and selected
Python commands passed in the current environment. It does not prove that
native targets compile, CTest passes, Gmsh/CGNS integrations execute, or CFD
results are physically valid. Mesh-quality diagnostics remain diagnostics,
not implicit blockers. Continue to run the applicable Windows/Linux build,
CTest, Python Tooling, and external-tool integration workflows.

Configured CTest evidence
-------------------------

After configuring and building a local tree, run CTest at least once, then
inspect the populated CTest registry without running tests:

.. code-block:: console

   python -m tools.verification.ctest_evidence --build-dir build --configuration Debug

The command uses ``ctest -N --show-only=json-v1``. Discovery proves only that
the configured build lists a test; it does not prove execution or success.
Actual outcomes belong in durable ExecutionRecords and CI results.

Recommended functional-slice workflow
-------------------------------------

#. Implement related code and tests.
#. Update the relevant contract, documentation, and evidence links.
#. Register Python tests in the existing execution registry (or document a
   deliberate exclusion).
#. Run the Change Gate and review its selected command list.
#. Run applicable CMake/CTest and full cross-platform CI workflows.
#. Commit implementation, tests, registry/CI wiring, and documentation
   together.

The tools deliberately reuse existing catalogs and registries. Do not add
parallel hand-maintained contract lists, Python test inventories, or evidence
counts.
