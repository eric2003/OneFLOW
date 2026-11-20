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

The command exits with status 0 when checks pass, 1 when consistency issues
are found, and 2 when the check cannot run (for example, Git changed-file
discovery fails).

Interpretation and boundaries
-----------------------------

A PASS means that the supported declarations and repository references are
internally consistent. It does not mean that C++ targets compile, CTest or
unittest cases execute successfully, CFD results are physically valid, or
contract evidence has passed. Runtime verification remains the responsibility
of the existing build and test workflows.

GoogleTest evidence is checked for a source declaration by the Contract
Catalog. This first preflight does not claim that every GoogleTest is included
in a CMake target or discovered by CTest. The static CMake test catalog and
native CMake/CTest jobs remain the authorities for those separate stages.

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
