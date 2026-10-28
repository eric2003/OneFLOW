Test Inventory and Contract Data
=================================

Purpose
-------

The Workbench should describe tests that actually exist in the repository, not
maintain a second manually curated list that can drift away from source code.
The first step is read-only discovery in tools/oneflow_gui/test_inventory.py:

* C++ tests are discovered from GoogleTest TEST, TEST_F and TEST_P declarations
  under tests.
* Python unittest methods are discovered from test_*.py files under tools.
* The Workbench inventory filters frameworks independently; GoogleTest and
  unittest are not presented as one undifferentiated test type.
* Each record includes a stable ID, framework, suite, name, source location,
  discovery method, execution hint, and separate input, expected and result
  fields.
* Scenario JSON is linked by explicit framework/suite/name identity. Only the
  linked contract test receives that scenario's input and expected tables.
  An integration test can link to the scenario without incorrectly inheriting
  the helper's input and expected data.

The inventory is a source index, not an execution engine. In the Workbench it
belongs in the repository/developer testing area, separate from a scenario's
execution details; it must not appear as a peer of the per-run Result view.
Frameworks are separate groups: GoogleTest (C++), Python unittest, and any
future framework each have their own section/filter. Do not merge GoogleTest
cases and Python test methods into one undifferentiated list.

Listing a test does not mean it has run. The result field stays empty until an
execution backend records an actual run. A scenario fixture backend may compare
fixed fixture output, but it must not populate a C++ test's result. A suite-level
summary such as "878/878 passed" is useful run evidence, but it is not a
per-test result record unless the runner can associate that result with the
specific test identity.

Contract model
--------------

The intended record shape is:

.. code-block:: json

   {
     "id": "gtest:FaceSearchBasicTest.Example",
     "framework": "GoogleTest",
     "suite": "FaceSearchBasicTest",
     "name": "Example",
     "source": "tests/geometry/FaceSearchBasicTest.cpp",
     "input": null,
     "expected": null,
     "result": null
   }

Input and expected are optional test-contract data, not fields to fill with
guesses. Many existing unit tests express behavior through C++ setup and
assertions rather than serializable tables. For those tests, source location
and the assertion remain authoritative until a meaningful data adapter exists.

Scenario bootstrap from GoogleTest
-----------------------------------

``tools/geometry_workbench/scenario_import.py`` provides the first write
path from the discovered GoogleTest inventory into a scenario JSON file.
It preserves the suite, test name, source path, and source line, and creates
editable ``input`` and ``expected`` objects with ``result: null``.
It refuses to overwrite an existing scenario unless the caller explicitly
opts in. This is a provenance-preserving scaffold, not a claim that arbitrary
C++ setup and assertions have been translated into data.

The dialog also exposes **Save scenario**. It persists the edited input and
expected JSON tables back to the selected scenario while retaining test identity
and provenance. Saving a contract does not write a run result: ``result`` remains
the existing value, and only an execution backend may record new execution
evidence. This separates authoring a hypothesis (expected data) from observing
the implementation (actual result).

The Geometry Contracts dialog now exposes **Create from GoogleTest...**.
It lists discovered GoogleTest declarations and writes a non-overwriting scenario
skeleton into the scenario directory. This is intentionally a user-triggered
import rather than an automatic bulk conversion: a generated skeleton is a
starting point for a reviewed contract, not an executable model of the C++
test. The source file and line are retained so a developer can inspect the
original setup and assertions while filling the data tables.

The distinction matters for cases such as
``CollectBoundaryInterfaceChildNodeMappingsTest.RejectsMissingParentRowInsteadOfUsingSharedChildGlobalRow``.
Its C++ test constructs ``LinkField``, ``InterfaceFaceSplitTable``, calls
the production function, and checks both an exception and unchanged output.
The current ``parent_local_mapping_required.json`` is a manually authored
contract derived from that test; source discovery does not generate its
specific input arrays or expected exception. The next increment should extract
contracts from an explicit source annotation or a small set of supported
assertion/setup patterns, rather than guessing that all C++ expressions are
serializable.

Incremental path
----------------

#. **Discover** existing tests from source and link explicit scenarios. This is
   the initial capability; it does not require a separate test database.
#. **Standardize** scenario-backed tests on input, expected and result. Preserve
   the distinction between fixture comparison and native C++ execution.
#. **Execute** selected tests through their real runners. CTest target layout
   must be inspected before promising per-GoogleTest-case execution; a CTest
   target is not necessarily the same thing as an individual GoogleTest case.
#. **Record** command, configuration, start/end time, exit status and actual
   output for each execution. Results should be generated by the runner, never
   manually typed into the expected-value field.
#. **Add adapters selectively** where a production algorithm has naturally
   serializable input/output, such as registered child-face resolution and
   parent-local mapping publication.

This keeps the source as the authority for test existence and uses scenario
files only for explicit, reusable data contracts. It avoids creating a large
hand-maintained database before the execution and data models are understood.
