Runtime Case Lifecycle and Production Entry
==============================================

Purpose
-------

This page describes the production path from the public ``Simulation::Run()``
entry point to case configuration, task dispatch, case teardown, and process
finalization. It is an architecture and verification contract, not a claim that
every solver or accelerator path has been validated.

The distinction between process and case ownership is central:

* Process resources (parallel runtime and accelerator backend) span all cases
  in one invocation.
* Case resources (control database values, project paths, solver registrations,
  fields, and case diagnostics) belong to one case.
* A case teardown must not finalize process resources when another case may
  follow.

Production sequence
-------------------

The full-simulation command line is parsed once. The process context is
initialized once, and each case gets a separate ``SimuImp`` / ``SimuContext``.
The sequence is:

.. code-block:: text

   Simulation::Run()
       |
       +-- RunImpl()
             |
             +-- ParseSimulationCommandLineOptions()
             +-- Prj::Init()
             +-- processContext.SetupProcessEnvironment()
             |      +-- SetUpParallelEnvironment()
             |      +-- InitializeAccelRuntime()
             |
             +-- for each case directory
             |      |
             |      +-- SimuImp(caseDir)
             |      +-- SimuImp::RunCase()
             |             |
             |             +-- PreProcess() -> SetupCaseEnvironment()
             |             |      +-- bind case root and log
             |             |      +-- ReadControlInfo(caseDir)
             |             |      +-- read typed CaseIOSettings
             |             |      +-- ApplyCaseIOSettings()
             |             |      +-- publish environment-ready state
             |             |
             |             +-- MainProcess() -> ResolveTaskFromControl()
             |             |                  -> TaskRegistry::Create()
             |             |                  -> task.Execute(context)
             |             |
             |             +-- PostProcess() -> TeardownCase()
             |                    +-- release case-owned backend state
             |                    +-- reset solver / field / zone state
             |                    +-- clear configuration and project paths
             |                    +-- clear context task state
             |
             +-- processContext.FinalizeEnvironment()
                    +-- release remaining context-owned backend state
                    +-- FinalizeAccelRuntime()
                    +-- HXFinalize()

Configuration ingestion contract
---------------------------------

``ReadControlInfo(caseDir)`` remains the production reader. It discovers
the real parameter documents under ``script/``, reads them in deterministic
filename order, commits their typed values to the global ``DataBase``, and
emits ``log/database.log``. A filename-only ``control.txt`` is checked
against the real files for consistency; it does not select files or define
parameter precedence. If ``control.txt`` contains parameter declarations,
it is loaded as an ordinary parameter document. ``SetupCaseEnvironment()`` then obtains the case I/O
settings and applies project, input, output, mesh, and flow roots through the
existing ``Prj`` compatibility API. The ``io_paths.txt`` artifact records the
resolved paths.

A test that substitutes a parser or directly populates ``DataPara`` does not
prove this production ingestion path. The production integration tests must
use real temporary files and the real reader.

Lifecycle invariants
--------------------

1. Readiness is published only after configuration loading and path application
   complete successfully.
2. A failed case setup is followed by the real ``RunCase()`` cleanup path;
   the original setup exception remains the reported failure.
3. Case teardown releases case state but leaves the process runtime active.
4. Process finalization releases backend-dependent state before finalizing the
   accelerator backend.
5. If backend finalization fails, process lifecycle ownership remains active so
   finalization can be retried.

Verification map
----------------

The runtime Contract Catalog entry is
``runtime/process-lifecycle``. Its evidence links the following production
integration boundaries:

* ``FailedProductionSetupRunsRealCaseTeardownAndKeepsProcessRuntimeActive``:
  failure after a partial control-file read, followed by case teardown.
* ``SuccessfulCasesRemainIsolatedWithinOneProcess``: two sequential cases,
  distinct typed values, per-case artifacts, and process-runtime continuity.
* ``FullSimulationEntryReadsControlFilesAndFinalizesProcessRuntime``: the
  public ``Simulation::Run()`` entry point, multi-file control ingestion, task
  dispatch, generated artifacts, and process finalization.
* ``FailedFullSimulationEntryTearsDownCaseAndFinalizesProcessRuntime``: a
  missing second control script after a partial read, with rollback verified
  through the public entry point and process finalization.
* Accelerator lifecycle tests: dependency-ordered release and retry after an
  injected backend finalization failure.

Evidence status is deliberately explicit. A declaration of a linked GoogleTest
is not an execution record. The Contract Workbench must continue to display
``linked_not_executed`` until a test run produces a result and evidence artifact.

Repository case setup-only evidence
-------------------------------------

``RepositoryRae2822ConversionCaseCanBeInitializedWithoutRunningGridTask``
copies the checked-in ``test/grid/rae2822convert/script`` directory to a
temporary case directory and invokes the production
``SimuContext::SetupCaseEnvironment()``. It checks typed values
(``dimension == 2``, ``simutask == "Grid"``, and
``ptTol == 1.0e-10``) and the materialized
``log/database.log`` and ``results/io_paths.txt`` artifacts. The test
then invokes the real ``TeardownCase()``, verifies that case readiness is
cleared while process runtime remains active, and only then finalizes the
process runtime. A scoped guard retries case teardown during early test exit.

The test deliberately does not resolve or execute the Grid task. It therefore
proves that this repository control-script set can be ingested and its case
environment can be torn down; it does not prove that the grid conversion
inputs are mutually consistent or that conversion succeeds.

Known limits and next checks
----------------------------

The current tests use a serial runtime and a registered lightweight task. They
do not establish multi-rank MPI behavior, physical CUDA/HIP driver behavior,
all solver-specific cleanup paths, or numerical correctness. The repository
case test is linked evidence until an exact-revision CTest run records its
execution. R1 remains Partial until the focused test, full suite, and
Windows/MSVC and Linux CI have passed on the same revision. Numerical
regression remains a separate evidence layer.

Architecture implication
------------------------

The lifecycle is best understood as two nested ownership scopes rather than a
single global cleanup routine:

.. code-block:: text

   Process scope:  [ setup -------------------------------- finalize ]
   Case scope:          [ setup -- execute -- teardown ]
   Case scope:                                  [ setup -- execute -- teardown ]

This separation is the foundation for safe multi-case execution and for
eventually presenting a single traceable workbench view:

``Contract -> Test -> Run -> Evidence -> Result``.

Each claim should identify its production boundary, its test, the latest
execution result, and the remaining unverified conditions.
