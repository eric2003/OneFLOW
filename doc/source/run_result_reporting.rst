Run-result status and evidence contract
========================================

Purpose
-------

OneFLOW tools must make outcomes obvious without requiring users to inspect
large JSON documents or infer success from output-file existence. The terminal
is the human-first summary. A versioned JSON report is the machine-readable
record consumed by GUIs, automation, tests, and later workflow stages.

Outcome layers
--------------

Do not collapse these independent claims into one generic validated flag:

* **Execution** — the process completed and its exit status is known.
* **Structure** — the artifact is readable and satisfies implemented
  structural invariants such as counts, connectivity, topology, and groups.
* **Mesh quality** — numerical quality metrics have been evaluated against
  stated acceptance criteria. Calculating or displaying a metric alone is
  not a quality pass.
* **CFD validation** — the mesh has been consumed by the intended OneFLOW
  runtime and passed declared solver/workflow checks. A smoke run is not a
  physical accuracy benchmark unless the case and reference establish that.

A successful mesh-generation process can therefore have successful execution
and structural checks while mesh quality remains NOT_CHECKED and CFD
validation remains NOT_RUN. The STL-to-Gmsh report explicitly records the
separate CFD validation stage as NOT_RUN; that does not invalidate a successful
mesh-generation result because CFD validation is outside this command.

Stage states
------------

Every stage uses exactly one of these case-sensitive values:

.. list-table::
   :header-rows: 1

   * - State
     - Meaning
   * - PASS
     - The stage ran and every declared assertion passed.
   * - FAIL
     - The stage ran and at least one declared assertion failed, or it raised
       an error.
   * - NOT_RUN
     - The stage was not reached because an earlier stage failed or execution
       was blocked.
   * - NOT_CHECKED
     - The capability or criterion is not implemented or was deliberately not
       evaluated in this run.

Never use PASS merely because a file exists, a subprocess returned zero, or
a report was written. Do not conflate NOT_RUN with NOT_CHECKED.

Versioned report envelope
-------------------------

The stable envelope is oneflow.run-report version 1:

.. code-block:: json

   {
     "schema": "oneflow.run-report",
     "version": 1,
     "workflow": "STL -> Gmsh",
     "status": "success",
     "exit_code": 0,
     "created_at_utc": "2026-11-16T00:00:00Z",
     "stages": {
       "input": {"status": "PASS"},
       "geometry": {"status": "PASS"},
       "gmsh": {"status": "PASS"},
       "mesh_structure": {"status": "PASS"},
       "physical_groups": {"status": "PASS"},
       "mesh_quality": {"status": "NOT_CHECKED"},
       "cfd_validation": {"status": "NOT_RUN"}
     },
     "metrics": {},
     "artifacts": {},
     "error": null
   }

The status is success only when the process-level run succeeds; otherwise it
is failed. The exit_code is the command result (the current STL-to-Gmsh
adapter uses 0 for success and 2 for handled failures). The error is null on
success and otherwise contains an error type and human-readable message.
Workflow-specific legacy fields may be retained during migration, but must
not redefine shared status semantics.

Human-readable terminal output
------------------------------

A tool should finish with a short summary containing:

* workflow name and an unmistakable SUCCESS or FAILED verdict;
* each stage and its state;
* a small set of useful counts and metrics;
* output artifact paths and report path;
* process exit code;
* the first actionable error on failure.

Do not dump the entire JSON report to the terminal by default. Keep detailed
diagnostics in the report and, when useful, a separate log. External-tool
progress output should not bury the final verdict.

Failure and stale-report behavior
---------------------------------

* Write a failure report when the destination is known, including the failed
  stage, the error, and later stages marked NOT_RUN.
* Remove or invalidate an earlier report before beginning a new run. A stale
  success report must never describe a failed new attempt.
* Write reports atomically in the destination directory so readers never see
  partially written JSON. Generated files and report temporaries belong in
  the designated external workspace, never in the checked-in source tree.
* The process exit code remains authoritative for shell/CI control flow. The
  report supplements it with stage detail; stdout text must not be parsed as
  the machine interface.

Current STL-to-Gmsh implementation
----------------------------------

grid/stl_to_gmsh/generate_volume_mesh.py is the first adapter implementing
this contract. Its stages are input, geometry, Gmsh execution, mesh structure,
and Physical Groups. Mesh quality is NOT_CHECKED: current surface triangle
diagnostics detect degeneracy and calculate descriptive metrics, but the
adapter does not yet apply a complete tetrahedron-quality acceptance policy.
CFD validation is a separate workflow and is not implied by this report.

The adapter retains its legacy validated field for compatibility. It means
only that the implemented geometry, MSH structure, and Physical Group checks
passed. It does not claim tetrahedron-quality or CFD validation.

Future Gmsh-to-CGNS, CGNS-to-OneFLOW, and Mesh-to-CFD adapters should preserve
the envelope, stage vocabulary, exit-code semantics, and distinction between
execution, structural, quality, and CFD evidence. Add domain-specific stages
and metrics rather than inventing incompatible status flags.


Gmsh-to-CGNS conversion report
------------------------------

The C++ Gmsh-to-CGNS adapter writes three complementary files in the Case
results/ directory:

* <stem>.gmsh-read-report.json — detailed reader diagnostics;
* <stem>.gmsh-report.md — human-readable mesh inspection and boundary mapping,
  headed by an explicit Gmsh -> CGNS: SUCCESS verdict;
* <stem>.run-report.json — the shared oneflow.run-report envelope.

The conversion report marks MSH input, Gmsh reading and CGNS writing PASS only
after the writer has completed and the output file exists. CGNS read-back is
NOT_CHECKED by the command (the independent C++ conversion tests reopen and
inspect the file); mesh quality is NOT_CHECKED until numerical acceptance
criteria are implemented; CFD validation is NOT_RUN until the separate
Mesh-to-CFD workflow is run. A successful conversion is therefore not a
claim of solver ingestion or physical accuracy.

Mesh-to-CFD Validation report
-----------------------------

The workflow report uses the same top-level envelope and stage vocabulary.
Its stages are mesh_generation, gmsh_to_cgns, cgns_to_oneflow, and cfd_smoke.
A blocked or unvisited stage is NOT_RUN. Mesh quality and CFD physical
accuracy remain NOT_CHECKED unless their own acceptance checks are
implemented and executed. The detailed per-Case messages remain in steps
and each stage log remains beside its staged Case outputs.

The GUI displays the concise terminal summary, including the process exit
status and each stage. JSON is the machine interface; the process exit code
remains authoritative for pass/fail control flow.


In the Mesh-to-CFD dialog, each completed run report is read after its process
exits and shown as a compact per-Case tree with stage, status and diagnostic
message. The detailed live log is retained separately; the report tree is not
reconstructed by guessing from log text.
