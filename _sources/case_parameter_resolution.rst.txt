Case parameter resolution and construction
=============================================

Purpose
-------

Parameter files are a serialization format consumed by the OneFLOW runtime.
Their physical split (``cfd.txt``, ``inflow.txt``, ``turb.txt``, and so on) is
not the identity or ownership model of a parameter. The Workbench should
therefore resolve requested changes by logical key and retain the physical
declaration path only as provenance.

A concrete example is ``maxSteps`` in
``test/cfd/plateuns2dslau2/script/cfd.txt``. The 50-step value is useful for a
fast regression test that checks residual behavior at every iteration. A
production-length run may instead require 35,000 steps, as demonstrated by
the related ``plateuns2dslau2_34950_35000`` case. Changing the number should
not require the user to know that the current example stores it in
``cfd.txt``. Another Case may place the same logical key in a different
manifest-listed file, or may not have a file named ``cfd.txt`` at all.

Resolution boundary
-------------------

``parameter_resolution.py`` introduces a Qt-independent resolution layer:

* callers provide loaded ``ParameterDocument`` objects and logical overrides;
* ``ParameterModel`` resolves the effective declaration using manifest order
  and the existing last-declaration-wins rule;
* the resolver produces a side-effect-free plan containing old value, new
  value, declared type, override origin, and actual source path;
* the plan can be previewed before it is applied;
* application renders every target file before writing, stages replacements
  beside their source files, and checks for external edits before replacement;
* if a later replacement fails, earlier replacements are rolled back from
  their original bytes before the plan reports failure.

The first slice deliberately updates existing declarations only. It rejects
missing keys and duplicate requests instead of silently inventing a file,
inserting a declaration with guessed syntax, or choosing an arbitrary target.
The writer retains comments, spacing, encoding, line endings, and unrelated
declarations. Multi-file replacement is not a single filesystem transaction:
the plan stages all files and performs best-effort rollback if an individual
replace fails. The GUI must report an incomplete rollback distinctly if the
filesystem also rejects restoration. A later schema-backed operation can add a new parameter only
when its type, serialization rule, applicability, manifest placement, and
runtime semantics are known.

This is not a solver feature and does not change the C++ runtime. It is a GUI
and test-infrastructure capability: it allows different front ends to express
the same logical change and makes the resulting configuration change auditable.

Multiple Case construction paths
--------------------------------

The same resolver is intended to serve several authoring paths:

#. **Derive from a validated Case** — copy a known profile and apply a small
   set of explicit differences, such as changing a regression profile from
   50 to 35,000 outer iterations.
#. **Select a versioned Case profile** — choose compatible task, dimension,
   equation system, physical model, grid, boundary conditions, and numerical
   methods, then resolve profile values through the same logical-key layer.
#. **Use a GUI parameter form** — edit values grouped by physical meaning,
   without exposing the legacy file layout as the user-facing organization.
#. **Import and normalize an existing Case** — load its ordered manifest,
   classify declarations, show duplicate and unsupported entries, and apply
   reviewed overrides.
#. **Inject workflow or experiment overrides** — generate reproducible
   parameter variants for parameter sweeps, regression matrices, or
   verification workflows.

These paths must converge on one resolution and validation pipeline rather
than implement independent file-editing logic.

Validation and evidence
-----------------------

A resolved parameter plan is not yet a valid simulation setup. The intended
validation stages are distinct:

* **Resolution** — every requested key maps to one effective declaration;
  conflicts and missing keys are reported.
* **Serialization** — only intended value tokens change; encoding, comments,
  manifest order, and unrelated files remain stable.
* **Production ingestion** — the generated Case is read through the same
  production manifest/configuration loader used by OneFLOW.
* **Initialization** — the selected task, solver, physical models, mesh, and
  boundary setup complete their real initialization path.
* **Numerical verification** — residual histories, conservation, reference
  solutions, and regression outputs are checked at an appropriate run length.

The 50-step flat-plate test and the 35,000-step production configuration are
different run profiles, not interchangeable evidence. A successful parameter
write or parser test must never be described as numerical validation. Runtime
acceleration and performance optimization are separate concerns; they should
not be mixed into the correctness of parameter resolution.

Initial tests
-------------

``test_parameter_resolution.py`` covers a ``maxSteps`` override from 50 to
35,000 when the declaration lives in an arbitrarily named
``run-control.txt``; case-insensitive logical lookup; duplicate declarations
with effective-source routing; missing-key rejection; duplicate override
rejection; preview without mutation; and external-edit detection before
application, successful updates across two files, and rollback when the second
file replacement is forced to fail. These tests exercise the GUI-side resolution boundary, not the
production C++ parser or a complete CFD run.

Versioned Case profiles
-----------------------

The first profile layer is implemented in ``case_profile.py``. A profile is a
versioned declaration of intent, source Case, purpose, logical-key overrides,
maturity, and limitations. It is deliberately independent of the source
parameter filenames. ``resolve_case_profile`` delegates to the common
parameter resolver; it does not introduce a second file-editing mechanism.

Two initial flat-plate profiles are provided:

* ``laminar-flat-plate-regression-50.json`` requests ``maxSteps = 50``
  for short regression checks.
* ``laminar-flat-plate-production-35000.json`` requests
  ``maxSteps = 35000`` for a production-length run configuration.

Both are marked ``proposed``. The production profile changes only the
iteration limit; it does not claim that a fresh 35,000-step simulation has
been run, converged, or numerically validated. In particular, the existing
``plateuns2dslau2_34950_35000`` case is a separate continuation/restart
scenario and must not be treated as evidence for a fresh production run.

Profiles are starting points, not rigid solver schemas. Future profile
composition must account for task, equation system, model compatibility,
grid and boundary requirements, restart state, output policy, and explicit
overrides. Incompatible combinations should be diagnosed rather than silently
merged. Each authoring path (template, GUI, import, workflow, batch study)
must converge on the same resolver and production validation pipeline.
