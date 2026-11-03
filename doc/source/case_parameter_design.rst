Case Parameter Design and GUI Editing
=====================================

Status and scope
----------------

This page records the current boundary between legacy case files, the C++ runtime
parameter database, and the optional Workbench GUI. The checked-in RAE2822 cases are
examples for understanding the format; their values are not recommendations for an
optimal mesh or CFD solution. A value appearing in a database dump proves that it was
parsed and committed, but does not prove that it affects the selected task.

Configuration path
------------------

``script/control.txt`` is an ordered manifest of parameter files. For example,
``test/grid/rae2822convert/script/control.txt`` lists ``simu.txt``, ``dimension.txt``,
``tolerence.txt``, ``grid.txt``, and ``part.txt``. The CFD case
``test/cfd/rae2822_roe_sa/script/control.txt`` additionally separates solver, model,
inflow, wall-distance, and multigrid settings into files such as ``cfd.txt``,
``vismodel.txt``, ``inflow.txt``, ``walldist.txt``, ``turb.txt``, and ``multigrid.txt``.

``ControlManifestReader.cpp`` reads the ordered manifest and resolves each
listed script path under the selected case directory. ``ControlFileReader.cpp``
then orchestrates the listed files in that order. ``ConfigLoader`` parses each
legacy file into a ``ConfigDocument``. ``ControlScriptParser.cpp`` invokes
``ConfigDatabaseAdapter::Commit()``, whose implementation is in
``ConfigDatabaseAdapter.cpp``. The adapter applies those entries through
``ProcessData`` to the legacy runtime database. The adapter preserves document order;
repeated names follow the runtime ``DataPara`` replacement policy. The current commit
operation is incremental, not transactional: an error in a later entry does not roll
back earlier entries. This is a useful migration boundary:
the parser handles text syntax, while task and physics code remain responsible for
runtime meaning.

The runtime database is not an authoring model. ``DataEntry`` stores a parameter name,
runtime type, array size, and value. ``DataPara`` indexes entries with an
``unordered_map``. Consequently, a database dump does not preserve the originating
file, declaration order, comments, whether a value was explicit or defaulted, units,
help text, constraints, or the conditions under which a parameter applies. A GUI must
not use the dump as a complete parameter editor schema or treat its order as stable.

Sample parameter domains
------------------------

The following groups describe where the example values are declared. They are a useful
starting point for navigation, not an exhaustive list of runtime consumers.

Grid conversion case
   ``test/grid/rae2822convert/script/simu.txt`` selects the Grid task and conversion
   mode. ``grid.txt`` declares source and target formats and paths, scale and
   translation, boundary-name maps, and format-specific switches. ``part.txt`` holds
   partition-related settings. ``dimension.txt`` and ``tolerence.txt`` provide shared
   dimensional and tolerance settings. Code in ``GridTypes.cpp``, ``Plot3D.cpp``, and
   ``DomainInp.cpp`` reads representative grid type, transform, format, source, and
   boundary parameters; ``SimuDef.cpp`` resolves the task and ``Tolerance.cpp`` reads
   ``ptTol``. Representative declared keys include ``gridObj``, ``sourceGridType``,
   ``targetGridType``, ``gridScale``, ``gridTrans``, ``sourceGridFileName``,
   ``targetGridFileName``, and ``parallelStrategy``.

RAE2822 Roe-SA CFD case
   ``simu.txt`` holds run control, files, and output cadence. ``cfd.txt`` selects
   equation/solver controls; ``vismodel.txt`` selects the viscosity model;
   ``inflow.txt`` contains inlet and reference-state settings; ``turb.txt`` contains
   turbulence controls; ``multigrid.txt`` configures the multigrid cycle; and
   ``walldist.txt`` configures wall-distance input. Representative consumers include
   ``NsCom.cpp`` and ``NsCtrl.cpp`` for flow controls, ``TurbCom.cpp`` and
   ``TurbCtrl.cpp`` for turbulence, ``Multigrid.cpp`` for multigrid controls, and
   ``MultiBlock.cpp`` for grid and wall-distance paths. Representative declared keys
   include ``time_integral``, ``rk_stage``, ``timestepModel``,
   ``global_dt``, ``cflst``, ``cfled``, ``ncfl``, ``max_time_ratio``,
   ``reynolds``, ``mach_ref``, ``vismodel``, ``inflow_intensity``,
   ``ref_sa``, ``turb_cfl_ratio``, ``nmg``, and ``mgCycleType``.

These source references confirm that representative parameters reach implementation
code. They do not classify every key in the two database dumps. Some settings are
conditional on task, model, or file format; some may be defaults, compatibility fields,
derived values, or historical leftovers. Establishing whether a parameter is unused
requires tracing its read path and conditions, not searching for one exact spelling
alone.

GUI parameter editor direction
------------------------------

The Workbench currently manages projects, case files, workspaces, runs, workflows, and
artifacts. Its workflow artifact bindings can update declared settings in staged
parameter files, but this is a narrow execution adapter, not a general-purpose editor
for every C++ runtime parameter.

A future parameter editor should be driven by a versioned, declarative schema that
records, per parameter:

* stable key and legacy spelling;
* value type and scalar/array shape;
* owning parameter file and display group;
* description, units, allowed choices, and validation rules when known;
* applicable task, solver, physics model, and format conditions;
* whether the field is required, optional, advanced, deprecated, or derived;
* evidence/source for its default and runtime consumer.

The editor should preserve unknown keys and avoid deleting comments or user ordering
when writing legacy files. It should distinguish values explicitly present in a file
from values inserted by defaults or derived during setup. Conditional visibility must
follow declared applicability rules rather than guesses from parameter names. A schema
entry with unknown units, range, or consumer should say so instead of inventing a
constraint.

Configuration adapter tests
----------------------------

``TextFileParser_test`` now also exercises the production
``ReadScriptFileNameList(caseDir, ...)`` operation with a temporary manifest.
It checks that blank lines are skipped, listed paths are rooted under the case,
and declaration order is preserved. This pins manifest discovery, not parsing
or committing the listed scripts.

``control_parameter_resolver_test`` also exercises the production
``ConfigDatabaseAdapter`` against the real global ``DataBase``. It checks typed
integer, real, and string entries and confirms that a later declaration of the same
name replaces the earlier value. These are adapter/database contract tests; they do
not exercise manifest ordering across multiple physical files, parallel distribution,
or the complete ``ReadControlInfo(caseDir)`` operation.

Current parser tests
--------------------

``tests/io/ConfigLoader_test.cpp`` specifies scalar type/name/value parsing, declaration
order, quoted scalar values, literal arrays spanning lines, callback-based named array
lengths, rejection when a named length has no resolver, and replacement of the document
when the loader parses another file. These tests exercise syntax translation only. They
do not assert runtime parameter validity, units, defaults, solver behavior, or editor
serialization.
First structured-editor foundation
-----------------------------------

The first implementation slice adds a Qt-independent legacy parameter document adapter
and source-backed metadata for a small set of scalar declarations.

* simu.txt: simutask selects Grid or Solve.
* cfd.txt: maxSteps, maxTime, time_integral, and ischeme describe iteration/time controls
  and documented numerical-method choices.
* vismodel.txt: vismodel and visname describe the viscous-model family and named choices
  documented by checked-in case files.

parameter_document.py detects UTF-8, UTF-8 BOM, and GBK; parses scalar declaration
value spans; and updates only changed value tokens. It retains comments (including GBK
Chinese comments), spacing, blank lines, declaration order, and line endings. It also
rejects a save if the file changed externally after it was loaded. ParameterSpec metadata
is deliberately limited to keys whose meaning and choices are supported by checked-in
configuration comments and representative cases.

This is the data and serialization foundation for the structured editor UI. It does not
yet edit arrays, nested/conditional fields, or every declaration in the control manifest.
Unknown declarations remain untouched. The current CFD examples place task selection in
simu.txt, iteration controls in cfd.txt, and viscous-model selection in vismodel.txt.
No generic parameter is assigned an invented unit, range, or solver meaning.

The parser tests cover typed values, GBK and UTF-8 BOM preservation, CRLF and comment
preservation, unknown syntax preservation, no-op formatting stability, and external-edit
conflict detection.

Logical parameter model
-----------------------

The physical files listed by ``control.txt`` are an input organization, not the
parameter model itself. The same logical CFD parameter set can be assembled from
``simu.txt``, ``dimension.txt``, ``tolerence.txt``, ``cfd.txt``, ``vismodel.txt``,
``inflow.txt``, ``turb.txt``, and other files. Grid cases follow the same pattern.
The Workbench should therefore present a unified parameter view and retain the
originating file for each value, rather than treating each file as an independent
parameter database.

The Qt-independent ``parameter_model.py`` is the first implementation of this
boundary. It combines multiple ``ParameterDocument`` instances in the order supplied
by the case's control manifest. It identifies ``simutask``, ``dimension``, and
``ptTol`` as shared parameters; recognizes initial CFD groups (run control, output,
numerical methods, reference conditions, physical models, turbulence, multigrid, and
wall distance); and recognizes initial Grid groups (grid definition, formats/files,
and partitioning). Classification is deliberately conservative: keys not yet in the
registry remain visible under ``Other`` or ``Unclassified`` rather than being dropped
or assigned invented metadata. The declared type and value always come from the source
declaration, not from a guessed schema type.

The model preserves declaration order and source paths. When the same key occurs in
multiple files, the last declaration is the effective value, matching the runtime
``DataPara`` replacement behavior; earlier declarations remain available as shadowed
entries for diagnostics. The model does not modify source files. ``control.txt`` stays
outside the parameter collection because it is the ordered manifest that tells the
loader which parameter files to read.

This is a logical consolidation only. It does not merge or rewrite the legacy files.
That distinction allows the Workbench to offer one coherent editor even if users later
choose to store a subset of parameters in one file or split them across several files.
Physical consolidation should be considered only after runtime loading order,
duplicate-key behavior, compatibility with existing cases, and round-trip preservation
have explicit tests. The current document parser still handles scalar declarations;
array declarations such as ``gridTrans[3]`` remain unsupported by the scalar model.
The Project Window conservatively detects declaration-like array statements and
other declaration types outside the scalar parser and lists them with source line
and reason. This scan is not exhaustive and must not be presented as complete
parameter coverage.

SQLite parameter catalog
------------------------

The Workbench stores a local database at ``%USERPROFILE%/.oneflow/parameter-catalog.sqlite3`` (or ``ONEFLOW_HOME/parameter-catalog.sqlite3``). This is an application index, not a replacement for legacy case files or the C++ runtime ``DataBase``.

The catalog stores project roots, relative file paths, source encoding, SHA-256 fingerprints, parsed typed declarations, schema-backed parameter definitions, and immutable content snapshots. Re-indexing an unchanged file is a no-op. A changed file adds a snapshot and refreshes the current-value view. Opening a project indexes supported parameter files without modifying them. Comments, formatting, line endings, and unknown declarations remain in the source file; the catalog stores parsed values and source identity.

The index is intentionally limited to scalar declarations understood by ``ParameterDocument``. Arrays and conditional declarations require explicit parser contracts before they can be represented safely. SQLite uses the Python standard library and adds no third-party runtime dependency.

Recommended migration sequence
------------------------------

#. Keep the existing parser and runtime database adapter stable while collecting
   source-backed parameter metadata for one task at a time.
#. Keep parser-contract tests for scalar metadata, declaration order, multiline arrays,
   named-length resolution, and repeated-parse replacement. Add schema-specific validation
   after a schema exists. Test serialization round trips at the file-writing boundary; do
   not encode numerical recommendations as parser rules.
#. Build a read-only categorized view in Workbench before enabling edits. Compare it
   with the originating files and report unsupported or unclassified entries.
#. Enable writing only for fields with a known source file, type, serialization rule,
   and validation contract. Keep the original file content recoverable.
#. Audit parameter consumers incrementally and mark evidence as confirmed, conditional,
   or unresolved. Do not infer that a field is unused from a successful run.

This sequence makes the catalog useful even before a full editor exists and lets the
cost of isolated schema and serialization tests expose coupling in the current runtime
configuration boundary.

Related documentation
---------------------

* :doc:`runtime_architecture` describes runtime configuration and task setup.
* :doc:`workbench_architecture` describes Workbench domain and service boundaries.
* :doc:`workbench_user_guide_zh` and ``tools/oneflow_gui/README.md`` describe current
  user-visible GUI behavior.



Runtime consumer evidence: first audited slice
--------------------------------------------

The Qt-independent parameter_consumer_catalog.py records a small set of source-level reads checked against the current hexin source. It is an evidence index, not a complete static-analysis engine. Each record names the parameter key, repository-relative source path, consuming function, read expression, read timing, applicability, and a short behavior/evidence note.

The initial audited keys are:

* maxSteps — required read in Iteration::Init. It initializes the iteration limit; SimuIterState::Running applies it when iexitflag == 0.
* time_integral — required reads in Ctrl::Init and NsCtrl::Init. These are distinct initialization paths, not duplicate declarations.
* mach_ref — required read in NsCom::Init and used with gama_ref to derive const_cp.
* vismodel — required read in VisModel::Init. The selected model subsequently gates viscous/turbulence validation and wall-distance behavior.
* ptTol — lazy read through Tolerance::Init on the first Tolerance::GetTol call; it is not necessarily read during every startup path.
* nmg — required read in MG::Init; values at or below one normalize the pre/post smoothing counts.

The catalog distinguishes required, optional, and lazy read kinds. A key absent from this curated list is not classified as unused: it is simply not audited in this slice. Likewise, a direct database read proves a consumer exists, not that the parameter affects every run or that its numerical recommendation has been verified. Conditional applicability and downstream numerical effects must be traced separately.

The focused unit tests verify case-insensitive lookup, multiple read sites, lazy-read classification, behavior notes, and the unknown-key boundary. The next extension should add consumer records only after checking the corresponding hexin implementation and should link them to task/model applicability and Formula & Knowledge evidence without inferring formulas from key names.

Time-integration parameter trace
---------------------------------

The consumer index also traces the pair ``time_integral`` and ``rk_stage``
beyond their declarations. Four control initializers read ``time_integral``
directly: ``Ctrl::Init``, ``NsCtrl::Init``, ``INsCtrl::Init``, and
``TurbCtrl::Init``. The corresponding four initializers read ``rk_stage``
and construct ``rk_coef[s] = 1 / (s + 1)``.

``TimeIntegral::Init`` maps the selector values declared in
``TimeIntegral.h`` as follows: ``MULTI_STAGE`` (1) selects ``RungeKutta``,
``SIMPLE`` (3) selects ``Simple``, and the remaining values take the ``Lusgs``
branch. The selector's ``LUSGS`` value is 2; the current dispatch code does not
show an explicit invalid-value rejection in this function. This is a description
of the implementation, not a recommendation to accept arbitrary values.

On the finest grid, ``TimeIntegral::RungeKutta`` iterates over ``rk_coef`` and
assigns each coefficient to ``ctrl.lhscoef`` before the residual/update task
sequence. The context-aware path may use the resident Euler backend only when
``EvaluateEulerRungeKuttaCapability`` accepts the solver, grid, equation count,
inviscid scheme, explicit time-integral selection, physical terms, limiter,
interface-exchange, and backend capability gates. Otherwise it falls back to
the legacy task sequence. This accelerated-path gate is narrower than the
general ``time_integral`` selector.

The Formula & Knowledge entry ``time_integral_rk_stage_weight`` links the
parameter, control initialization, stage loop, and implementation documentation.
It intentionally has no theory curve: the current source establishes the stored
stage-weight schedule, but does not by itself establish a complete Butcher
tableau, formal temporal order, stability region, convergence rate, or CFD
accuracy. Those require separate numerical analysis and case-level verification.

The contract schema also distinguishes the stable navigation ID from the
legacy runtime key. For example, ``parameter/max-steps`` is the Workbench
contract ID while ``parameter_key`` is ``maxSteps``, the exact key read by
``GetDataValue<int>``. This explicit mapping avoids guessing camel-case runtime
identifiers from human-readable slugs.

The Parameter contract view exposes these records in a dedicated **Runtime
Consumers** tab. It displays audited read sites, read timing, applicability,
behavior notes, downstream source/symbol links, and related Formula & Knowledge
IDs. Consumer source paths are also added to the implementation file navigator.
The source audit checks recorded paths and symbols; it is not a runtime execution
record. A parameter without an entry is marked not yet audited, never unused.

The consumer catalog's source-audit test checks that every recorded direct read,
consumer symbol, and named downstream behavior symbol still resolves in the
current source tree. It does not claim to discover every consumer in the repository.

Maximum iteration steps
-----------------------

The parameter contract ``parameter/max-steps`` documents ``maxSteps`` as an
outer-iteration limit, not as a convergence tolerance or a physical time
limit.

``Iteration::Init`` reads ``maxSteps`` as an integer and copies it to both
``Iteration::maxSteps`` and ``Iteration::maxIterSteps``. The actual stop
condition is selected in ``SimuIterState::Running``:

* when ``ctrl.iexitflag == 0``, the loop stops once
  ``Iteration::outerSteps >= Iteration::maxSteps``;
* when ``ctrl.iexitflag != 0``, the loop instead stops when
  ``ctrl.currTime >= ctrl.maxTime``.

Consequently, changing ``maxSteps`` does not establish a time limit, a
residual-convergence target, or a bound on every inner/sub-iteration. The
current source path also does not establish an admissible input range or an
explicit positivity check. Those policies should be documented only after
their actual parser, setup, and runtime validation paths are traced.

The Contract Workbench places this parameter in the **Run Control** group.
Its evidence is source-declared and marked ``declared_not_executed``; it is
not a recorded test result or a claim that a case converges within the limit.

Maximum simulation time
-----------------------

The related parameter contract ``parameter/max-time`` records the time-based
outer-loop stop condition. ``Ctrl::Init`` reads ``maxTime`` as a ``Real`` and
initializes ``currTime`` to zero. In ``SimuIterState::Running``, the comparison
``currTime >= maxTime`` is active when ``iexitflag`` is nonzero. When
``iexitflag`` is zero, the loop instead uses the ``maxSteps`` iteration limit.

The audited code establishes this comparison but does not, by itself, define
the physical unit or nondimensional scale of ``maxTime`` for every solver and
case. Nor does it establish how each solver advances ``currTime`` or guarantee
residual convergence by the time threshold. The parameter should therefore be
interpreted together with the case's time convention and solver configuration.

Both ``maxSteps`` and ``maxTime`` appear in the **Run Control** group. They are
alternative outer-loop termination controls selected by ``iexitflag``; neither
is a residual-based convergence criterion.

Evolving parameter classification policy
----------------------------------------

The logical parameter model classifies parameters along independent axes. These
axes are deliberately separate so that the current task selector is not mistaken
for evidence of a parameter's ownership.

* **Domain / ownership** identifies the subsystem that owns the setting:
  ``Shared``, ``CFD``, ``Grid``, or ``Other``.
* **Functional group** describes the behavior being controlled. The initial
  groups include ``General``, ``Run Control``, ``Time stepping``, ``Numerical
  methods``, ``Flow and reference conditions``, ``Physical models``,
  ``Turbulence``, ``Multigrid``, ``Wall distance``, ``Output``, ``Grid
  definition``, ``Grid formats and files``, and ``Partitioning``.
* **Task applicability** records whether a setting is known to apply to
  ``Both``, ``Solve``, ``Grid``, or remains ``Unknown``. This is not inferred
  from the currently selected task.
* **Source provenance** remains attached to each declaration: source file,
  manifest/declaration order, and any shadowed declaration are preserved.

A known key is assigned only when the current catalog has semantic evidence
for that assignment. An unknown key is retained under ``Other /
Unclassified`` with applicability ``Unknown``; it is not promoted to CFD merely
because ``simutask`` currently selects ``Solve``. Source-file location may help
guide a later audit, but a filename alone is not treated as proof of ownership.

The first classification refinement separates **Time stepping** from
**Numerical methods**. ``global_dt``, ``timestepModel``, ``idualtime``,
``ncfl``, ``cflst``, ``cfled``, and the currently catalogued time-step model
selectors belong to the time-stepping group. The time-integration selector
(``time_integral``), Runge-Kutta stage count (``rk_stage``), and spatial
discretization selectors remain under Numerical methods. This is a navigation
and semantic grouping decision; it does not assert that CFL values, global
time-step settings, and local/global time-step algorithms are interchangeable.
Their formulas, units, defaults, validation ranges, and runtime effects must
continue to be audited independently.

Classification is intentionally versionable. As source consumers and real case
files are audited, a key may move between groups or domains, and a new group may
be introduced when it represents a distinct physical or numerical concept.
Such changes should include a source-backed rationale and focused model tests.
Unrecognized declarations and unsupported syntax remain visible throughout
that process; completeness is increased incrementally rather than simulated by
assigning every key to a known category.
