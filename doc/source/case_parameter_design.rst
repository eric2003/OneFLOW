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

Representative case inventory
-----------------------------

The manifest-driven classification test now loads all parameter files listed by
the following checked-in cases, rather than validating only one CFD case and
one grid case:

* ``test/cfd/plateuns2dslau2``: laminar flat plate with SLAU2;
* ``test/cfd/turbplateuns2droe_sa``: turbulent flat plate with Roe-SA;
* ``test/cfd/rae2822_roe_sa``: transonic airfoil with Roe-SA;
* ``test/cfd/m6wingroe_sa``: three-dimensional wing with Roe-SA;
* ``test/grid/rae2822convert``: grid conversion and partitioning.

The inventory test fails if any scalar or unsupported declaration name in
these manifests is assigned to ``Other / Unclassified``. It also pins selected
functional distinctions: turbulence inflow inputs are separated from model
controls; physical inflow and wall thermal conditions remain separate from
turbulence inflow; output cadence is separate from result/restart filenames;
solver sub-iteration limits remain run-control parameters; and all ``mg*``
settings, including turbulence-viscosity switches and factors, are grouped
under ``Multigrid`` rather than split across physics groups.

This is a growing representative corpus, not a claim that all ``test/`` cases,
all optional solver modes, or all legacy configuration files have been
enumerated. New case families should be added when their manifests introduce
new keys or materially different applicability. Existing classifications may
be refined as source consumers and case comments provide stronger evidence.

These source references confirm that representative parameters reach implementation
code. They do not classify every key in every repository case manifest. Some settings are
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
* maxTime — required read in Ctrl::Init. It initializes the simulation-time threshold and currTime; SimuIterState::Running applies it when iexitflag != 0, as the alternative to maxSteps.
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
  groups include ``General``, ``Run Control``, ``Solver system``, ``Time
  stepping``, ``Numerical methods``, ``Flow and reference conditions``,
  ``Thermodynamic reference state``, ``Transport properties``, ``Boundary and
  inflow conditions``, ``Physical models``, ``Turbulence inflow conditions``,
  ``Turbulence``, ``Multigrid``,
  ``Wall distance``, ``Output control``, ``Output files``, ``Input and output
  files``, ``Grid definition``, ``Grid formats and files``, and ``Partitioning``.
* **Task applicability** records whether a setting is known to apply to
  ``Both``, ``Solve``, ``Grid``, or remains ``Unknown``. This is not inferred
  from the currently selected task. ``Solve`` identifies the top-level solver
  task family; it does not claim that every solver model or physical regime
  consumes the parameter.
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

The next catalog expansion was checked against the repository's representative
case manifests and declarations: ``test/cfd/rae2822_roe_sa/script/control.txt``
and ``test/grid/rae2822convert/script/control.txt``. This exposed several
distinct concepts that should not be hidden inside broad buckets:

* Solver equation counts and equation-model selectors are grouped as
  ``Solver system``.
* Iteration/output intervals and display/dump switches are ``Output control``;
  result and restart filenames are ``Output files``.
* Flow input/output roots and the grid reference are ``Input and output files``.
* Inflow selection and wall thermal boundary settings are ``Boundary and inflow
  conditions``.
* Gas-state setup and dimensional reference quantities are ``Thermodynamic
  reference state``; Prandtl and Schmidt numbers are ``Transport properties``.

These groups are still an initial classification, not a declaration that every
case-specific option has been inventoried. The manifest is the authoritative
starting point for discovering files; actual declarations and C++ consumers
remain the evidence for adding keys. Arrays such as ``gridTrans[3]`` continue
to be reported as unsupported by the scalar editor even when their names can
be recognized for navigation.

CFL ramp and time-step strategy: source trace
---------------------------------------------

The RAE2822 solver case declares ``cflst``, ``cfled``, ``ncfl``,
``timestepModel``, and ``global_dt`` together, but these keys control
different layers of the time-step process.

* ``Iteration::Init`` reads ``cflst`` (ramp start), ``cfled`` (ramp end),
  and ``ncfl`` (ramp length). It also reads ``idualtime`` when the key is
  present; legacy files without it retain ``dualtime = 0``.
* ``TimeStep::CalcCfl`` selects an iteration index through
  ``CflIterationIndex``. Physical-time mode uses ``outerSteps``; dual-time
  mode uses ``innerSteps``. ``CflAtIteration`` linearly interpolates while
  ``iteration < ncfl`` and returns ``cfled`` afterwards:

  .. math::

     \mathrm{CFL}(i) = \mathrm{CFL}_{start} +
       \frac{\mathrm{CFL}_{end} - \mathrm{CFL}_{start}}{n_{CFL}} i

  The source contract is in ``TimeStepPolicy.h/.cpp``. It defines the ramp
  schedule; it does not establish a universally optimal CFL value or a
  convergence guarantee for a particular mesh and flow.
* ``NsCom::Init`` reads ``timestepModel``. The unstructured solver dispatches
  0 to local time stepping, 1 to global time stepping, and the remaining
  branch to the local/global hybrid path (the representative case comment
  documents 2 as ``lg``).
* Solver control initialization reads ``global_dt`` into ``Ctrl::pdt``.
  ``UTimeStep::CalcGlobalTimeStep`` applies this value in global mode. The
  hybrid path computes local steps, takes their minimum, and writes that
  value back to ``Ctrl::pdt``. Thus ``global_dt`` is not interchangeable
  with the CFL ramp and is not the controlling value in every strategy.

``currTime`` is runtime state, not a declared input parameter. The solver
control initializers set it to zero; ``SimuIterState::Running`` compares it
with ``maxTime`` when ``iexitflag`` selects time-based termination. When
``iexitflag`` is zero, the outer iteration count is instead compared with
``maxSteps``. This separates the time-step calculation from the outer-loop
termination policy.

Output cadence is another independent group. ``Iteration::Init`` loads
``nFieldSave``, ``nVisualSave``, ``nResSave``, and ``nForceSave``. The current
``Iteration::ResOk`` and ``Iteration::ForceOk`` implementations show that
residual and force output cadence can be indexed by inner steps for
``startStrategy == 2`` and by outer steps otherwise. Representative NS,
incompressible NS, turbulence, and multigrid paths use ``nFieldSave`` to
schedule restart output; NS and incompressible NS paths use ``nVisualSave``
to schedule visualization tasks. This is a source-audited slice, not proof
that every solver and every output backend has identical scheduling semantics.

Turbulence inflow inputs have a separate functional group from general
turbulence-model controls. ``TurbCom::Init`` reads ``inflow_intensity`` and
``inflow_viscosity``. Their audited consumer is the
``transition_model == ITReGama`` branch inside ``TurbCom::InitInflow``,
not every two-equation turbulence initialization path. In that branch,
``inflow_viscosity`` supplies ``mutoo``; the initializer first forms
``k_e = mutoo * kwoo / (cmu * Re * rho_ref)``. For positive
``inflow_intensity``, it replaces ``k_e`` with

.. math::

   k_e = \frac{3}{2}\,(0.01 I)^2

and recomputes ``kwoo = k_e * (cmu * Re * rho_ref) / mutoo``.
For non-positive intensity it keeps the initial viscosity-based ``k_e``
and obtains the empirical intensity through ``trans.CalcIntensity``.
The source's ``0.01`` scaling is recorded literally; this does not establish
a universal external unit convention, valid range, or recommended value.
The non-transition two-equation initialization follows a separate
``ref_sst``/default path.

The time-step strategy and configured global step are also exposed as independent Parameter contracts: parameter/time-step-strategy (timestepModel) and parameter/global-time-step (global_dt). The six inputs are therefore represented by four searchable Parameter entries: the grouped parameter/cfl-ramp schedule (cflst, cfled, ncfl), plus idualtime, timestepModel, and global_dt as separate controls.

The optional dual-time selector is also exposed as its own Parameter contract:
parameter/dual-time-index (runtime key idualtime). The existing parameter/cfl-ramp
contract remains the combined algorithm entry for cflst, cfled, and ncfl, while
timestepModel and global_dt remain separately source-audited within that time-step
strategy trace. This avoids treating a selector for the CFL iteration counter as
equivalent to the time-step strategy or the physical time-step input.

Iteration::Init checks for idualtime before reading it. If the declaration is
absent, Iteration::dualtime remains zero; when it is one, the CFL schedule uses
innerSteps rather than outerSteps. The consumer catalog records this as an optional
read, not a required parameter. This does not prove every solver path implements
dual-time stepping or establish a temporal accuracy order.

The consumer catalog now links both keys to the direct typed reads and the
``TurbCom::InitInflow`` implementation. This is code-path evidence, not a
claim that every turbulence model uses the same inflow closure.
The Parameter runtime-consumer catalog now links the three CFL ramp inputs
to ``Iteration::Init``, ``TimeStep::CalcCfl``, ``CflIterationIndex``, and
``CflAtIteration``. It separately links ``timestepModel`` to the local/global
dispatch, ``global_dt`` to its four solver-control read sites and the
unstructured time-step strategy, and the output intervals to their iteration
and representative solver scheduling consumers. These links are source evidence, not
execution records. Units, admissible ranges, default values, and behavior
outside the audited solver path remain unresolved until separately verified.


Turbulence CFL ratio: unsteady spectrum scaling
------------------------------------------------

The parameter ``turb_cfl_ratio`` is read as a ``Real`` in
``TurbCom::Init``. Its downstream use is narrower than a global CFL
controller: ``UTurbSpectrum::CalcUnsteadySpectrum`` multiplies the
selected per-cell time-step by this ratio before assembling the implicit
unsteady contribution to the turbulence equations.

For the single-time-step branch (``idualtime == 0``), the code forms

.. math::

   t_s = \Delta t_{turb} \, r_{turb}, \qquad
   S_{unsteady} = \frac{sp_2}{t_s} V,

where ``r_turb`` is ``turb_cfl_ratio``, ``V`` is cell volume,
and ``sp_2`` is the time-integration coefficient. In dual-time mode,
the turbulence time-step is scaled in the same way, while the physical-time
term ``sp_2 / pdt1`` is added separately. Therefore this parameter is not
the same control as the mean-flow CFL ramp (``cflst/cfled/ncfl``) or
the solver time-step strategy (``timestepModel``).

The current source audit establishes the read and the two assembly branches.
It does not establish a universal recommended value, admissible range,
default, or numerical stability guarantee; those remain unverified until
case-level sensitivity and convergence studies are recorded.


Spalart-Allmaras production limiter and viscosity-ratio audit
-------------------------------------------------------------

The ``rprod`` input is read by ``TurbCom::Init`` and passed along the
Spalart-Allmaras source-term path to ``SaLimitedDiffusion``. The policy
implementation in ``TurbSaPolicy.cpp`` is

.. math::

   D_{lim} = \min(r_{prod} P, D),

where ``P`` and ``D`` are the production and diffusion contributions passed
by the SA source assembly. This documents the implementation relation only.
The current audit does not establish a recommended value, admissible range,
or applicability to SST and other turbulence models.

The ``max_vis_ratio`` input is read by ``TurbCom::Init`` and passed
through ``TurbSrcFlux`` to both the SA and SST eddy-viscosity calculations.
The shared ``TurbCappedEddyViscosity`` policy returns

.. math::

   \mu_t = \min(\mu_{t,\mathrm{uncapped}}, r_{\mu,\max} \mu_l),

where ``r_{\mu,\max}`` is ``max_vis_ratio`` and ``\mu_l`` is
the local laminar viscosity. Pure-policy tests cover applying the cap and
preserving an uncapped value below it. These tests do not establish a
recommended value, admissible range, or end-to-end CFD accuracy.


Parameter worked examples and evidence scope
--------------------------------------------

Each parameter contract now carries at least one structured worked example
with an explicit input, calculation, result, and limitations. The examples
are explanatory source-derived arithmetic or pure-policy checks, not
completed CFD case executions. In particular:

* ``turb_cfl_ratio`` demonstrates how the scaled turbulence time step
  changes the implicit unsteady spectrum coefficient;
* ``rprod`` follows the minimum-selection rule in
  ``SaLimitedDiffusion`` and its pure-policy unit test;
* ``max_vis_ratio`` follows the shared SA/SST eddy-viscosity cap and its
  pure-policy tests;
* ``inflow_intensity`` and ``inflow_viscosity`` are restricted to the
  audited ``ITReGama`` transition-inflow branch.

A worked example is explanatory evidence only. A linked test declaration,
a hand calculation, a parsed parameter value, and a full solver run are
different evidence types and must not be promoted into one another. Parameter
defaults, valid ranges, dimensional conventions, and case-level sensitivity
remain open wherever the contract explicitly says they are not yet verified.

Parameter catalog completeness and Case-template readiness
-----------------------------------------------------------

The Parameter entry in Contract Workbench and the Project Window parameter
tree are related views, but they have different responsibilities:

* The Project Window tree shows declarations found in the opened Case's ordered
  ``control.txt`` manifest. It preserves source files, declaration order,
  effective values, and overridden declarations. It is case-local.
* The Contract Workbench Parameter catalog is a repository-level knowledge
  index. Each contract has a stable ID and an explicit ``parameter_key``
  mapping to the legacy runtime key. It can describe source consumers,
  formulas, aliases, examples, and evidence independently of one Case's
  current values.
* The runtime database remains the authority for the values consumed by a
  particular run. A catalog definition is not itself a runtime assignment.

At the current inventory checkpoint, the Parameter contract directory
contains 150 entries covering 152 individual parameter keys. The original
14 source-audited contracts cover 16 keys (the CFL ramp contract groups
``cflst``, ``cfled``, and ``ncfl``). A further 136 entries were added
from a declaration inventory of five representative CFD Cases and one Grid
Case. These new entries are deliberately marked ``proposed``: they establish
names and functional groups in Contract Workbench, but do not claim a verified
runtime consumer, formula, unit, default, admissible range, or numerical
validation.

The current directory contains 150 JSON contract definitions. The 14
source-audited entries retain their existing source mappings and knowledge;
the 136 proposed entries are inventory placeholders intended to make the
full set navigable and auditable. Proposed entries intentionally have no
fabricated implementation path or verification result. In particular,
``gridTrans`` is represented in the catalog but its array declaration is
still reported by the scalar Project Window parser as unsupported; catalog
presence does not imply that the editor can safely rewrite it.

The inventory groups include run control, solver system, time stepping,
output control and files, numerical methods, flow/reference conditions,
thermodynamic state, transport properties, boundary/inflow conditions,
physical models, turbulence, multigrid, wall distance, grid definition,
grid formats/files, and grid partitioning. Shared task/dimension and geometry
tolerance entries are also included. The exact runtime key spelling is
preserved as ``parameter_key``.

This is a meaningful increase in Parameter Workbench coverage, not proof of
repository-wide completeness. The first inventory pass is based on the
representative Case set; it does not yet enumerate every Case manifest or every
``GetDataValue``/optional/lazy consumer in production code. A parameter
entry being present in the catalog does not guarantee that its consumer has
been audited or that its formula, units, default, and valid range are known.
Unknown declarations and array declarations must remain visible rather than
being silently discarded.

The next inventory pass must enumerate declarations from every repository
Case manifest and cross-check them against production parameter reads. It
must also distinguish:

* a key declared in a Case file;
* a key directly read by production code;
* a key whose downstream behavior has been traced;
* a key with a documented physical or numerical definition;
* a key with a tested admissible range or default;
* a key validated by a representative executable Case.

These are independent facts. The inventory should preserve the exact spelling
of the runtime key, declaration type, source file, manifest order, aliases,
classification, and evidence state. A missing consumer audit means
``not_yet_audited``, not unused. A key with no known formula remains a
parameter entry with an explicit explanation, not an invented equation.

Case parameter-file generation is a later consumer of this catalog, not a
shortcut around it. A generated Case should be assembled from a compatible
profile of task, dimension, solver/equation system, physical models, grid,
boundary conditions, and numerical methods. The generator must:

#. start from a known Case template or an explicitly versioned parameter
   profile, preserving required legacy file syntax and the ordered manifest;
#. include only parameters applicable to the selected task and model, while
   retaining required defaults and dependencies;
#. report unresolved required values, unsupported arrays, unknown keys, and
   conflicting declarations before writing files;
#. preserve existing file encoding and formatting where editing an existing
   Case; never rewrite source-tree examples in place;
#. validate the generated manifest and parameter documents with the production
   parser, then run a real setup/initialization smoke test;
#. record the generated profile version, source Case, parameter provenance,
   validation result, and any user modifications.

A generated file is therefore a candidate Case configuration until the
production reader accepts it. Successful parsing is not proof of a stable
simulation, numerical accuracy, or physical validity. The first usable
profiles should be derived from existing, understood Cases (for example,
laminar flat plate and RAE2822 Roe-SA), with small, explicit differences
recorded against their source profiles. Arbitrary combinations of individually
known parameters must not be treated as a valid solver configuration.
