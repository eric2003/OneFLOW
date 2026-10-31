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
   include ``time_integral``, ``rk_stage``, ``reynolds``, ``mach_ref``, ``vismodel``,
   ``inflow_intensity``, ``ref_sa``, ``turb_cfl_ratio``, ``nmg``, and ``mgCycleType``.

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

