OneFLOW Runtime Environment and Dependency Information
========================================================

Purpose
-------

OneFLOW exposes runtime and build-environment information through the
OneFLOW executable. The first implementation focuses on CGNS because the
actual CGNS library loaded by a process is important evidence when diagnosing
mesh conversion and file I/O behavior.

The Workbench does not link to CGNS. It asks the selected OneFLOW executable
for a versioned JSON document and displays the information returned by the
runtime. This keeps the GUI independent of third-party library headers and
ABI details.

Running OneFLOW
---------------

The examples below assume a Windows Debug build from the repository root.
Adjust the executable path and configuration for Release or an installed
build.

Informational commands do not require a Case directory:

.. code-block:: powershell

   & .\build\bin\Debug\OneFLOW.exe --version
   & .\build\bin\Debug\OneFLOW.exe --environment
   & .\build\bin\Debug\OneFLOW.exe --environment --format json

Equivalent Linux commands (from the repository root):

.. code-block:: console

   ./build/bin/OneFLOW --version
   ./build/bin/OneFLOW --environment
   ./build/bin/OneFLOW --environment --format json

The version command prints the OneFLOW application version. The environment
command prints a human-readable report. The JSON form writes a machine-readable
document to standard output. Its schema_version field versions the JSON
contract independently of the application version.

To run a normal CFD Case, pass the r command and the Case directory:

.. code-block:: powershell

   & .\build\bin\Debug\OneFLOW.exe r .\test\cfd\plateuns2dslau2

Use the actual Case directory and executable selected for the run. Informational
commands are handled before the simulation lifecycle starts; they do not
initialize a Case, MPI process environment, solver, or mesh.

Saving and consuming the JSON report
------------------------------------

The JSON command writes to standard output and does not create a file by
itself. To save a point-in-time report, redirect standard output:

.. code-block:: powershell

   & .\build\bin\Debug\OneFLOW.exe --environment --format json |
       Set-Content -Encoding utf8 .\environment-report.json

On Linux:

.. code-block:: console

   ./build/bin/OneFLOW --environment --format json > environment-report.json

A Python caller can invoke OneFLOW directly and parse the response:

.. code-block:: python

   import json
   import subprocess

   executable = "build/bin/Debug/OneFLOW.exe"
   completed = subprocess.run(
       [executable, "--environment", "--format", "json"],
       check=True,
       capture_output=True,
       text=True,
   )
   environment = json.loads(completed.stdout)

The Workbench uses this same process boundary. It invokes the configured
OneFLOW executable with the JSON command and parses standard output. The GUI
does not link to CGNS and does not infer a runtime version from an installation
directory.

CGNS version evidence
---------------------

The report deliberately separates three kinds of information:

* runtime_version is populated only when the configured CGNS headers and
  linked library expose the optional cg_get_cgns_version and
  cg_get_cgns_version_string APIs. CMake checks declaration and linkability
  so OneFLOW can compile the corresponding call. This is a build-capability
  check; it is not itself evidence of the version loaded at runtime.
* header_version is derived from the CGNS_VERSION macro in the header used to
  compile OneFLOW. It is compile-time evidence, not proof of the version of a
  different library that the operating system may load.
* runtime_library.module_path and runtime_library.linkage identify the module
  containing the actual cg_open function address. On Windows this identifies
  the loaded DLL or the executable when CGNS is linked into it. On Linux and
  macOS it uses dynamic-loader module information. A static build may
  therefore report the OneFLOW executable as the containing module.
* version_compatibility summarizes the evidence without confusing inference
  with an independent runtime query. If both versions are available, OneFLOW
  reports verified_match or mismatch. If only the header version is available,
  it reports assumed_match and explicitly records header_version as the basis.
  If neither side provides enough evidence, it reports unknown.

The public CGNS 4.5.2 header does not declare the two optional version-query
functions. For that version the runtime API field remains unavailable, while
the header version and loaded-module information remain visible. The
compatibility result is therefore assumed_match, not a runtime verification.
The default assumption is a reporting convention, not a guarantee of ABI
compatibility.

Do not confuse the library version with cg_version(fn, &version). The latter
queries the CGNS file-format version stored in an opened file; it does not
identify the version of the CGNS library used by the process.

Workbench usage
---------------

Start the GUI from the repository root:

.. code-block:: powershell

   python -m pip install -r tools/requirements-gui.txt
   python -m tools.oneflow_gui

Open Runtime Environment... on the home screen. The dialog queries the
currently selected executable when opened/refreshed. If the executable is
changed in the Configuration dialog, the report describes that selected
executable, not necessarily the executable used by an earlier Case run.

The table columns mean:

* Dependency: dependency name.
* Runtime version: version returned by a runtime API, including its source.
  Unavailable means no reliable runtime version was returned.
* Header version: version macro from the header used at compile time.
* Version consistency: verified_match when independently obtained versions
  agree, mismatch when they differ, assumed_match when only the header version
  is available, or unknown when evidence is insufficient. The basis and
  diagnostic detail are included in the JSON report.
* Status: whether the dependency is enabled/available and, where applicable,
  the runtime-version query status.
* Loaded module: actual module path containing the queried library function,
  plus linkage classification when known.

Status available together with Runtime version Unavailable means the CGNS
library is loaded and usable, but this CGNS build does not expose the optional
library-version query functions. It is not a load failure.

Current persistence behavior
----------------------------

The standalone JSON command writes to standard output; it does not create a file
by itself. The Project Workbench and Case Workbench queue capture the selected
executable's environment at run start and store it under
``~/.oneflow/workflow-data/executions/<run-id>/`` (or the equivalent
``ONEFLOW_HOME`` root). The snapshot is written even if the query fails: a
failed query produces a small JSON diagnostic document and does not prevent
the CFD run from starting. A user can still save a standalone point-in-time
report using the redirection examples above.

Workbench settings, workflow definitions, run history, and execution
workspaces have separate storage responsibilities under the OneFLOW user data
root (``~/.oneflow`` by default; ``ONEFLOW_HOME`` can relocate that root).
Workflow definitions and history are stored under
``~/.oneflow/workflow-data/definitions`` and
``~/.oneflow/workflow-data/runs``. Project run artifacts and environment
evidence are stored under ``~/.oneflow/workflow-data/executions``.

Persisted Workbench run evidence
--------------------------------

For each Project Workbench run, the background runner captures the environment
from the exact selected executable before launching the Case or workflow. The
Case Workbench queue captures the same information once before launching its
selected Case batch. Snapshots are stored under the user data root, independently
of whether Case files run in the default workspace, a custom workspace, or in place.

The current layout is:

.. code-block:: text

   ~/.oneflow/
     workflow-data/
       executions/
         <run-id>/
           run-manifest.json
           environment.json
           <project-name>/             (Project Workbench default mode)
             ... project and Case files ...
           <case-id>/                  (Case Workbench default mode)
             ... staged Case files ...

``environment.json`` is the JSON document returned by
``OneFLOW --environment --format json``. It is written once at run start and is
not modified when the run finishes. If the environment query fails, the file
contains a diagnostic JSON object with ``capture_status: failed`` and the
error message; this does not prevent the CFD run from starting.

``run-manifest.json`` is created before execution and finalized afterward. It
records its schema version, Run ID, status, start and finish times, executable
path, command template, working directory, source project directory,
execution mode, relative environment report path, environment schema/query
status, and process return code. The Project Workbench run-history record
links to the absolute paths of both evidence files in its project-step
metadata. The Run History dialog resolves those links when a record is
selected and displays the parsed Environment Report and Run Manifest JSON in
the run-details pane. Legacy history records without evidence links remain
readable and are marked as not linked.

A retry creates a new Run ID and a new evidence directory. Existing snapshots
are not overwritten. The environment report is evidence, not configuration;
editing it does not change the executable or reproduce a previous runtime.

This persistence integration covers Project Workbench runs and Case Workbench
queues. A queue captures one environment snapshot for its selected executable
and Run ID; all Cases in that queue reference the same snapshot. Project
workflows may launch multiple Case processes and currently record the common
executable and command template. A future enhancement may record per-process
environment snapshots when process-specific launch environments differ.

External interface and schema policy
------------------------------------

The JSON command is the first supported external information interface. It
can be used by Python tools, scripts, diagnostics collectors, and other
applications without including OneFLOW or CGNS headers.

Consumers should validate required fields and ignore unknown fields. A
dedicated stable C ABI may be added later if embedding OneFLOW in third-party
applications becomes a concrete requirement. Until then, the CLI and its
versioned JSON schema avoid exposing internal C++ classes or ownership rules.

Dependency expansion policy
---------------------------

Additional dependencies should be registered through the same typed report,
not through GUI-specific probes. Each dependency entry should distinguish
runtime API version, loaded module path and linkage, compile-time/header
version, build configuration and enabled/disabled state, query status
with diagnostic detail when unavailable or failed, and a version-consistency
assessment.

The consistency assessment uses four states: verified_match, mismatch,
assumed_match, and unknown. OneFLOW compares independently obtained version
strings when both are present. When only a usable header version is available,
it defaults to assumed_match and records header_version as the basis. This
assumption must never be presented as an independent runtime observation or
as proof of binary/ABI compatibility. If only a runtime version is known, or
neither version is known, the assessment remains unknown.

HDF5 is now inspected without introducing an undeclared compile-time or
link-time dependency. At runtime OneFLOW searches the modules already loaded
into its process for the public ``H5get_libversion`` symbol. On Windows it uses the Windows Process Status API (PSAPI) to enumerate
loaded process modules, then resolves the symbol with ``GetProcAddress``.
The PSAPI calls are isolated inside the runtime symbol-discovery helper; the
report-building and JSON layers use platform-neutral C++ data structures.
The Windows executable links ``Psapi`` for ``EnumProcessModules``. On Linux
and macOS the helper uses the dynamic loader's default symbol scope
(``dlsym``), and obtains the containing module with ``dladdr``. These
platform-specific loader calls are kept behind the same helper boundary. If found,
OneFLOW calls the function and reports the returned major, minor, and release
numbers together with the module containing the symbol. If the symbol is not
exported by any loaded module, the report uses ``status: not_detected`` and
``runtime_version_status: runtime_api_symbol_not_found``. A successful lookup
reports the actual runtime API result; it does not infer HDF5 presence merely
because CGNS may use HDF5 internally.

This is a process-level symbol probe, not a guarantee that every private or
statically embedded HDF5 build can be identified. A static copy hidden inside
another library may not export ``H5get_libversion`` to the process symbol
lookup and will therefore remain not detected. OneFLOW does not add HDF5
headers or a direct HDF5 link dependency for this probe. The platform dynamic
loader library is linked where required.

MPI should use MPI_Get_library_version after MPI initialization or an
explicitly documented implementation-safe query path. METIS and JSON support
should report only evidence available from their public APIs or build
metadata. Compiler and accelerator-backend information should be reported
separately from third-party library versions.

Verification
------------

The C++ runtime_environment_test target verifies the JSON schema, CGNS and
HDF5 dependency entries, text-report output, and command dispatch without
running a simulation. The Python
test_runtime_environment_dialog suite verifies executable invocation, JSON
validation, process errors, missing executables, and timeout handling.

Actual dependency values are process- and installation-specific. HDF5 may be
reported as available only when its runtime API symbol is visible in the
process; otherwise it is explicitly not detected. Tests verify the response
contract and status semantics, not a developer machine's CGNS/HDF5 path or
version.
