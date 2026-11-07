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

* runtime_version is populated only when the CGNS library linked into this
  build exposes cg_get_cgns_version and cg_get_cgns_version_string. CMake
  checks declaration and linkability against the configured CGNS headers and
  library before enabling the calls.
* header_version is derived from the CGNS_VERSION macro in the header used to
  compile OneFLOW. It is compile-time evidence, not proof of the version of a
  different library that the operating system may load.
* runtime_library.module_path and runtime_library.linkage identify the module
  containing the actual cg_open function address. On Windows this identifies
  the loaded DLL or the executable when CGNS is linked into it. On Linux and
  macOS it uses dynamic-loader module information. A static build may
  therefore report the OneFLOW executable as the containing module.

The public CGNS 4.5.2 header does not declare the two cg_get_cgns_version
functions. For that version the runtime API field is reported as unavailable,
while the header version and loaded-module information remain visible.
OneFLOW must not relabel the header macro as a runtime query result.

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
by itself. The Project Workbench now captures the selected executable's
environment at run start and stores it under
``~/.oneflow/workflow-data/executions/<run-id>/`` (or the equivalent
``ONEFLOW_HOME`` root). The snapshot is written even if the query fails: a
failed query produces a small JSON diagnostic document and does not prevent
the CFD project from running. Catalog Case-queue integration is a separate
follow-up and is not yet covered by this persistence slice. A user can still
save a standalone point-in-time report using the redirection examples above.

Workbench settings, workflow definitions, run history, and execution
workspaces have separate storage responsibilities under the OneFLOW user data
root (~/.oneflow by default; ONEFLOW_HOME can relocate that root). Catalog Case
runs staged by the Workbench use
~/.oneflow/workflow-data/executions/<run-id>/<project-name> by default. The runtime report is not
currently part of that run artifact contract.

Planned per-run evidence storage
--------------------------------

The Project Workbench persistence slice captures the environment from the exact
executable selected for a run and stores an immutable JSON snapshot under the
user data execution directory. The current layout is:

.. code-block:: text

   ~/.oneflow/
     executions/
       <run-id>/
         run-manifest.json
         environment.json
         <case-id>/
           input/
           output/
           logs/

For a multi-case or dependency workflow, environment.json belongs at the run
root when all Cases share one OneFLOW process. If Cases execute as separate
processes, each process must receive its own snapshot, or the manifest must
reference the snapshot captured for that exact process.

The run manifest should reference the environment snapshot by relative path
and record the executable path, command arguments, working directory, run
start time, and environment schema version. The snapshot is evidence, not
configuration: it must not be edited to change a past run. A retry creates a
new run ID and a new snapshot. Existing run-history records and artifacts
must not be deleted or rewritten by this addition.

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
version, build configuration and enabled/disabled state, and query status
with diagnostic detail when unavailable or failed.

HDF5 is the next candidate because H5get_libversion returns the version of
the HDF5 library linked to the application. It should only be queried when
OneFLOW has a direct, correctly configured HDF5 link dependency; CGNS using
HDF5 internally does not by itself prove that OneFLOW may safely call HDF5
symbols directly. MPI should use MPI_Get_library_version after MPI
initialization or an explicitly documented implementation-safe query path.
METIS and JSON support should report only evidence available from their public
APIs or build metadata. Compiler and accelerator-backend information should
be reported separately from third-party library versions.

Verification
------------

The C++ runtime_environment_test target verifies the JSON schema and command
dispatch without running a simulation. The Python
test_runtime_environment_dialog suite verifies executable invocation, JSON
validation, process errors, missing executables, and timeout handling.

Actual dependency values are process- and installation-specific. Tests must
verify the response contract, status semantics, and storage isolation; they
must not assert a developer machine's CGNS installation path or version.
