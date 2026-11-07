OneFLOW Runtime Environment and Dependency Information
========================================================

Purpose
-------

OneFLOW exposes runtime and build-environment information through the
``OneFLOW`` executable. The first implementation focuses on CGNS because the
actual CGNS library loaded by a process is important evidence when diagnosing
mesh conversion and file I/O behavior.

The Workbench does not link to CGNS. It asks the selected OneFLOW executable
for a versioned JSON document and displays the information returned by the
runtime. This keeps the GUI independent of the solver's third-party library
headers and ABI.

Command-line interface
----------------------

Run these commands from a terminal:

.. code-block:: console

   OneFLOW --version
   OneFLOW --environment
   OneFLOW --environment --format json

``--version`` prints the OneFLOW application version.
``--environment`` prints a human-readable report.
``--environment --format json`` writes a machine-readable document to
standard output. The JSON document has ``schema_version: 1`` so external
clients can validate the response before consuming it.

The informational commands are handled in the executable entry point before
the simulation lifecycle starts. They do not initialize a case, MPI process
environment, solver, or mesh.

CGNS version evidence
---------------------

The report deliberately separates three kinds of information:

* ``runtime_version`` is populated only when the CGNS library linked into this
  build exposes the ``cg_get_cgns_version`` and
  ``cg_get_cgns_version_string`` APIs. CMake checks both declaration and
  linkability against the configured CGNS headers and library before enabling
  the calls.
* ``header_version`` is derived from the ``CGNS_VERSION`` macro in the header
  used to compile OneFLOW. It is compile-time evidence, not proof of the
  version of a different library that the operating system may load.
* ``runtime_library.module_path`` and ``runtime_library.linkage`` identify the
  module containing the actual ``cg_open`` function address. On Windows this
  identifies the loaded DLL or the executable when CGNS is linked into it. On
  Linux and macOS it uses the dynamic-loader module information. A static
  build may therefore report the OneFLOW executable as the containing module.

The public CGNS 4.5.2 header does not declare the two
``cg_get_cgns_version*`` functions. For that version the runtime API field is
reported as unavailable, while the header version and loaded-module
information remain visible. OneFLOW must not relabel the header macro as a
runtime query result.

Do not confuse the library version with ``cg_version(fn, &version)``. The
latter queries the CGNS file-format version stored in an opened file; it does
not identify the version of the CGNS library used by the process.

Workbench integration
----------------------

The Workbench home screen provides **Runtime Environment...**. It invokes the
configured executable with ``--environment --format json`` and displays the
OneFLOW build identity, dependency status, runtime version (when the library
supports a runtime query), header version, and loaded module path.

The GUI treats the JSON as data returned by OneFLOW. It does not independently
load CGNS, infer versions from a configured installation directory, or
substitute a CMake cache value for a runtime result.

External integration
--------------------

The JSON command is the first supported external information interface. It
can be used by Python tools, scripts, diagnostics collectors, and other
applications without including OneFLOW or CGNS headers.

A dedicated stable C ABI may be added later if embedding OneFLOW in third-party
applications becomes a concrete requirement. Until then, the CLI and its
versioned JSON schema avoid exposing internal C++ classes or ownership rules.

Verification
------------

The C++ ``runtime_environment_test`` target verifies the JSON schema and
command dispatch without running a simulation. The Python
``test_runtime_environment_dialog`` suite verifies executable invocation,
JSON validation, process errors, missing executables, and timeout handling.

The actual dependency values are process- and installation-specific. A test
must verify the response contract, not assert a developer machine's CGNS
installation path or version.
