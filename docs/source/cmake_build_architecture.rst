CMake Build Architecture
========================

Purpose and scope
-----------------

OneFLOW's CMake infrastructure is part of the software architecture. It
defines how source files are discovered, how modules expose headers, how
targets acquire dependencies, and how changes are verified across platforms.
It must evolve without forcing a simultaneous source-tree migration.

This first stage introduces a reusable source-discovery contract and tests it
independently. It deliberately does not move geometry files, rename existing
targets, or change production link relationships.

Design principles
-----------------

Compatibility before migration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current tree contains legacy modules and modules using the
include/ + src/ convention. Both are valid inputs during the transition.
A module may add nested submodules without requiring every existing module
to be reorganized first.

Discovery is not dependency inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finding a .cpp file does not establish which libraries it needs. Source
discovery and target dependency declarations are separate responsibilities.
For example, a Gmsh reader may be independent of CGNS while a conversion
writer depends on CGNS. The build system must not infer that every geometry
component should link every optional dependency.

Production source boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The source-discovery helper excludes generated build, CMake-generated, and VCS
directories. It does not exclude directories merely because they are named
test or tests. The former codes/test module has been retired: JSON behavior is
covered by focused GoogleTest cases under tests/io, while Vencat analytic
profile generation is maintained under tools/verification. Test/production
boundaries are established by the selected module root and target registration,
not by guessing from a directory name. Generated files belong in the binary
tree, never in the source tree.

Deterministic results
~~~~~~~~~~~~~~~~~~~~~

Discovered source, header, and include-directory lists are deduplicated and
sorted. Stable ordering makes generated project files, diagnostics, and
architecture tests easier to compare.

Current module discovery contract
----------------------------------

cmake/modules/OneFLOWSourceDiscovery.cmake provides:

.. code-block:: cmake

   oneflow_discover_module_sources(
       MODULE_ROOT
       source_list
       header_list
       include_dir_list
   )

The function recursively discovers C++ implementation files (.cpp) and
headers (.h, .hpp). Header-containing directories are returned as
include directories. It supports legacy implementation directories and
the conventional include/ + src/ layout, including nested source folders.

Directories named build, CMakeFiles, or .git are excluded. Directories named
test or tests remain discoverable when they are inside the selected module
root. The fixture locks this generic rule independently of the current
production module layout.
Implementation files placed under a directory named include are not treated
as production sources.

The helper uses CONFIGURE_DEPENDS so additions and removals are noticed by
CMake's generated build-system checks. A configure-time test also explicitly
reconfigures a fixture project after adding and removing files.

The production runtime source list in codes/CMakeLists.txt now uses
oneflow_discover_module_sources. The legacy ConstructSolutionDirTree function
is retained as an independent compatibility oracle for the parity test, not
as the production discovery path. Visual Studio source-group hierarchy is
applied separately after discovery, keeping IDE presentation out of the
source enumeration contract.

Testing the build system
------------------------

The CMakeSourceDiscovery CTest runs a child CMake project with
project(... NONE). It does not need a C++ compiler and does not build
OneFLOW. The fixture is created below the CMake binary directory and verifies:

* simultaneous discovery of legacy and include/ + src/ layouts;
* recursive discovery in nested source directories;
* header and include-directory discovery;
* preservation of valid sources under a directory named test;\n* exclusion of generated build and CMakeFiles sources;
* discovery of a source added after the initial configure;
* removal of a source deleted before a later configure.

Run it with:

.. code-block:: console

   ctest --test-dir build -R CMakeSourceDiscovery --output-on-failure

The CMakeLegacyDiscoveryParity CTest additionally invokes both the existing
ConstructSolutionDirTree implementation and the new discovery helper against
the real codes/ tree. It compares sorted source, header, and include-directory
sets. This is a compatibility gate: production callers must not migrate until
this comparison passes.

Both tests are independent of GoogleTest and third-party CFD libraries. This
makes failures attributable to the CMake contract rather than to compiler,
CGNS, METIS, or solver configuration.

Evolution plan
--------------

The infrastructure should evolve in controlled steps:

#. Establish and test source-discovery semantics while retaining the
   legacy implementation as a compatibility oracle.
#. Migrate the production runtime source enumeration to the tested helper,
   preserving IDE source groups and comparing source, header, and
   include-directory sets.
#. Adapt additional modules one at a time, comparing old and new discovery
   results before changing target construction.
#. Introduce target-level module boundaries only where dependency direction
   and ownership are understood.
#. Move STL, Gmsh, CGNS, and conversion components into clearer submodules
   after their build and test dependencies are explicit.
#. Validate Windows Visual Studio and Linux Ninja/Makefiles configurations,
   then run the full CTest and Python suites.

A successful configure is not sufficient evidence for a migration. The
accepted source set, target link interface, optional dependency behavior,
generated files, and regression results must all remain observable.
