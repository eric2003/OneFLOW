# OneFLOW

**An open-source CFD codebase for mesh generation, flow simulation, and reproducible engineering workflows.**

OneFLOW is an actively developed computational fluid dynamics (CFD) project. Its goal is not only to run flow cases, but also to make the surrounding engineering work clearer and more reproducible: define a mesh or case, run it in a controlled workspace, inspect generated artifacts, validate results, and connect cases when one calculation depends on another.

This repository is evolving. This README describes the code and tools that exist today, distinguishes working features from longer-term direction, and provides build and test instructions for the current CMake project. OneFLOW is the open-source CFD and scientific-computing project in this repository; it is unrelated to similarly named commercial contract-lifecycle-management products.

[简体中文](./README_zh_CN.md) · [Online documentation](https://oneflow-cfd.readthedocs.io/en/latest/) · [Source code](https://github.com/eric2003/OneFLOW) · [Issues](https://github.com/eric2003/OneFLOW/issues)

## Project vision

OneFLOW is developed as a general-purpose, unstructured-mesh CFD code and a scientific-computing engine for broader multidisciplinary research. Its scope includes PDE computation across one, two, and three dimensions, learning and validating numerical methods, solving engineering flow cases, and providing a modular codebase that researchers can inspect and extend. The equations, mesh types, and physical models available in practice depend on the solver and case implemented in the current revision. The longer-term vision is to make the scientific-computing workflow more programmable and automated: from defining a problem and running a case to checking results, recording evidence, and composing dependent tasks. The phrase *scientific-computing operating system* describes this long-term direction, not a claim that a complete operating-system-like platform already exists.

The repository contains capabilities at different maturity levels. Build options or names for MPI, CUDA, HIP, or other backends do not by themselves establish production readiness, scalability to a particular core count, or validated performance on a given device. Check the implementation, configuration output, and matching test evidence before relying on such capabilities.

## What is in the project?

- **CFD runtime:** a C++20 application with solver, task, data, geometry, and I/O code organized into modules.
- **Mesh and grid workflows:** mesh conversion and grid-generation cases, structured and unstructured mesh work, boundary and connectivity handling, and tools for inspecting mesh metadata.
- **Case-based verification:** GoogleTest/CTest unit tests, existing CFD regression cases, and grid-conversion checks. A successful process exit alone is not sufficient evidence for every regression case.
- **Engineering workbench:** an optional Python/PySide6 desktop application for opening projects, running cases, editing mesh layouts, managing multi-case workflows, and reviewing logs and run history.
- **Incremental backend work:** a configurable compute-backend boundary, with CPU as the practical baseline. Accelerator integrations have different levels of readiness; do not assume every backend is production-ready.

## Design principles

### Keep the solver independent of the workbench

The C++ runtime is the authority for executing OneFLOW cases. Python tools and the desktop workbench orchestrate that runtime; they should not become hidden prerequisites for the solver or duplicate its numerical implementation.

### Treat cases as projects and runs as artifacts

Inputs and case definitions belong in the source repository. Generated grids, reports, logs, and flow results should live in a run workspace. The workbench stages built-in cases away from the source tree to avoid accidentally overwriting checked-in inputs.

### Separate data extraction, representation, and presentation

Mesh and geometry work is moving toward small, independently testable interfaces: extract a summary from a runtime mesh, represent the snapshot without depending on a GUI, then serialize or visualize it. A mesh summary is metadata (for example, dimensions, counts, and bounds), not a replacement for the mesh itself.

### Make dependencies explicit

Some cases are independent; others consume a grid or artifact produced by an earlier case. The workflow model records execution mode and dependency edges so orchestration can validate and order work instead of relying on undocumented manual steps. Some workbench paths are still being integrated; consult its dedicated guide for current behavior.

### Prefer evidence over assumptions

CTest, unit tests, mesh expectations, and CFD regression baselines provide different kinds of evidence. They are not interchangeable. Hardware-dependent accelerator validation is separate from the normal CPU build and should be reported as such.

## Architecture at a glance

OneFLOW is best understood as a layered engineering codebase rather than a solver executable surrounded by unrelated scripts.

- **Numerical runtime:** C++ owns the actual case execution and numerical work. Command-line use remains possible without starting the desktop workbench.
- **Mesh data pipeline:** distinguish the live runtime mesh from a stable, independently testable snapshot and from output formats such as JSON. A summary describes metadata such as dimension, entity counts, and bounds; it is not the mesh itself. Geometry, topology, connectivity, and boundary information should be inspectable through focused interfaces rather than one all-purpose exporter.
- **Case and artifact boundary:** case definitions and reference inputs are versioned; generated meshes, flow results, logs, and reports belong to a run workspace. Explicit input/output roots allow a producer case and a consumer case to exchange artifacts without relying on accidental working-directory state.
- **Workflow orchestration:** a case may run alone, in a sequence, independently alongside other cases, or as part of a dependency graph. Dependencies should name the artifact-producing relationship and be validated before execution, rather than being hidden in manual steps.
- **Tooling and presentation:** CLI utilities, Python services, and the optional PySide6 workbench discover and orchestrate work, inspect outputs, and present results. They should reuse shared case/workflow contracts where practical, not duplicate solver behavior or make the runtime depend on GUI code.
- **Testable evolution:** unit tests, mesh invariants, and numerical regression tests cover different contracts. Refactoring should follow real call paths, make ownership and interface boundaries clearer, and preserve compatibility where required; tests are evidence for behavior, not a reason to rewrite code without a concrete need.

These are architectural boundaries and an ongoing direction, not a claim that every workflow or mesh interface has already reached its final form. The case catalog and profiles, for example, are not yet the authoritative source for every legacy runtime parameter; see the [grid guide](./grid/README.md) for current migration boundaries. For a code-backed map of runtime entry points, configuration layers, task dispatch, and refactoring rules, see the [runtime architecture and development guide](./doc/source/runtime_architecture.rst).

## Repository map

| Path | Purpose |
| --- | --- |
| codes/ | C++ runtime, solver, geometry, mesh, and I/O implementation |
| ThirdParty/ | Vendored dependencies used by the build, including GoogleTest and JSON libraries |
| tests/ | C++ unit tests registered with CTest |
| grid/ | Mesh-generation and mesh-conversion cases, profiles, and case catalog |
| test/ | CFD/grid regression runners, suite manifests, and numerical baselines |
| tools/ | Case validation, mesh utilities, and the optional desktop workbench |
| doc/ | Official documentation root: architecture, plans, runbooks, evidence, reports, and Sphinx sources |
| docs/ | Existing supporting notes; new project documentation belongs under `doc/` |
| cmake/ | CMake helper functions |

Start with the [online documentation](https://oneflow-cfd.readthedocs.io/en/latest/). For repository-specific details, see the [grid case guide](./grid/README.md), [regression test guide](./test/README.md), [OneFLOW Workbench guide](./tools/oneflow_gui/README.md), and [documentation organization plan](./doc/plans/oneflow-documentation-consolidation-plan.md). For AI-assisted development or a fast handoff between sessions, start with the [AI project context and current refactoring state](./doc/AI_CONTEXT.md), then follow its links to architecture notes, acceptance criteria, tests, and CI evidence.

## Build from source

### Requirements

- CMake **3.30 or newer**.
- A C++20-capable compiler.
- Git.
- MPI, METIS, and CGNS development libraries for the default build configuration.
- HDF5 is normally required by CGNS packages.
- Python 3 is useful for case tools, regression tests, and the optional desktop workbench.

Dependency discovery differs by platform. The current build does not use a single universal package-manager configuration, so verify the include and library paths that CMake reports when configuring.

### Windows: Visual Studio 2022

Install Visual Studio 2022 with the **Desktop development with C++** workload, CMake, and Git. Install compatible MPI, METIS, and CGNS libraries. The current CMake defaults expect MS-MPI under its standard SDK location, METIS under C:/dev/METIS/METIS-VS2022-STATIC, and CGNS under C:/dev/cgns/4.4.0. Override these paths for your installation.

From PowerShell at the repository root:

~~~powershell
cmake -S . -B build `
  -G "Visual Studio 17 2022" -A x64 `
  -D CMAKE_INSTALL_PREFIX="C:/dev/OneFLOW" `
  -D METIS_ROOT="C:/dev/METIS/METIS-VS2022-STATIC" `
  -D CGNS_ROOT="C:/dev/cgns/4.4.0"

cmake --build build --config Debug --parallel
ctest --test-dir build -C Debug --output-on-failure
~~~

For Release, use --config Release for both build and test commands. Visual Studio is a multi-configuration generator; select Debug or Release at build/test time rather than relying on CMAKE_BUILD_TYPE.

The executable is typically under build/bin/Debug/ or build/bin/Release/, depending on the selected configuration.

### Windows: Ninja / single-configuration build

If Ninja and the MSVC developer environment are available, a single-configuration build can be convenient:

~~~powershell
cmake -S . -B build-ninja -G Ninja `
  -D CMAKE_BUILD_TYPE=Debug `
  -D CMAKE_INSTALL_PREFIX="C:/dev/OneFLOW" `
  -D METIS_ROOT="C:/dev/METIS/METIS-VS2022-STATIC" `
  -D CGNS_ROOT="C:/dev/cgns/4.4.0"

cmake --build build-ninja --parallel
ctest --test-dir build-ninja --output-on-failure
~~~

Use a separate build directory when changing generators or compilers. CMake does not support switching the generator in an existing build tree.

### Linux: GCC and system libraries

On Debian/Ubuntu, install the compiler and common development packages:

~~~bash
sudo apt update
sudo apt install build-essential cmake git \
  libhdf5-dev libcgns-dev libmetis-dev libmpich-dev
~~~

The current CMake configuration accepts explicit library paths through environment variables. The following paths match a common Ubuntu/MPICH layout; adjust them for your distribution:

~~~bash
export MPI_HOME_INC=/usr/include/x86_64-linux-gnu/mpich
export MPI_HOME_LIB=/usr/lib/x86_64-linux-gnu/libmpich.so
export CGNS_HOME_INC=/usr/include
export CGNS_HOME_LIB=/usr/lib/x86_64-linux-gnu/libcgns.so
export METIS_HOME_INC=/usr/include
export METIS_HOME_LIB=/usr/lib/x86_64-linux-gnu/libmetis.so

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_INSTALL_PREFIX="$PWD/Run"
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure
~~~

If your distribution installs HDF5 or MPI headers in a different directory, update the corresponding paths. If CMake configuration or linking fails, inspect the detected include directories and libraries before changing source code.

### Run a case

The command-line entry point uses this case-directory form:

~~~text
OneFLOW.exe r <case-directory>
~~~

For example, after building on Windows:

~~~powershell
.\build\bin\Debug\OneFLOW.exe r .\grid\minimal_rectangle
~~~

Run cases in a writable workspace rather than directly in the checked-in source case when they generate files. The workbench can create isolated run copies for supported task types. Case inputs and outputs are described in the relevant case directory and the [grid guide](./grid/README.md).

## Tests and verification

### C++ unit tests

CTest discovers the GoogleTest-based C++ tests configured by CMake:

~~~powershell
# Visual Studio generator
ctest --test-dir build -C Debug --output-on-failure
~~~

For a single-configuration generator, omit -C Debug and use the build directory configured above. When investigating a test failure, run the relevant test executable or filtered GoogleTest case as well as the complete CTest suite.

### CFD and grid regression

The separate test/ tree contains regression manifests, case inputs, and accepted numerical baselines. See [test/README.md](./test/README.md) for supported commands, suite names, expected outputs, and tolerance policy. Grid case discovery and profile validation are documented in [grid/README.md](./grid/README.md).

Do not treat these checks as equivalent:

- **Unit tests** check individual contracts and components.
- **Grid checks** validate generated mesh properties and expected topology.
- **CFD regression** compares solver outputs with accepted numerical references.
- **Accelerator tests** require the appropriate compiler, runtime, and hardware; they are not implied by a CPU-only pass.

## Optional desktop workbench

The desktop application is a Python/PySide6 tool layered over the C++ executable. It is not required to build or run the command-line solver.

~~~powershell
python -m pip install -r tools/requirements-gui.txt
python -m tools.oneflow_gui
~~~

The workbench supports project discovery, run configuration, mesh layout editing, test entry points, multi-case/workflow definitions, logs, and run history. Its feature set is evolving; the detailed [workbench README](./tools/oneflow_gui/README.md) documents current behavior and known integration boundaries.

## Compute backends: current status

The CMake cache exposes ONEFLOW_ACCEL_BACKEND with AUTO, CPU, HIP, CUDA, and KOKKOS values, but these names do **not** mean that all implementations are equally complete.

- **CPU** is the baseline for ordinary development and regression.
- **HIP** can be enabled explicitly when the HIP package, compiler, and target architecture are available. Hardware-dependent validation must be performed on a matching system.
- **CUDA** follows a legacy opt-in source path and requires ONEFLOW_ENABLE_CUDA=ON; it should not be assumed to be enabled by default.
- **Kokkos** is currently a reserved extension point, not an enabled port.

Example CPU configuration:

~~~bash
cmake -S . -B build-cpu -D ONEFLOW_ACCEL_BACKEND=CPU
~~~

For HIP, consult the CMake output and set the compiler and architecture explicitly when cross-compiling. Do not enable an accelerator option unless the required toolchain and libraries are installed.

## Contributing

Contributions are welcome, especially focused changes with tests and reproducible evidence. A useful contribution usually:

1. Identifies the runtime or tool boundary it changes.
2. Adds or updates a test for the intended contract.
3. Keeps generated artifacts out of source case directories.
4. Runs the relevant unit or regression checks and reports any hardware or platform limitations.
5. Updates documentation when commands, parameters, or observable behavior change.

The project is actively refactoring its mesh, geometry, test, and workbench boundaries. Prefer incremental changes that improve a real call path or make behavior easier to verify over broad rewrites without a concrete use case.

## Developer and contact

**Primary developer:** Dr. He Xin ([GitHub profile: @eric2003](https://github.com/eric2003)).

- **Email:** [fantasy_2003_@hotmail.com](mailto:fantasy_2003_@hotmail.com)
- **Source code:** [github.com/eric2003/OneFLOW](https://github.com/eric2003/OneFLOW)
- **Bug reports and feature requests:** [GitHub Issues](https://github.com/eric2003/OneFLOW/issues)
- **Design and implementation discussion:** open an issue with the relevant case, expected behavior, observed behavior, and reproducible logs or test results. This helps keep technical decisions and follow-up work visible to other contributors.

Email is also available for direct contact. The address above was already published in an earlier version of the project README.

## Project and license

OneFLOW is maintained by He Xin and contributors. Please check the repository's license file and contribution history for the terms and process that apply to your intended use.

---

*OneFLOW is an evolving open-source CFD project. Documentation describes the checked-in implementation and its current limitations; it is not a claim that every planned capability is complete.*
