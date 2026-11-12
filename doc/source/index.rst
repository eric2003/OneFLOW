OneFLOW Scientific Computing
==============================

OneFLOW is an open-source codebase for computational fluid dynamics (CFD) and scientific computing. It brings together a C++ runtime, mesh generation and conversion, numerical methods, case-based verification, and engineering workflow tools.

This documentation serves three audiences: people getting OneFLOW built and running, users preparing and checking cases, and developers working on the runtime or engineering tools. OneFLOW is evolving; the code, tests, and runbooks for the revision in use determine which capabilities are currently supported. Roadmap items and build options alone do not establish implementation completeness.

Project links
-------------

* `Repository and source code <https://github.com/eric2003/OneFLOW/>`_
* `Online documentation <https://oneflow-cfd.readthedocs.io/en/latest/>`_
* `Issues and bug reports <https://github.com/eric2003/OneFLOW/issues>`_
* `Scientific computing resource collection <https://eric2003.github.io/Prometheus/>`_

Getting started
---------------

Start with the project overview, then follow the installation and build instructions for your platform.

* :doc:`info` — project scope and documentation map.
* :doc:`install/index` — installation and environment setup.
* `Repository README <https://github.com/eric2003/OneFLOW/blob/hexin/README.md>`_ — current build, test, and command-line instructions.

Using OneFLOW
-------------

* :doc:`user` — existing articles and scientific-computing resources.
* :doc:`cmd` — command-line and tool reference.
* :doc:`cfd/index` — CFD topics and numerical application references.
* :doc:`test` — test and verification references.

Mesh and engineering workflows
------------------------------

* :doc:`mesh_architecture` — mesh information levels, snapshots, and format adapters.
* :doc:`stl_mesh_io` — native C++ STL ASCII/Binary I/O, facet semantics, and reader/writer verification.
* :doc:`gmsh_interface` — Gmsh MSH 4.1-first input adapter, version support, semantic mapping, and conversion validation.
* :doc:`stl_to_gmsh_2d` — planar STL boundary extraction, explicit CFD boundary groups, and Gmsh 2D generation.
* :doc:`stl_to_gmsh_3d` — cube STL volume reconstruction, fixed HEXA_8/TETRA_4 topology, and automatic tetrahedral meshing.
* :doc:`stl_to_gmsh_sphere` — unit-sphere STL validation, Gmsh surface reconstruction, tetrahedral volume generation, commands, and evidence.
* :doc:`boundary_condition_architecture` — physical boundary semantics, Gmsh/CGNS/OneFLOW mappings, required state data, and verification boundaries.
* :doc:`grid_conversion_cases` — runnable mesh conversion Cases, topology coverage, and format-route matrix.
* :doc:`development_case_methodology` — how to design, implement, document, and verify reusable 2D/3D Development Cases.
* :doc:`grid_topology` — topology contracts, boundary normalization, and derived connectivity.
* :doc:`face_stage_pipeline` — staged face write path, BC export, reorder contract, identity vs orientation.
* `Grid case guide <https://github.com/eric2003/OneFLOW/blob/hexin/grid/README.md>`_ — mesh cases and conversion workflows.
* :doc:`workbench_user_guide_zh` — Chinese user guide for the optional Workbench GUI.
* :doc:`case_parameter_design` — case parameter files and GUI editing design.
* :doc:`case_parameter_resolution` — file-independent parameter overrides and Case construction plans.
* :doc:`engineering_nomenclature` — source-backed cross-software CFD terminology, variable aliases, units, and semantic caveats.
* `Workbench module reference <https://github.com/eric2003/OneFLOW/blob/hexin/tools/oneflow_gui/README.md>`_ — detailed current behavior, paths, and service tests.

Developer resources
-------------------

* :doc:`test_inventory` — source-discovered test catalog, separated by framework, with input/expected/result contracts.
* :doc:`cmake_build_architecture` — module source discovery, compatibility rules, dependency boundaries, and CMake self-tests.
* :doc:`runtime_architecture` — runtime entry points, configuration boundaries, execution flow, and refactoring rules.
* :doc:`runtime_environment` — runtime dependency introspection, CGNS version evidence, CLI JSON, and Workbench display.
* :doc:`runtime_case_lifecycle` — production CLI-to-case execution sequence, configuration ingestion, cleanup ownership, and integration evidence.
* :doc:`database_architecture` — parameter ownership, type-erased storage, const access, field pointer registration, and serialization contracts.
* :doc:`evolutionary_computation_architecture` — computation-unit contracts, intermediate representations, and the long-term algorithm/dataflow evolution model.
* :doc:`architecture_review` — whole-system CFD/CAE architecture checklist, major cross-cutting risks, and AI-assisted development review prompts.
* :doc:`architecture_priority_roadmap` — ranked investment sequence by system-wide return, verification cost, and cross-subsystem leverage.
* :doc:`verification_evidence_lifecycle` — how claims, tests, runs, artifacts, numerical results, and references form an auditable evidence chain.
* `AI project context and refactoring handoff <https://github.com/eric2003/OneFLOW/blob/hexin/doc/AI_CONTEXT.md>`_ — current focus, verified findings, open gaps, and next actions for continuing work.
* :doc:`workbench_architecture` — Workbench concepts, service boundaries, run lifecycle, and shared Parameter contracts.
* :doc:`contract_workbench_information_architecture` — shared Contract Workbench navigation with typed Formula, Case, Runtime, Workflow, Mesh, Geometry, and Algorithm knowledge.
* `Verification Contract Workbench proposal <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/verification-contract-workbench.md>`_ — proposed shared model for numerical, geometry, runtime, workflow, and case evidence; not yet an implementation claim.
* `Regression test guide <https://github.com/eric2003/OneFLOW/blob/hexin/test/README.md>`_ — CFD and grid regression checks.
* `Documentation organization plan <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/oneflow-documentation-consolidation-plan.md>`_ — document ownership and consolidation status.
* `Runtime refactoring roadmap <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/runtime-refactoring-roadmap.md>`_ — staged goals, acceptance criteria, and implementation evidence.

Reference library
-----------------

The existing reference library covers programming languages, build systems, parallel computing, numerical methods, visualization, and related tools. It remains available below; the curated entry points above are intended to help readers find a starting point without replacing that material.

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Get started

   /info
   /install/index

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: User Guides

   /user
   /cmd
   /test
   /workbench_user_guide_zh

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Reference

   /bash
   /cfd/index
   /cgns
   /chatgpt
   /cpp/index
   /cmake
   /cuda
   /design_patterns
   /eigen
   /git
   /fortran
   /hdf5
   /hpc
   /julia
   /jupyter
   /latex
   /makefile
   /metis
   /mkdocs
   /modules
   /mpi
   /numerical/index
   /oneapi
   /oneflow
   /openacc
   /opencl
   /opengl
   /openmp
   /powershell
   /python/index
   /qt
   /regex
   /science
   /sphinx
   /system/index
   /taichi
   /tcltk
   /tecplot
   /unreal_engine
   /vcpkg
   /yaml

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: OneFLOW Design and Development

   architecture_review
   cmake_build_architecture
   verification_evidence_lifecycle
   runtime_architecture
   runtime_case_lifecycle
   database_architecture
   evolutionary_computation_architecture
   mesh_architecture
   stl_mesh_io
   gmsh_interface
   stl_to_gmsh_2d
   stl_to_gmsh_3d
   stl_to_gmsh_sphere
   grid_conversion_cases
   development_case_methodology
   grid_topology
   boundary_condition_architecture
   face_stage_pipeline
   workbench_architecture
   mesh_to_cfd_validation
   contract_workbench_information_architecture

.. toctree::
   :hidden:
   :maxdepth: 1

   case_parameter_design
   engineering_nomenclature
