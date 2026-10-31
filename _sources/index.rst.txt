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
* :doc:`grid_topology` — topology contracts, boundary normalization, and derived connectivity.
* :doc:`face_stage_pipeline` — staged face write path, BC export, reorder contract, identity vs orientation.
* `Grid case guide <https://github.com/eric2003/OneFLOW/blob/hexin/grid/README.md>`_ — mesh cases and conversion workflows.
* :doc:`workbench_user_guide_zh` — Chinese user guide for the optional Workbench GUI.
* :doc:`case_parameter_design` — case parameter files and GUI editing design.
* `Workbench module reference <https://github.com/eric2003/OneFLOW/blob/hexin/tools/oneflow_gui/README.md>`_ — detailed current behavior, paths, and service tests.

Developer resources
-------------------

* :doc:`test_inventory` — source-discovered test catalog, separated by framework, with input/expected/result contracts.
* :doc:`runtime_architecture` — runtime entry points, configuration boundaries, execution flow, and refactoring rules.
* :doc:`database_architecture` — parameter ownership, type-erased storage, const access, field pointer registration, and serialization contracts.
* :doc:`evolutionary_computation_architecture` — computation-unit contracts, intermediate representations, and the long-term algorithm/dataflow evolution model.
* :doc:`architecture_review` — whole-system CFD/CAE architecture checklist, major cross-cutting risks, and AI-assisted development review prompts.
* :doc:`architecture_priority_roadmap` — ranked investment sequence by system-wide return, verification cost, and cross-subsystem leverage.
* `AI project context and refactoring handoff <https://github.com/eric2003/OneFLOW/blob/hexin/doc/AI_CONTEXT.md>`_ — current focus, verified findings, open gaps, and next actions for continuing work.
* :doc:`workbench_architecture` — Workbench concepts, service boundaries, and run lifecycle.
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
   runtime_architecture
   database_architecture
   evolutionary_computation_architecture
   mesh_architecture
   grid_topology
   face_stage_pipeline
   workbench_architecture

.. toctree::
   :hidden:
   :maxdepth: 1

   case_parameter_design
