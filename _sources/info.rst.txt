OneFLOW
=======

OneFLOW is an evolving open-source codebase for computational fluid dynamics (CFD) and scientific computing. This page is a short guide to the project and the documentation.

The repository includes a C++ runtime, mesh generation and conversion workflows, numerical methods, tests, and optional developer tools. Solver models, numerical schemes, platforms, and acceleration backends have different levels of maturity. Check the implementation, matching tests, and relevant runbooks for the revision you use rather than inferring support from a design goal or build option.

Choose a starting point
-----------------------

* **Build and run:** use the `README.md` in the `hexin` branch for current prerequisites, CMake build commands, tests, and command-line usage: `Repository README <https://github.com/eric2003/OneFLOW/blob/hexin/README.md>`_.
* **Installation:** see :doc:`install/index` for the existing installation documentation.
* **Mesh generation and conversion:** start with the `grid/` guide in the repository: `Grid case guide <https://github.com/eric2003/OneFLOW/blob/hexin/grid/README.md>`_.
* **CFD and grid regression:** use the `test/` guide for test suites and verification expectations: `Regression test guide <https://github.com/eric2003/OneFLOW/blob/hexin/test/README.md>`_.
* **Desktop Workbench:** see the `tools/oneflow_gui/` guide for the current project, case, and workflow tools: `Workbench guide <https://github.com/eric2003/OneFLOW/blob/hexin/tools/oneflow_gui/README.md>`_.

Architecture and design
-----------------------

* :doc:`mesh_architecture` — mesh summaries, snapshots, and format-adapter boundaries.
* :doc:`grid_topology` — topology semantics, boundary normalization, and derived connectivity.
* :doc:`workbench_architecture` — project/case concepts, execution modes, persistence, and process lifecycle.

Documentation is organized under the repository's `doc/` directory. The `doc/source/` tree contains the curated Sphinx site; engineering plans, architecture records, runbooks, evidence, and reports remain in their respective directories. See the `documentation consolidation plan <https://github.com/eric2003/OneFLOW/blob/hexin/doc/plans/oneflow-documentation-consolidation-plan.md>`_ for the inventory and staged cleanup status.
