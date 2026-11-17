# OneFLOW Documentation

[English](README.md) | [简体中文](README_zh_CN.md) | [Project README](../README.md) | [Online documentation](https://oneflow-cfd.readthedocs.io/en/latest/)

The `doc/` directory is the repository's documentation root. It contains engineering records and the source files for the published Sphinx manual. These have different purposes:

- **Engineering documentation** supports development, review, handoffs, reproducibility, and verification.
- **`source/`** contains the curated Sphinx documentation site and its reader-oriented navigation.
- **Not every engineering record belongs on the public documentation site.** The manual should summarize stable concepts and link to detailed records instead of duplicating them.

English is the primary language for public project documentation and contributions. Chinese documentation is provided as a reading aid for selected topics; it should preserve the same concepts, terminology, and navigation without requiring every technical note to be translated line by line.

## Directory map

- `architecture/`: architecture descriptions, design boundaries, and review findings.
- `plans/`: maintained TODOs, roadmaps, improvement plans, and documentation-governance plans.
- `handoff/`: project handoffs, milestone summaries, and continuation notes.
- `runbooks/`: repeatable build, test, debug, cluster, and porting procedures.
- `evidence/`: verifiable build, MPI, and runtime evidence.
- `reports/`: formal reports organized by topic.
- `_drafts/`: temporary drafts to review, promote, or remove; not publication entry points.
- `source/`: Sphinx documentation source files.
- `Makefile`, `make.bat`: Sphinx build entry points.
- `DevelopManual-chinese.doc`, `UserManual-chinese.doc`: historical Word/WPS manuals awaiting content review.

## Documentation consolidation

A small number of technical documents still live in the repository root or the top-level `docs/` directory. Do not move, delete, or rename them in bulk before comparing content and checking links. First identify the authoritative version for each topic, then incorporate stable material into the Sphinx manual where appropriate.

## Read from the system down to implementation details

The documentation follows this progression:

1. **Project entry and goals.** Start with the repository [README](../README.md) and this page to understand OneFLOW's scope and documentation boundaries.
2. **System architecture.** Read the [runtime architecture and development guide](source/runtime_architecture.rst) for runtime execution, case lifecycle, mesh data, the Workbench, and verification boundaries.
3. **Subsystem architecture.** Follow the topic relevant to the task: mesh, topology, Workbench, or numerical components.
4. **Algorithms and interfaces.** Explain mathematical definitions, important equations, data flow, actual code entry points, and current implementation limits.
5. **Verification evidence.** Link contract tests, integration tests, mesh checks, and CFD regressions. Distinguish verified behavior from assumptions that still need validation.

Top-level documents should answer what the system contains, why it is layered this way, how its main control and data flows work, and where the major risks remain. Topic documents should explain how an algorithm works, what an interface promises, and where the implementation lives. Avoid duplicating function-level details in the system overview, and avoid technical notes that list code locations without explaining intent.

For a broad architecture review, read [system architecture review](source/architecture_review.rst) first, then the [runtime architecture guide](source/runtime_architecture.rst). Continue with [mesh architecture](source/mesh_architecture.rst), [grid topology](source/grid_topology.rst), or [Workbench architecture](source/workbench_architecture.rst), depending on the task. The review page is a checklist; it does not claim that every general capability discussed there has already been implemented.

For the normative [continuous evolution and scientific verification development specification](AI_TASK_TEMPLATE.md), use the English version or its [Chinese counterpart](AI_TASK_TEMPLATE_zh_CN.md). It defines single-workstream prioritization, implementation constraints, scientific verification and validation, references, reproducible plots, and delivery evidence.

For a quick continuation of ongoing development, read [AI project context and handoff](AI_CONTEXT.md). It summarizes branch constraints, recently inspected implementation baselines, confirmed behavior, open verification gaps, and suggested next steps. It is a navigation aid, not a substitute for source inspection and verification against the target commit.

For the initial inventory and planned consolidation of existing documents, see the [documentation consolidation plan](plans/oneflow-documentation-consolidation-plan.md). The plan records intended ownership and sequencing; it does not imply that every listed migration is complete.

New engineering documents belong under `doc/`; do not grow the top-level `docs/` directory into a second documentation system. Follow the [document naming rules](DOCUMENT_NAMING.md) for new files in `architecture/`, `plans/`, `handoff/`, `runbooks/`, `evidence/`, `reports/`, and `_drafts/`. Existing files should not be renamed in bulk because historical links may depend on their paths.

## Build the HTML documentation

From the repository root:

```bash
cd doc
make html
```

The output is written to `doc/build/html/`. CI and Read the Docs use `doc/source/` as the Sphinx source directory. When changing the published manual, check the relevant `toctree` entries, cross-references, and build output. An engineering note does not need to be added to the website navigation unless it is intended for publication.
