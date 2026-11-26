# OneFLOW Documentation

[English](README.md) | [Simplified Chinese](README_zh_CN.md) | [Project README](../README.md) | [Online documentation](https://oneflow-cfd.readthedocs.io/en/latest/)

The `docs/` directory is the documentation root for OneFLOW. It contains
engineering records and the source of the published Sphinx manual. These
serve different audiences: engineering records support implementation,
review, handoff, reproducibility, and verification; `source/` presents
stable concepts in a reader-oriented manual. Not every engineering record
belongs on the public documentation site.

English is the canonical baseline for public project documentation and
contributor-facing guidance. The Chinese index is a companion: keep its
navigation, technical meaning, and policy aligned, while retaining useful
Chinese-specific explanatory context where it adds value. When the two
versions conflict, follow the English version and update the Chinese version
to match. Neither index requires a sentence-by-sentence translation of every
engineering record.

## 1. Directory structure and retention

| Directory | Contents | Maintenance policy |
|---|---|---|
| [`architecture/`](architecture/) | Stable architecture descriptions and boundary reviews | Keep one current document per topic. Delete superseded reviews; Git history is the archive. |
| [`plans/`](plans/) | Living roadmaps and long-running plans | Keep living documents without date suffixes. Delete completed one-off checklists. Upstream architecture snapshots mapped into this tree are documented below. |
| [`handoff/`](handoff/) | Living project handoff | Maintain one `oneflow-current-handoff.md` in place; do not create dated handoff copies. |
| [`runbooks/`](runbooks/) | Repeatable build, cluster, porting, and operational procedures | Keep one maintained runbook per process; update it in place. |
| [`evidence/`](evidence/) | Reviewable environment and build facts | Keep evidence tied to a specific conclusion; remove it when that conclusion is superseded. |
| [`reports/`](reports/) | Formal performance, regression, and delivery reports | Do not delete published reports. Mark obsolete reports as historical; keep Markdown and HTML files with the same stem synchronized. |
| [`source/`](source/) | Sphinx manual sources for Read the Docs | Maintain according to the reader-facing manual structure. |
| `_drafts/` | Unfinished drafts | Move into the proper destination or delete after review; do not use as a permanent archive. |

Document naming rules are defined in [DOCUMENT_NAMING.md](DOCUMENT_NAMING.md).
Choose the destination and artifact type before writing. Adding, deleting, or
renaming a document requires updating this index and all affected references.

## 2. Current entry points

| Topic | Document |
|---|---|
| Phase A–F roadmap and stable identifiers | [Development plan](plans/oneflow-development-plan.md) |
| Current work slice, acceptance steps, and evidence | [Development TODO](plans/oneflow-development-todo.md) |
| Living handoff: environment, fixed paths, startup order, and pitfalls | [Current handoff](handoff/oneflow-current-handoff.md) |
| Fixed Kunshan workspace paths, job rules, and cleanup | [Kunshan workspace runbook](runbooks/oneflow-kunshan-workspace-runbook.md) |
| DCU porting method: contract → oracle → gate → benchmark | [DCU porting playbook](runbooks/oneflow-dcu-porting-playbook-20260903.md) |
| Main solver accelerator architecture and integration boundaries | [DCU architecture audit](architecture/oneflow-1d-euler-dcu-architecture-audit-20260904.md) |
| Published reports | [Reports index](reports/README.md) |

## 3. Fork and upstream document ownership

The repository contains both fork-maintained documents and documents shared
with upstream. Before editing, use these ownership rules:

- **Fork-maintained:** `plans/oneflow-development-*.md`, `handoff/`, this
  index, and fork-specific additions to `DOCUMENT_NAMING.md`.
- **Shared with upstream:** `architecture/`, `evidence/`, `reports/`,
  operational documents under `runbooks/`, and `source/`. Content may be
  improved, but avoid broad moves or rewrites that create structural merge
  conflicts. If a move is necessary, make it a separate, deliberate change and
  update references in `ports/`, `ci/`, and elsewhere.
- Do not recreate upstream's `doc/reports/architecture/` tree in this fork.
  The fork maps those documents into its established structure.

When merging new upstream files from `doc/reports/architecture/`, use these
mappings and complete the move and reference updates in the same merge slice:

| Upstream path | Fork destination |
|---|---|
| `reports/architecture/oneflow-dcu-porting-playbook-*.md` | `runbooks/` (same filename) |
| `reports/architecture/oneflow-gpu-backend-delivery.*` | `reports/delivery/` (same filename) |
| `reports/architecture/oneflow-euler-optimization-plan.md` | `plans/` (same filename) |
| `reports/architecture/<other>.md` (design notes or snapshots) | `plans/` (same filename) |

## 4. Reading order and technical entry points

For a broad architecture review, start with [the system architecture review](source/architecture_review.rst), then read [the runtime architecture guide](source/runtime_architecture.rst). Continue with [mesh architecture](source/mesh_architecture.rst), [grid topology](source/grid_topology.rst), or [Workbench architecture](source/workbench_architecture.rst), depending on the task.

The normative [continuous evolution and scientific verification specification](AI_TASK_TEMPLATE.md) defines workstream prioritization, implementation constraints, scientific verification and validation, references, reproducible plots, and delivery evidence. Its [Chinese counterpart](AI_TASK_TEMPLATE_zh_CN.md) is maintained as a companion.

For a quick continuation of ongoing development, read [AI project context and handoff](AI_CONTEXT.md). It summarizes branch constraints, inspected implementation baselines, confirmed behavior, open verification gaps, and suggested next steps. It is a navigation aid, not a substitute for source inspection and verification against the target commit.

For the initial inventory and planned consolidation of existing documents, see the [documentation consolidation plan](plans/oneflow-documentation-consolidation-plan.md). A plan records intended ownership and sequencing; it does not imply that every listed migration is complete.

New engineering documents belong under `docs/`. Follow [the document naming rules](DOCUMENT_NAMING.md) for files in `architecture/`, `plans/`, `handoff/`, `runbooks/`, `evidence/`, `reports/`, and `_drafts/`. Do not bulk-rename existing files: historical links may depend on their paths.

For topic-oriented reading, use this progression:

1. Project goals and documentation boundaries: the repository [README](../README.md) and this index.
2. System architecture: [runtime architecture and development guide](source/runtime_architecture.rst).
3. Subsystem architecture: mesh, topology, Workbench, or numerical component guides.
4. Algorithms and interfaces: mathematical definitions, equations, data flow, implementation entry points, and known limits.
5. Verification evidence: contract tests, integration tests, mesh checks, and CFD regression. Distinguish verified behavior from assumptions still awaiting validation.

The architecture review is a checklist, not a claim that every discussed capability has already been implemented. Top-level documents should explain system structure, rationale, major control/data flows, and risks; topic documents should explain algorithms, interface contracts, and implementation locations. Avoid duplicating function-level details in system overviews.

## 5. Build the HTML manual

From the repository root:

```bash
cd docs
make html
```

The output is written to `docs/build/html/`. CI and Read the Docs use
`docs/source/` as the Sphinx source directory. When editing published manual
pages, validate their toctrees, cross-references, and build output. Add an
engineering record to the website navigation only when it is intended for
publication.

Formal reports submitted in both Markdown and HTML must use the same filename
stem (for example, `oneflow-euler-performance-current.md` and
`oneflow-euler-performance-current.html`) and be updated together.

## 6. End-of-slice documentation checks

1. Update the checked items, status snapshot, and evidence in
   `plans/oneflow-development-todo.md`.
2. Update the task pointer and new pitfalls in
   `handoff/oneflow-current-handoff.md` in place. Do not add dated handoffs;
   remove temporary checklists and drafts when finished.
3. After adding, deleting, or renaming documents, validate links in this
   index, `docs/reports/README.md`, and `AGENTS.md`. Run
   `python3 scripts/check-doc-paths.py`; it should report no broken paths or
   obsolete `doc/` references.
4. Run `git diff --check`. Do not publish raw CI/Slurm logs, hostnames,
   account identifiers, job IDs, or private absolute paths in documentation.

## 7. Documentation layers and related contracts

The documentation tree contains two layers:

- **Engineering documentation** supports development, review, handoff,
  reproducibility, and verification, following the lifecycle-oriented
  directory policy above.
- **`source/`** is the reader-facing Sphinx manual. It should summarize
  stable concepts and link to detailed engineering records rather than copy
  them wholesale.

Useful related entry points:

| Topic | Document |
|---|---|
| Case lifecycle and incremental verification contract | [Case lifecycle contract](architecture/oneflow-case-lifecycle-contract.md) |
| System architecture review | [Architecture review](source/architecture_review.rst) |
| Runtime architecture | [Runtime architecture](source/runtime_architecture.rst) |
| Mesh, grid topology, and Workbench architecture | [Mesh architecture](source/mesh_architecture.rst), [grid topology](source/grid_topology.rst), [Workbench architecture](source/workbench_architecture.rst) |
| AI collaboration specification (English/Chinese) | [AI task template](AI_TASK_TEMPLATE.md), [Chinese counterpart](AI_TASK_TEMPLATE_zh_CN.md) |
| AI context and handoff | [AI context](AI_CONTEXT.md) |
| Documentation ownership and consolidation | [Documentation consolidation plan](plans/oneflow-documentation-consolidation-plan.md) |
| Change preflight and Change Gate | [Change preflight](source/change_preflight.rst) |

The top-level index answers what the system contains, why it is organized this
way, how its main flows work, and where the risks are. Topic documents explain
the algorithms, contracts, and implementation details.
