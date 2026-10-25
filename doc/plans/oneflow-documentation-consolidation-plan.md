# OneFLOW Documentation Consolidation Plan

Status: In progress
Last updated: 2026-10-26
Scope: Documentation assets on the `hexin` branch. This plan does not authorize deleting or moving files until each content group has been reviewed.

## 1. Decisions

1. Keep `doc/` as the single canonical documentation root.
2. Keep Sphinx and Read the Docs as the publishing toolchain for now. Do not migrate to MkDocs solely to reorganize content.
3. Separate the repository's working documents from the curated Sphinx site:
   - `doc/architecture/`, `doc/plans/`, `doc/runbooks/`, `doc/evidence/`, and `doc/reports/` hold engineering records by purpose.
   - `doc/source/` holds the curated, navigable documentation site.
4. Treat root-level Markdown files and the top-level `docs/` directory as legacy/discovery locations to review, not as destinations for new documentation.
5. Use Markdown for new engineering notes where practical. Keep existing reStructuredText pages working; evaluate MyST in a small isolated trial before changing the Sphinx source format.
6. Do not mechanically move or rename existing files. Preserve useful history and links; consolidate only after content comparison.

## 2. Current inventory and proposed disposition

| Current path | Subject / role | Initial assessment | Proposed action |
|---|---|---|---|
| `README.md`, `README_zh_CN.md` | Project entry points | Keep at repository root | Keep concise; link to published docs and contributor/developer guidance |
| `doc/README.md` | Documentation policy and directory entry | Canonical policy already exists | Update with this plan and clarify canonical locations |
| `doc/DOCUMENT_NAMING.md` | Engineering-document naming rules | Canonical naming policy | Extend only if a new document category is approved; do not create competing naming rules |
| `doc/source/` | Existing Sphinx site with `.rst` pages and a large topic reference | Existing publication source; not yet a modern OneFLOW user/developer handbook throughout | Preserve current pages and toctrees; improve navigation incrementally |
| `doc/UserManual-chinese.doc`, `doc/DevelopManual-chinese.doc` | Historical Chinese manuals | Useful legacy assets, but freshness and editability need review | Retain until contents are reviewed; later mark as historical or extract still-valid material |
| `doc/architecture/oneflow-1d-euler-dcu-architecture-audit-20260904.md` | Dated implementation audit | Correct fit for engineering architecture records; includes a historical branch baseline | Keep as a dated audit; do not present it as the current architecture without a fresh review |
| `doc/evidence/*` | Build and runtime evidence | Fits existing evidence category | Keep; link from relevant reports/design pages where useful |
| `doc/reports/*` | Performance, regression, delivery, and project reports | Mostly fits report taxonomy; some architecture and handoff records are nested under reports/architecture | Review directory semantics later; avoid bulk moves in the first pass |
| `doc/runbooks/*` | Repeatable operational procedures | Fits runbooks | Keep; link stable procedures from the developer guide |
| `docs/MESH_DATA_AND_VISUALIZATION.md` and root `MESH_DATA_AND_VISUALIZATION.md` | Earlier mesh reconstruction/visualization design notes | Compared; overlapping goals and incompatible draft information-level numbering | Both now point to the curated Sphinx pages; retain as legacy entry points, not parallel authorities |
| `docs/architecture/grid-data-contracts-and-visualization.md` | Mesh summary/schema and adapter contracts | Contains implementation-specific detail absent from broad design notes | Retain as the detailed engineering record and link from the curated mesh overview |
| `GRID_TOPOLOGY.md` (root) | Mesh topology semantics, invariants, cache refresh contract | More specific than the general mesh interchange document | Keep as source material until remaining unique content and inbound links are checked |
| `docs/ONEFLOW_WORKBENCH_ARCHITECTURE.zh-CN.md` | Workbench architecture memo in Chinese | Useful architecture memo, currently outside canonical root | Review against implementation; keep as source material during alignment |
| `docs/ONEFLOW_WORKBENCH_ARCHITECTURE.en.md` | Workbench architecture memo in English | Corresponding translation, but section detail is not perfectly symmetric with Chinese | Keep as source material during alignment; do not assume the translations are identical |
| `AGENTS.md`, `CLAUDE.md` | Repository/agent instructions | Operational instructions, not user documentation | Keep at root; link to relevant developer conventions where appropriate without duplicating their full contents |

## 3. Target information architecture

This is a target for gradual convergence, not a mandate to create every directory immediately.

### Engineering records

- `doc/architecture/`: stable design notes and dated architecture audits.
- `doc/plans/`: living development plans and documentation consolidation tasks.
- `doc/runbooks/`: repeatable build, test, debug, and operations procedures.
- `doc/evidence/`: verifiable test/build/runtime evidence.
- `doc/reports/`: published project, performance, regression, and delivery reports.
- `doc/handoff/`: only if recurring handoff documents justify this category; avoid creating empty categories just to match a template.
- `doc/_drafts/`: temporary drafts only; promote or discard them deliberately.

### Curated Sphinx site

Organize pages by reader intent, while retaining useful existing topic references:

- Getting started: install, configure, build, and run a first case.
- User guide: cases, project/workspace concepts, mesh workflows, and result inspection.
- Developer guide: code layout, build configurations, tests, debugging, and contribution workflow.
- Design: runtime architecture, mesh geometry/topology/boundaries, workflow execution, and workbench boundaries.
- Reference: parameter files, data schemas, command-line arguments, file formats, and stable interface contracts.
- Existing scientific-computing reference: preserve useful CFD, numerical, language, and tool pages; reorganize only when navigation benefits are clear.

Avoid duplicate copies between engineering records and the published site. A detailed engineering note can remain in its working location while a curated design page summarizes the stable conclusions and links to the note.

## 4. Mesh documentation comparison and consolidation

The comparison identified distinct content that should remain accessible without competing as canonical design:

- `doc/source/mesh_architecture.rst` is the curated overview: information levels, canonical versus derived data, snapshot/adapter boundary, and validation strategy.
- `doc/source/grid_topology.rst` preserves topology-specific invariants: 1D/2D/3D meaning, `UnsGrid::NormalizeBc()` ghost-ID convention, explicit `c2f`/`c2c` rebuild responsibilities, face-array reorder alignment, and the limitation caused by publicly mutable source arrays.
- `docs/architecture/grid-data-contracts-and-visualization.md` retains concrete implementation details for `GridSummary`, `oneflow.grid-summary`, the separate runtime `oneflow.mesh-info` output, and the proposed package layout.
- The root and `docs/` copies of `MESH_DATA_AND_VISUALIZATION.md` have been changed into legacy entry points linking to those authoritative pages. Their earlier taxonomies are explicitly described as drafts, not established serialized schemas.
- Root `GRID_TOPOLOGY.md` remains untouched until its remaining unique content and inbound links are checked.

Do not infer a serialized package schema from the historic L0-L4/L0-L5 labels. Any future interchange schema must be explicitly versioned and tested. A summary is metadata, not a reconstructable mesh package.

## 5. Recommended review order

1. Finish the inventory and check internal links and Sphinx toctrees.
2. Check inbound links to root `GRID_TOPOLOGY.md` and compare any remaining unique material before replacing or retiring it.
3. Review the Chinese and English Workbench architecture memos for translation parity and implementation status.
4. Review the two legacy Word manuals and decide whether they are current, archival, or sources for extracted material.
5. Improve `doc/source/index.rst` navigation so new users, users, and developers have clear entry points without hiding the existing reference library.
6. Evaluate a minimal MyST proof of concept in a separate page; adopt it only if cross-references, navigation, and Read the Docs builds work cleanly.
7. Add lightweight documentation checks (broken links and Sphinx build) to CI when the current build path is understood and stable.

## 6. Rules for ongoing development

- A code change that changes a module boundary, data contract, user-visible workflow, or developer procedure should update the relevant document or explain why no document change is needed.
- Important architectural decisions should record context, decision, alternatives, consequences, and verification evidence. Keep these records concise.
- Label planned architecture as proposed; describe current behavior only when supported by code/tests.
- Tests and evidence should support claims, not be replaced by prose.
- Keep code and documentation changes on `hexin`; do not target `master`.
- Do not claim a documentation build or test passed unless it was run against the changed revision.

## 7. Project-description style

For repository-authored README files, manuals, architecture notes, and new documentation, describe OneFLOW through its technical purpose, implemented capabilities, architecture, numerical methods, workflows, reproducibility, and verified limitations. Do not frame the project as "国产", "domestic", or by national origin. Existing external article links may retain their original titles verbatim; do not rewrite those titles solely to remove wording from the linked article. Distinguish implemented and tested behavior from roadmap or long-term vision.

## 8. Workbench documentation alignment

A curated page, `doc/source/workbench_architecture.rst`, consolidates stable concepts from the Chinese and English Workbench memos: project, case, workspace/run directory, workflow definition, workflow run, artifact, execution modes, persistence, and process lifecycle.

The source memos are not identical. The Chinese memo distinguishes cases and run modes more explicitly; the English memo gives more detail on the versioned workflow model and incremental implementation. The curated page combines these ideas and adds a clear status rule: planned background-task recovery, visualization, diagnostics, and related features are not claims of implementation until code and tests demonstrate them. The source memos remain available as engineering records during review.

The Sphinx `info.rst` page is a concise documentation entry point. Existing external article titles in `doc/source/user.rst` remain unchanged as requested.

## 9. Verification status

The documentation edits in this pass were committed to `hexin`. A Sphinx build, broken-link check, and CI run have not been performed in this session; do not treat the documentation as build-verified until those checks run against the resulting branch head.
