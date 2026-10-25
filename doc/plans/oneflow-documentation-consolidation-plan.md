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
| `docs/MESH_DATA_AND_VISUALIZATION.md` | Mesh data, reconstruction, and visualization target architecture | Rich design note with status and test strategy | Candidate source for a consolidated mesh architecture page; compare against similarly named root file before choosing canonical content |
| `docs/architecture/grid-data-contracts-and-visualization.md` | Mesh summary/schema and adapter contracts | Related to the mesh architecture note, but may contain more concrete schema decisions | Compare and merge unique details into one canonical mesh design; avoid parallel authoritative descriptions |
| `MESH_DATA_AND_VISUALIZATION.md` (root) | Earlier/alternate mesh architecture note | Overlaps with `docs/MESH_DATA_AND_VISUALIZATION.md`; headings and information levels differ | Compare claims and history; then keep one canonical document and replace other copies with links or remove only after links are checked |
| `GRID_TOPOLOGY.md` (root) | Mesh topology semantics, invariants, cache refresh contract | More specific than the general mesh interchange document | Incorporate into the mesh design family as a distinct topology chapter/page; do not merge away unique invariants |
| `docs/ONEFLOW_WORKBENCH_ARCHITECTURE.zh-CN.md` | Workbench architecture memo in Chinese | Useful architecture memo, currently outside canonical root | Review against implementation; move/copy into the curated design documentation only after source and link strategy is settled |
| `docs/ONEFLOW_WORKBENCH_ARCHITECTURE.en.md` | Workbench architecture memo in English | Corresponding translation, but section detail is not perfectly symmetric with Chinese | Establish one conceptually authoritative version and keep translations aligned; don't silently assume the two files are identical |
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

## 4. Mesh documentation comparison: first consolidation pass

The first content review found three overlapping mesh design notes with distinct value:

- `docs/MESH_DATA_AND_VISUALIZATION.md` provides the broad reconstruction and visualization goals and an information-level model.
- `docs/architecture/grid-data-contracts-and-visualization.md` records concrete L0 summary and JSON-adapter boundaries, including the distinction between `GridSummary` and runtime mesh-info output.
- Root-level `GRID_TOPOLOGY.md` contains topology-specific invariants, boundary ghost-cell normalization, derived-adjacency rebuild behavior, and mutation/cache limitations.
- Root-level `MESH_DATA_AND_VISUALIZATION.md` overlaps the broad design note and remains to be compared for unique implementation details before it can be retired.

A curated Sphinx overview now provides a consistent information vocabulary and points readers toward the architecture goals, snapshot/adapter boundary, topology invariants, and validation strategy. This does not delete or relocate the source notes. In particular, the different level numbering in earlier drafts must not be treated as an established serialized schema; the schema remains versioned design work until implemented and tested.

## 5. Recommended review order

1. Finish the inventory and check internal links and Sphinx toctrees.
2. Compare the three mesh-related documents semantically and identify one canonical design page plus any separate topology/reference pages.
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


## 8. Mesh documentation: topology page published in source tree

The first mesh-documentation pass now has two curated Sphinx pages:

- `doc/source/mesh_architecture.rst` provides the information-level model, canonical-versus-derived data rules, snapshot/adapter boundary, and validation strategy.
- `doc/source/grid_topology.rst` preserves the unique topology contracts from the root engineering note: 1D/2D/3D meaning, `UnsGrid::NormalizeBc()` ghost-ID convention, explicit `c2f`/`c2c` rebuild responsibilities, face-array reorder alignment, and the known limitation caused by publicly mutable source arrays.

The root `GRID_TOPOLOGY.md` remains in place as engineering source material. It should only be retired or replaced with a pointer after its remaining unique content and inbound links have been checked. The two mesh overview notes still need a final line-by-line comparison before any duplicate is removed.

The Sphinx home page description was also rewritten to emphasize the current repository shape and long-term workflow direction rather than making broad capability claims without revision-specific evidence. External article titles in existing link lists are not changed.
