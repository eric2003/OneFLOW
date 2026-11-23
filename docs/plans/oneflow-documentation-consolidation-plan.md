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
- Design: whole-system architecture review and development principles, runtime architecture, mesh geometry/topology/boundaries, workflow execution, and Workbench boundaries.
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

1. Continue the inventory by checking internal links and Sphinx toctrees against an actual build; the landing page now has task-oriented entry points and retains the existing reference library.
2. Check inbound links to root `GRID_TOPOLOGY.md` and compare any remaining unique material before replacing or retiring it.
3. Review the two legacy Word manuals and decide whether they are current, archival, or sources for extracted material.
4. Evaluate a minimal MyST proof of concept in a separate page; adopt it only if cross-references, navigation, and Read the Docs builds work cleanly.
5. After the new `hexin` Sphinx workflow runs, address real build errors first; add a separate broken-link checker only after defining how external-link failures and legacy reference pages should be handled.

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

The Chinese and English source memos were reviewed and aligned around the same core concepts: Project, Case, Workspace, Run location, Run mode, Workflow, WorkflowRun, and Artifact. The Chinese memo's duplicate Workspace/运行空间 entry was corrected to distinguish the actual run location. Both language versions now distinguish independent multi-case execution from a dependency workflow and clarify that a multi-case run may mean separate process launches or one process receiving multiple directory arguments. The curated page remains the stable navigation target and explicitly states that planned background-task recovery, visualization, diagnostics, and related features are not claims of implementation until code and tests demonstrate them.

The new `doc/source/workbench_user_guide_zh.rst` is the Chinese user-manual entry for the optional GUI and is included in the Sphinx User Guides navigation. `tools/oneflow_gui/README.md` remains the detailed module-level operational reference; `doc/source/workbench_architecture.rst` records current implementation boundaries and gaps. Changes to visible GUI behavior should update the user guide and README, while service ownership or execution-lifecycle changes also update the architecture page.The Sphinx landing page `doc/source/index.rst` now presents task-oriented entry points for getting started, using OneFLOW, mesh/workbench workflows, and developer resources while retaining the existing reference library in hidden toctrees. `doc/source/info.rst` now links readers directly to the current `hexin` README, installation pages, grid and regression guides, and Workbench guide. Existing external article titles in `doc/source/user.rst` remain unchanged as requested. A dedicated `.github/workflows/docs-check.yml` now requests a Sphinx HTML build on `hexin` pushes and pull requests when documentation paths change; this is a build check, not a deployment workflow.

## 9. Verification status

The mesh-note comparison, legacy-pointer updates, bilingual Workbench memo alignment, Sphinx landing-page improvements, Workbench README formatting fix, and `hexin` documentation-check workflow were committed to `hexin`. The runtime architecture page and README navigation links have now been added. The documentation-check workflow completed successfully for commit `f3219d9a6024e8f8e70700f4e499363b728f3162`, which included the new page and its configuration details. The later README link-only updates are covered by the normal Windows/Linux CI, which is still running at the current branch head. A separate automated broken-link check has not been run.


## 10. System architecture and code-alignment pass

The next documentation pass broadens from asset consolidation to a code-backed system map.
The new curated page `doc/source/runtime_architecture.rst` documents:

- design intent and the boundary between the C++ runtime, case inputs, orchestration, and presentation;
- the current application entry point and the `Simulation -> SimuImp -> SimuContext -> TaskRegistry -> ISimuTask` execution path;
- the distinction between build configuration, case runtime parameters, and Workbench project/workflow/run data;
- patterns visible in current code, including task commands, registry/factory construction, context objects, parsing adapters, and explicit lifecycle management;
- the limits that remain in place, including legacy free-function task bodies and shared runtime state;
- a development loop that ties behavior contracts, focused tests, implementation changes, documentation updates, and validation evidence together.

This page is an initial source-backed map, not a substitute for module-specific contracts. It should evolve alongside real refactors. Keep the detailed mesh and Workbench pages authoritative for their respective domains rather than duplicating their content here.

For the optional Python GUI, `doc/source/workbench_user_guide_zh.rst` is the user-facing Sphinx guide, `tools/oneflow_gui/README.md` is the detailed module and operations reference, and `doc/source/workbench_architecture.rst` records service boundaries and implementation gaps. Update the user guide when visible behavior changes, the README when module workflows or commands change, and the architecture page when ownership or execution lifecycle changes.

### Next code-alignment order

1. Trace the control-file path from case setup through task selection and parameter lookup; document which component parses syntax, validates values, and applies runtime meaning.
2. Trace one representative grid-generation task and one solve task end to end, identifying global state, lifecycle ownership, and seams that can be tested without booting the full runtime.
3. Compare the documented Workbench run-location and workflow contracts with the actual execution services and tests.
4. Use those findings to choose one small code refactor with a behavior test, then update the architecture page only where the implementation contract changes.
5. Record gaps as explicit current limitations or proposed boundaries; do not document planned behavior as implemented.

The first architecture page has been added to the Sphinx navigation. The documentation-check workflow and the full Windows/Linux CI runs must still be checked for the resulting branch head; do not claim the page or build is validated until those checks complete.


## 11. Refactoring progress as a maintained engineering record

The runtime architecture page describes the system map and the rules for code changes. The
living plan `doc/plans/runtime-refactoring-roadmap.md` records the staged implementation
direction, acceptance criteria, and verified progress for runtime refactoring.

Keep this record outcome-oriented rather than a transcript of every thought or command:
capture the observed problem, current call path, intended contract, smallest useful change,
validation evidence, and remaining gap. Update a stage only when the source, tests, and
current CI evidence support the status. The first selected investigation is case I/O root
resolution in `SimuContext::SetupCaseEnvironment()`; it is not yet marked implemented.
