# AI-Assisted Engineering Task Specification

This document defines the standard task contract for substantial AI-assisted or human-led
engineering work in OneFLOW. It is intended to make requests actionable, reviewable, and
reproducible. Use it together with `doc/source/architecture_review.rst` and the maintained
subsystem roadmap.

## 1. Project and execution context

Record the context before changing files:

| Field | Required value |
|---|---|
| Repository | `eric2003/OneFLOW` |
| Target branch | `hexin`, unless explicitly overridden |
| Base revision | Latest remote `hexin` HEAD, recorded as a full commit SHA |
| Local validation baseline | User-reported test result, e.g. `977/977`; identify it as pre-change evidence |
| Working environment | OS, compiler/toolchain, build configuration, and relevant optional backends |
| Scope | Runtime, mesh, numerical methods, I/O, Workbench, tests, documentation, or a named cross-cutting concern |

Do not infer the current branch revision from an earlier conversation. Fetch and inspect the
current remote HEAD. Do not target `master`, and do not assume that a user-reported test count
describes the post-change tree.

## 2. Objective and expected outcome

State one primary objective as an observable engineering outcome.

A valid objective identifies at least one of the following:

- a behavior or invariant that must be established;
- a dependency, ownership, or lifecycle boundary that must be clarified;
- a numerical method whose mathematical and implementation contracts must be verified;
- a reproducibility, diagnostic, or documentation capability that must be delivered.

Avoid objectives phrased only as activities, such as "continue refactoring", "clean up code",
or "add more tests". Explain what becomes safer, clearer, more reusable, or more verifiable.

## 3. Whole-system prioritization

Select one primary work item. Do not distribute effort uniformly across unrelated modules.

Rank candidate work using these criteria:

1. **System impact:** number and importance of workflows or subsystems affected.
2. **Risk reduction:** severity and likelihood of defects or invalid state prevented.
3. **Architectural leverage:** improvement to ownership, dependency direction, reuse, or testability.
4. **Scientific impact:** effect on conservation, stability, accuracy, convergence, or model validity.
5. **Evidence quality:** whether the change can be validated with representative tests and reproducible runs.
6. **Cost and disruption:** implementation effort, compatibility risk, and migration burden.

Prefer changes that reduce a shared architectural risk or establish a reusable contract. Defer
isolated renames, cosmetic restructuring, and low-impact helper extraction unless they are
necessary to complete a higher-value boundary.

The current project-wide priority order is maintained in
`doc/source/architecture_review.rst` and the relevant subsystem roadmap. Reassess it when
new source evidence or test results materially change the risk profile.

## 4. Required investigation before implementation

Before editing, establish the following from the target revision:

- actual entry points, callers, consumers, and relevant build targets;
- data representation, ownership, mutability, lifetime, and cleanup responsibilities;
- process-, case-, task-, and object-scoped state;
- normal, repeated-use, and failure execution paths;
- existing tests, regression cases, and known CI limitations;
- current documentation and any existing design decision that constrains the change.

Separate findings into **verified facts**, **unverified hypotheses**, and **proposed design**.
Do not convert a roadmap statement or a previous AI response into evidence without checking
the current source.

## 5. Implementation principles

- Preserve established behavior unless a behavior change is an explicit objective with tests.
- Introduce abstractions only when they establish a real contract, remove an undesirable
  dependency, clarify ownership, or enable meaningful verification.
- Prefer small, cohesive changes with explicit inputs, outputs, side effects, and failure
  semantics.
- Keep parsing, type validation, domain interpretation, state application, and execution
  distinct where the domain requires those separate responsibilities.
- Make resource ownership and teardown order explicit. Case teardown must not finalize
  process-scoped resources.
- Preserve compatibility through deliberate adapters when migration cannot be atomic.
- Improve names and comments when they clarify the responsibility being changed; do not
  make renaming a substitute for architectural progress.
- Keep C++ source files ASCII-only.
- Do not introduce a parallel implementation of production algorithms in Python or GUI code
  merely to make a test or demonstration easier.

## 6. Verification requirements

Tests must correspond to the claim being made.

| Evidence level | What it establishes |
|---|---|
| Unit / contract test | A local invariant, input rule, output rule, or failure condition |
| Component test | Collaboration between adjacent production components |
| Integration test | The real production path, including relevant state changes and side effects |
| Regression test | Compatibility of representative application behavior |
| Numerical verification | Mathematical properties such as consistency, conservation, stability, and convergence |
| Platform CI | Build and test behavior on the named OS, compiler, and configuration |

For each test, state the contract it verifies and what it does **not** verify. A passing
unrelated test suite does not establish correctness of the changed path. A test count is not
a quality metric without meaningful assertions.

For numerical algorithms, use appropriate combinations of analytic solutions, manufactured
solutions, limiting cases, conservation checks, grid/time refinement, and benchmark cases.
Record error norms, refinement ratios, tolerances, and assumptions where relevant.

Never claim a build, test, benchmark, or CI job passed unless its result was observed for the
exact commit being reported. Distinguish user-reported local results from tool-observed CI
results.

## 7. Documentation and scientific provenance

Documentation is a maintained engineering artifact, not a final packaging task.

### Source comments and generated API reference

Use concise, structured comments for stable public interfaces and non-obvious invariants.
Document intent, parameter meaning and units, ownership/lifetime, preconditions,
postconditions, side effects, failure behavior, and compatibility constraints where relevant.
Use stable cross-reference identifiers such as `@doc-algorithm`, `@doc-stage`, or
`@doc-theory` only when they resolve to maintained documentation.

Evaluate Doxygen as the source for API reference generation. A Doxygen/Sphinx integration
(such as Breathe) should begin as a small pilot on a stable module. Generated API reference
must complement, not replace, curated pages in `doc/source/`.

### Curated Sphinx documentation

Keep the following material in maintained documentation:

- architecture rationale, module responsibilities, dependency and lifecycle diagrams;
- mathematical models, derivations, assumptions, discretization, and boundary treatment;
- algorithm flowcharts, data-flow diagrams, geometric and mesh schematics;
- implementation mapping, limitations, verification strategy, and known failure modes;
- references to papers, textbooks, standards, and authoritative software manuals.

Use Sphinx mathematics for equations. Define symbols, units, sign conventions, assumptions,
and the conditions under which a formula applies. Use Mermaid or Graphviz for process and
dependency diagrams, and SVG or generated figures for geometric and numerical illustrations.
Figures must have captions and explain the information the reader should extract.

### References and external software

Cite original papers, established textbooks, standards, and authoritative manuals when they
support a mathematical or engineering claim. A citation records provenance; it does not
prove that OneFLOW implements the cited method correctly.

Other open-source and commercial CFD/CAE manuals may be used for terminology, workflow
comparison, and publicly documented methods. Distinguish published theory, vendor claims,
independent analysis, and OneFLOW's verified implementation. Do not imply access to or
knowledge of proprietary internal algorithms.

## 8. Documentation-to-code traceability

For important algorithms and architectural contracts, maintain a navigable chain:

`theory / requirement -> mathematical or behavioral contract -> design page
-> source API / implementation -> tests -> representative verification evidence`.

Use stable identifiers where they provide real value. Avoid introducing a large metadata
system before a few end-to-end examples demonstrate the need. At milestones, check for
orphaned documentation links, undocumented public contracts, and tests whose names or
assertions no longer match implementation behavior.

## 9. Scope control and decision records

For each change, explicitly state:

- what is included;
- what is intentionally excluded;
- compatibility constraints;
- migration or follow-up work;
- alternatives considered when the decision has durable architectural consequences.

Create an Architecture Decision Record only for decisions with meaningful alternatives,
long-term consequences, or cross-subsystem impact. Do not create an ADR for every local edit.

## 10. Completion criteria and delivery record

A task is complete only when the applicable criteria below are satisfied:

- implementation matches the stated objective and preserves required compatibility;
- tests exercise the changed contract at the appropriate evidence levels;
- relevant documentation and diagrams are updated;
- source comments and generated API references do not contradict curated documentation;
- build, test, numerical, and CI evidence is recorded against exact revisions;
- remaining limitations and follow-up work are explicit;
- the whole-system priority assessment is revisited.

The delivery summary must include:

1. base and resulting commit SHAs;
2. problem statement and system-level rationale;
3. changed files and responsibility changes;
4. tests added or modified and their exact contracts;
5. commands and observed results, including failures;
6. exact-commit CI links and status;
7. documentation and reference updates;
8. unresolved risks and the next highest-value action.

Do not report a task as complete when a required integration path or platform check remains
unverified. Use **Complete**, **Partial**, **Blocked**, or **Deferred** with a short evidence-based
explanation.

## 11. Reusable task request

Copy and complete this block when starting a substantial task:

```text
REPOSITORY: eric2003/OneFLOW
BRANCH: hexin
BASE: latest remote HEAD; record full SHA
LOCAL_BASELINE: <user-reported result and environment>

PRIMARY_OBJECTIVE:
<observable engineering outcome>

SYSTEM_CONTEXT:
<affected workflow, subsystem, upstream/downstream dependencies>

PRIORITY_RATIONALE:
<system impact, risk reduction, architectural/scientific leverage>

REQUIRED_INVESTIGATION:
<entry points, ownership, lifecycle, failure paths, build/test targets>

CONSTRAINTS:
- Preserve compatibility unless explicitly changed and tested.
- Keep C++ source ASCII-only.
- Do not introduce abstractions without a concrete contract and consumer.
- Do not claim unobserved build/test/CI results.
- Follow the current hexin branch and repository contribution conventions.

IMPLEMENTATION_SCOPE:
<planned code changes and explicit exclusions>

VERIFICATION:
<unit/component/integration/regression/numerical/platform evidence>

DOCUMENTATION:
<doc/source pages, equations, diagrams, references, API comments>

ACCEPTANCE_CRITERIA:
<objective, testable completion conditions>

DELIVERY:
<commit SHAs, changed files, test evidence, CI links, limitations, next priority>
```

For current project priorities and subsystem status, consult
`doc/source/architecture_review.rst`, the relevant roadmap, and `doc/AI_CONTEXT.md`.
This template defines how to express and execute a task; it does not replace source inspection.
