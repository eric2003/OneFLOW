# Evolutionary scientific-computing architecture

Status: architectural proposal and design vocabulary; not a claim that OneFLOW already has a general plugin runtime.

## 1. The objective is evolution, not a file-reader feature

A scientific code must preserve more than the ability to execute today's input. It must preserve the meaning of data, the assumptions of algorithms, the ownership of state, the evidence supporting numerical results, and the ability to replace one implementation without rewriting every consumer.

The long-term objective is an ecosystem in which algorithms and data transformations can be described, inspected, tested, composed, replaced, and compared independently. A control-file reader is a useful first specimen because it already contains a sequence of distinct transformations:

```text
legacy text
  -> lexical/syntax parsing
  -> configuration document
  -> typed database entries
  -> runtime parameter state
  -> distributed state
  -> diagnostic snapshot
```

This is not merely a chain of functions. Each arrow changes representation and therefore needs a contract: accepted input, produced output, invariants, side effects, failure behavior, and verification evidence.

## 2. The DNA/RNA factory analogy, made operational

The biological analogy is useful if each term has a precise engineering counterpart:

| Analogy | OneFLOW architectural meaning |
| --- | --- |
| DNA | Stable, versioned description of algorithm/data contracts, assumptions, units, and provenance |
| Transcription | Convert a human-authored or legacy representation into a structured intermediate representation |
| RNA processing | Normalize, validate, enrich, or adapt the intermediate representation without changing its declared meaning |
| Ribosome / assembly | Compose compatible units into an executable pipeline or solver workflow |
| Protein | A concrete executable algorithm, data structure, mesh, field, or solver result |
| Quality control | Preconditions, postconditions, invariants, tests, numerical verification, and provenance checks |
| Mutation / evolution | A reviewed change to an algorithm, schema, adapter, or implementation, with compatibility and evidence recorded |

A “small algorithm robot” is not necessarily a class, a plugin, or a dynamically loaded library. It is first a unit of responsibility with an explicit contract. A free function can qualify if it has clear inputs, outputs, state effects, and tests. Do not introduce a factory or plugin manager merely to make the analogy look complete.

## 3. The computation-unit contract

The proposed machine-readable descriptor is defined in `doc/specs/computation-unit.schema.json`. It captures:

- identity and version;
- role in a transformation pipeline;
- typed input and output ports, including units and ownership;
- preconditions and postconditions;
- side effects such as filesystem writes, global-state mutation, collectives, device state, and logging;
- failure and rollback policy;
- determinism and lifecycle scope;
- compatibility adapters;
- verification evidence and known limits;
- scientific provenance and references.

The schema is a design contract, not an executable runtime ABI. Its first purpose is to make architecture review precise and machine-checkable. Later tooling may validate descriptors, build dependency graphs, generate test matrices, or render pipeline diagrams. Such tooling must not silently infer missing scientific assumptions.

## 4. Intermediate representations are the key to flexibility

The system should avoid direct coupling between every input format and every consumer. Prefer a small number of explicit intermediate representations (IRs), each with a documented semantic scope.

Examples of candidate IRs:

1. `ConfigDocument`: parsed configuration syntax before it is committed to runtime storage.
2. Typed parameter snapshot: validated values with declared types and defaults.
3. Mesh snapshot: canonical geometry/topology values independent of runtime solver globals.
4. Topology and boundary views: derived connectivity with explicit invalidation/rebuild rules.
5. Algorithm/workflow description: selected algorithms, dependencies, data requirements, and execution scope.
6. Result/diagnostic snapshot: outputs with units, provenance, and verification context.

An IR is valuable when multiple producers or consumers need to interoperate, or when it freezes a meaningful semantic boundary. It is not valuable when it simply duplicates a data structure or adds conversions with no decoupling benefit.

Adapters should be small and directional. A legacy input adapter may translate old syntax into a current IR; a serializer may translate an IR into JSON; neither should own solver lifecycle or numerical policy. This makes it possible to replace a parser, storage backend, or output format without rewriting numerical algorithms.

## 5. Architecture principles

### 5.1 Separate representation from meaning

Parsing answers “what was written?” Validation answers “is it well-formed and type-correct?” Domain interpretation answers “what does this value mean for this case?” Application answers “which runtime state is changed?” Execution answers “what computation is performed?”

These stages may share optimized implementation details, but their contracts should remain distinguishable.

### 5.2 Make state and side effects explicit

A component that reads a file, mutates a global database, performs a parallel collective, or changes device state is not a pure transformation. Declare the effect and its lifecycle scope. Prefer value snapshots at boundaries. Where legacy globals remain necessary, isolate them behind adapters and explicit setup/teardown contracts.

### 5.3 Preserve semantics while allowing implementation replacement

Keep stable external contracts while moving implementation behind module boundaries. Version schemas and IRs when meaning changes. Use compatibility adapters for old representations; do not scatter version checks throughout solver algorithms.

### 5.4 Make failure behavior composable

Every stage must define whether it leaves state unchanged, rolls back its own mutations, or leaves diagnostics. Cleanup ownership must be unambiguous. A case-scoped operation must not finalize process-scoped resources. A failed stage must not silently contaminate the next case.

### 5.5 Tests are executable architecture constraints

Tests do more than detect defects:

- unit tests pin a local invariant;
- component tests pin collaboration between adjacent transformations;
- integration tests pin the real production pipeline and its side effects;
- regression tests pin scientific behavior against known cases;
- manufactured-solution, conservation, and convergence tests provide evidence about numerical correctness.

If a tiny contract test requires linking the entire solver executable, that is evidence to inspect the dependency boundary. It is not automatic proof that another abstraction is required. A green test only establishes the property it actually asserts.

### 5.6 Scientific claims require provenance and verification

For numerical methods, retain a traceable chain:

```text
published theory/reference
  -> OneFLOW mathematical formulation
  -> implementation and configuration
  -> focused property tests
  -> manufactured/benchmark case
  -> regression evidence and known limits
```

Citations establish provenance, not correctness. Unit tests establish local properties, not universal model validity. Numerical baselines must be reproducible and reviewed before being changed.

## 6. Current specimen: control-file ingestion

The current production flow is:

```text
script/control.txt
  -> ControlFileReader: resolve ordered script list
  -> ConfigLoader: parse legacy text
  -> ConfigDocument: structured representation
  -> ConfigDatabaseAdapter: commit typed entries
  -> DataPara: runtime parameter storage
  -> ControlDatabaseIO: broadcast and diagnostic dump
```

The first extraction moved reader orchestration into `ControlFileReader.cpp` and database broadcast/dump operations into `ControlDatabaseIO.cpp`. Public declarations in `ParaFile.h` remain stable. The legacy parser helpers remain in `ParaFile.cpp`.

This is deliberately a source-level boundary first. It does not yet establish that the reader is independently linkable in a focused integration target. In particular, `ReadMultiScriptFiles()` still calls `ReadOneFLOWScriptFile()`, which is currently implemented with legacy/configuration helpers in `ParaFile.cpp`. The next dependency review must determine whether the script-parse adapter should move to a separate translation unit or whether a narrower parser library target already exists.

## 7. Expected problems and how to respond

| Expected problem | Why it occurs | Preferred response |
| --- | --- | --- |
| New source file is not compiled | Source collection may be explicit or cached | Inspect CMake source discovery and build graph; add explicit source registration if required |
| Link test still pulls legacy parser | Reader calls a parser entry point in the broad legacy TU | Extract the parser adapter only after mapping its dependencies; do not duplicate it in tests |
| Global state contaminates a second case | Database, paths, logs, parallel or solver registries outlive one case | Test A then B in one process; assert B's defaults and teardown contract |
| Parallel reader hangs or reads no files | Reader assumes server rank and initialized process metadata | Initialize the real serial backend; test collective boundaries separately |
| Error cleanup hides the original failure | Teardown throws while handling setup failure | Preserve the first exception; assert cleanup and exception identity/message |
| A generic plugin framework grows too early | Architecture is designed before stable contracts exist | Start with ordinary C++ functions and value objects; add runtime registration only when a real extension use case requires it |
| IRs multiply without value | Every layer creates a duplicate representation | Require multiple consumers/producers or an independently testable semantic boundary |
| Machine-readable metadata drifts | Documentation descriptors are not checked | Add schema validation and source/test links gradually; keep unverified claims explicit |
| Numerical evidence is mistaken for software tests | A passing unit test is treated as model validation | Separate unit, integration, regression, manufactured-solution, conservation, and convergence evidence |

## 8. Evolution roadmap

### Short term: make today's boundaries real

1. Complete the control-reader source split without changing parameter semantics.
2. Confirm automatic source discovery and compile the production executable on Windows and Linux.
3. Extract only the parser adapter dependencies needed to link the real reader in a focused integration test.
4. Build a temporary multi-file fixture; initialize the real serial parallel runtime; assert typed values and `log/database.log`.
5. Add sequential-case and failure-cleanup tests at the correct lifecycle ownership boundary.
6. Keep `doc/AI_CONTEXT.md`, the runtime architecture page, and this investigation synchronized with exact evidence.

### Medium term: establish reusable data contracts

1. Stabilize `ConfigDocument` and typed parameter snapshots as versioned representations.
2. Apply the same pattern to mesh snapshots, topology views, boundary data, and result diagnostics.
3. Define explicit ownership, units, invalidation, and lifecycle rules for each representation.
4. Add schema validation, generated dependency diagrams, and test-to-contract traceability.
5. Introduce compatibility adapters where legacy formats or APIs must coexist with new representations.

### Long term: an algorithm and dataflow ecosystem

1. Describe solver pipelines as compositions of compatible computation units.
2. Allow alternative implementations to satisfy the same mathematical and data contracts.
3. Support CPU, accelerator, serial, and distributed implementations through declared capabilities and effects.
4. Generate test plans and provenance reports from unit descriptors and pipeline graphs.
5. Compare algorithm variants using reproducible verification suites and controlled numerical baselines.
6. Consider dynamic plugin loading only when external extension, runtime selection, or independent release cycles justify the complexity.

The target is not a universal abstraction layer over every line of code. It is a small set of stable semantic contracts that lets the system change shape without losing correctness, traceability, or scientific meaning.

## 9. Review questions for every future refactor

Before introducing a new unit or IR, answer:

1. Which concrete producer and consumer are being decoupled?
2. What representation crosses the boundary, and what does it mean?
3. Which invariants, units, ownership, and lifecycle rules apply?
4. What side effects occur, and who owns rollback?
5. Which test proves the contract, and which higher-level evidence remains missing?
6. Can an alternative implementation be substituted without changing consumers?
7. Is this a real reusable boundary, or merely another wrapper?

If these questions have no concrete answers, keep the implementation simple and defer the abstraction.
