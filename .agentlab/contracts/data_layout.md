# Data Layout Contract (must be filled by a human; agents must not guess)

- Field layout (AoS / SoA / other): TODO
- Conserved variable order and count: TODO
- Index convention (i,j,k order, ghost included?): TODO
- Ghost width: TODO
- Mesh type (structured/unstructured) and partitioning: TODO
- Non-dimensionalization and units: TODO
- Memory ownership (who allocates/frees Field): TODO
- Precision (double/mixed): TODO
- Parallel model (MPI / OpenMP / GPU): TODO

## v0.1 recommended minimum (replace with your decision)

- Structured 1D, AoS, one conserved scalar, ghost width 1, `double`, serial `Comm` no-op.
- If you accept this, delete the TODO lines above and state the decision here.
