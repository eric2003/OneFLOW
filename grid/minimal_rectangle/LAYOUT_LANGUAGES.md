# Progressive Mesh Description

The rectangle case starts from a high-level description and compiles to the
existing OneFLOW atomic layout language. It does not replace the runtime grid
generator or change the established atomic parser.

## Three levels of description

1. **Intent level**: a human-readable statement such as "a unit square with
   four boundary edges, three points along each edge, and one boundary region."
2. **Topology level**: \`rectangle.json\` declares four ordered corners, four
   boundary edges implied by their order, and an independent point count for
   each edge.
3. **Atomic level**: \`compile_rectangle.py\` emits \`Point\`, \`Line\`, \`Dim\`,
   \`Ds\`, \`Boundary\`, \`Add L2F\`, and \`Add F2B\` statements for the current
   OneFLOW layout reader.

The high-level JSON is the editable input for this example. The generated
atomic file is an intermediate artifact; do not hand-edit both descriptions
independently. Recompile after changing the corners or \`edgePointCounts\`.

## Compile and inspect

From the repository root:

\`\`\`sh
python tools/mesh_layout/compile_rectangle.py \\
    grid/minimal_rectangle/rectangle.json \\
    --output /tmp/minimalRectangleLayout.txt \\\n    --svg-output /tmp/minimalRectangleTopology.svg
\`\`\`

Use the generated file with the existing \`GridMachine\` workflow. For this
rectangle, \`edgePointCounts: [3, 3, 3, 3]\` retains the current nine-node,
four-cell structured grid. Each edge count can be changed independently; the
current structured rectangle reader/generator determines how these counts
interact at corners and interior cells.

## Incremental direction

The levels should evolve from one shared geometric/topological model rather
than from independent hand-maintained files:

- add a natural-language description as documentation first, then consider a
  constrained command vocabulary rather than unrestricted language parsing;
- add a topology model for blocks, oriented edges, boundary regions, and
  per-edge spacing;
- compile that model to OneFLOW atomic syntax and validate it with the existing
  generator;
- visualize the same model as an SVG topology sketch, labeling points, edge
  IDs, directions, point counts, and boundary groups;
- later add an editor where changing a point count updates the model,
  generated atomic text, and preview together.

The SVG preview is generated from the same JSON model as the atomic layout, so changing a point count updates both outputs on the next compile. This first preview is intentionally read-only; a linked interactive editor can be added after the shared topology model is generalized.\n\nThe next generalization should extract the rectangle model from this example
into a format-neutral topology representation. Geometry/topology extraction
and JSON serialization should remain separate from this layout compiler.
