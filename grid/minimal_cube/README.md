# Minimal Cube Case

This case describes one unit cube using OneFLOW's existing atomic grid-layout
syntax. It is intended as a small, editable case template, parallel to
`grid/minimal_rectangle`.

## Geometry and topology

- Points 1-4 are the lower square; points 5-8 are the upper square.
- Lines 1-8 are the two square perimeters.
- Lines 9-12 connect corresponding lower and upper corners.
- Each `Add L2F` group assigns four ordered edges to one of the six faces.
- `Add F2B` assigns the six faces to the six local positions of one block.
- `Dim = 3` gives three nodes along each edge, producing a nominal 3 x 3 x 3
  structured grid (8 hexahedral cells).

The `Add F2B` positions in the input are one-based. The generator converts
them to zero-based indices when selecting the block's six face slots.

## Run

Use `minimalCubeLayout.txt` as the layout input in the existing grid-layout
generation workflow, with the spatial dimension set to 3D. Configure the
workflow's normal source-grid and boundary-output paths as usual.

This is deliberately a direct layout sample rather than a unit-test fixture.
When the topology or generated coordinates need to be investigated, edit this
case and run it through the regular application workflow.
