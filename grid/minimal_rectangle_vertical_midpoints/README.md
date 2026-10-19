# Minimal Rectangle with Vertical Midpoints

This case is a variant of `minimal_rectangle`. It adds one midpoint to each
vertical side while preserving the rectangle's four topological boundary edges.

The two vertical edges are each represented by two atomic line segments:
- right edge: lines 2 and 3, sharing point 5;
- left edge: lines 5 and 6, sharing point 6.

The layout still has one face with four side positions. Both line segments on
each vertical side map to the same face-side position, so six atomic line
segments do not mean six topological edges.

## Run

Copy this case to a writable workspace and run it through the standard
workflow:

```sh
python tools/mesh_workflow.py run \
    --case grid/minimal_rectangle_vertical_midpoints/case.json \
    --workdir /path/to/workspace/grid/minimal_rectangle_vertical_midpoints \
    --oneflow /path/to/build/bin/OneFLOW
```

The expected structured grid remains 3 by 3 nodes (9 nodes and 4 quadrilateral
cells). The midpoint points are shared endpoints of adjacent line segments,
not additional grid nodes beyond the original three-point vertical sides.
All generated artifacts stay in the external workspace under `results/`.
