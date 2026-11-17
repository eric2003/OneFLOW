# Two Blocks Sharing an Edge with One Interior Node

This is the second two-block 2D layout example. It keeps the same two adjacent
unit squares as the endpoint-only example, but adds one interior node to the
shared edge.

- Block 1: `[0, 1] x [0, 1]`
- Block 2: `[1, 2] x [0, 1]`
- Shared edge: `x = 1`, represented by two endpoints and one interior node.
- Each block has a 3 x 3 node lattice and four quadrilateral cells.

Expected combined topology: 15 unique nodes, 8 quadrilateral cells, 21 unique
edges, 12 boundary edges, and 9 internal edges.

**This case specifically verifies matching a shared edge that contains an
interior discretization node.** It is not the minimum two-block connection;
see `two_blocks_shared_edge_endpoints_only` for the two-endpoint baseline.

## Run

Copy this directory to a writable workspace outside the repository and run:

```text
OneFLOW.exe r <case-directory>
```

The case profile and layout are checked with:

```sh
python tools/mesh_cases.py validate two_blocks_shared_edge_one_interior_node
python tools/test_mesh_cases.py
```
