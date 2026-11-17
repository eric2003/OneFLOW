# Two Blocks Sharing an Edge with Endpoints Only

This is the minimum two-block 2D layout case.

- Two adjacent unit squares.
- One quadrilateral cell per block.
- The shared edge has exactly two nodes: its endpoints only.
- No interior node is inserted on the shared edge.

Expected combined topology: 6 unique nodes, 2 quadrilateral cells, 7 unique
edges, 6 boundary edges, and 1 internal edge.

**This case only verifies the basic shared-edge connection between two blocks.**
It does not test subdivision or interior interface nodes. Those are covered by
`two_blocks_shared_edge_one_interior_node`.

## Run

Copy this directory to a writable workspace outside the repository and run:

```text
OneFLOW.exe r <case-directory>
```

Validate its case profile with:

```sh
python tools/mesh_cases.py validate two_blocks_shared_edge_endpoints_only
python tools/test_mesh_cases.py
```
