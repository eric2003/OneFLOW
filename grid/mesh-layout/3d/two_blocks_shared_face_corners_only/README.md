# Two Blocks Sharing a Face with Four Corner Nodes

This 3D case generates two adjacent unit cubes. Their shared face is represented
by its four corner nodes only; there are no interior nodes on that face.

Expected combined topology: 12 unique nodes, 2 hexahedral cells, 11 unique
faces, 10 boundary faces, and 1 internal face.

**This case verifies the basic shared-face connection between two 3D blocks.**
It does not test subdivision or interior nodes on the shared face.

## Run

Copy this directory to a writable workspace outside the repository and run:

```text
OneFLOW.exe r <case-directory>
```
