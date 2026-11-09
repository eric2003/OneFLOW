# STL to Gmsh: 2D Circle External Flow

This problem has two independent implementations:

- `cpp/`: native OneFLOW Grid Case using the C++ STL reader and Gmsh conversion pipeline
- `python/`: Python reference using the independent planar STL geometry implementation

Both use byte-identical 32-segment circular STL input, the same rectangular outer domain and matching mesh-size controls. The STL body is an internal `wall` boundary; the rectangular outer boundary has `inflow`, `outflow` and `farfield` groups. The fluid surface is the rectangle with the circular body loop subtracted as a hole and is named `fluid_domain` in both implementations.

The two meshes need not have identical node or element counts. Validation compares the geometry/domain definition, Physical Group semantics and mesh validity rather than requiring identical unstructured triangulations.
