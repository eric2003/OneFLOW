# Mesh Conversion Templates

This directory contains reusable, runnable Case templates for mesh-format
conversion routes. It does not contain the C++ converter implementations;
those remain in the corresponding reader, writer, and mesh-construction
modules.

## Route templates

| Template | Route | Responsibility |
|---|---|---|
| `cgns_to_oneflow` | CGNS -> OneFLOW | Read a CGNS mesh and write the native OneFLOW grid plus boundary/grid metadata |
| `gmsh_to_cgns` | Gmsh -> CGNS | Read Gmsh MSH input and write a CGNS interchange mesh with boundary semantics |
| `gmsh_to_oneflow` | Gmsh -> CGNS -> OneFLOW | Compose the Gmsh reader/writer and CGNS reader into a complete native-grid conversion |

Each template owns its runtime parameters, required input contract, and
instructions. Preserve these directories as independently stageable Cases.
Do not place generated outputs in the checked-in template directories.

## Composition in larger workflows

A workflow such as Mesh-to-CFD Validation can stage these templates into an
external workspace and bind artifact paths between stages:

```text
STL -> Gmsh producer
       |
       v
conversion/gmsh_to_cgns
       |
       v
CGNS interchange artifact
       |
       v
conversion/cgns_to_oneflow
       |
       v
native OneFLOW grid + boundary metadata
       |
       v
CFD consumer and validation evidence
```

The composed `gmsh_to_oneflow` template is also available when a workflow
needs the complete route as one Case. It uses the same Gmsh-to-CGNS mapping
and CGNS interchange boundary; it is not a second direct Gmsh converter.

## Focused examples

Examples that exercise a specific conversion contract belong beneath the
route that owns their execution. For example,
`gmsh_to_cgns/examples/convert_2d/two_quads` tests two adjacent QUAD_4 cells,
their shared internal edge, and Physical Group membership. It remains a
Gmsh-to-CGNS example rather than a generic conversion fixture.

## Inventory and implementation boundaries

- `grid/case_catalog.json` locates runnable Cases for GUI and CLI discovery.
- `grid/interoperability/conversion_matrix.json` records supported routes,
  templates, implementation status, and coverage.
- `grid/interoperability/` owns format descriptions and shared semantic
  contracts; it does not duplicate runnable route templates.
- C++/Python converter implementations remain in their implementation
  modules. Templates configure and exercise those implementations.

Tests and workflows must copy/stage templates into temporary or user-data
workspaces. Source fixtures and template parameter files remain unchanged.
