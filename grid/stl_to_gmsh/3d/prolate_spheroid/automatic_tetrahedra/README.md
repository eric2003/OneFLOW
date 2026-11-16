# 6:1 Prolate Spheroid - Automatic Tetrahedral Volume Mesh

This geometry pairs the native C++ OneFLOW Grid Case with an independent Python reference Case. Both use the same byte-identical 1280-facet STL and equivalent layout settings.

- `cpp/`: native OneFLOW Grid task, producing `results/prolate_spheroid_automatic_tetrahedra.msh`.
- `python/`: independent `converter.py` reference path, producing the same mesh artifact contract.

The spheroid is centered at the origin with semi-axes (3, 0.5, 0.5), full length 6, maximum diameter 1 and aspect ratio 6:1. Its analytic volume is pi. The STL is a triangulated affine image of the unit-sphere reference surface.

Both Cases use the generic closed, connected, consistently oriented genus-zero STL contract and are registered for the Mesh-to-CFD smoke route. The wall-only Physical Group does not define a complete external-flow benchmark.
