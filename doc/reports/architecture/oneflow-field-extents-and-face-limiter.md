# Field Extents and Face-Limiter Reconstruction

## Scope

This note records the size-type contract behind the `MRField` face-state buffers and the role those buffers play in the finite-volume reconstruction path. The change is intentionally narrow: it aligns the second `Marray` constructor extent with `HXSize_t`; it does not change limiter mathematics or buffer reuse.

## Storage contract

`HXSize_t` is the unsigned size/count type used by OneFLOW and aliases `std::size_t`. `MRField` is an alias for `Marray<Real>`, whose storage is a vector of equation fields:

- The outer extent is the number of equations (`nEqu`).
- Each inner extent is the number of values stored for one equation (`nValues`); the values may represent cells, faces, or another domain-specific collection. The generic container should not assume that every inner extent is a cell count.
- For face-state buffers, the inner extent is the number of faces.

The `Marray` constructor previously accepted `HXSize_t` for the equation count but `int` for the inner extent. Passing a face count from `HXSize_t` therefore narrowed a size value to `int`, producing MSVC warning C4267 at the `std::make_unique<MRField>` call in `ULimiter.cpp`. The constructor now uses `HXSize_t` for both extents, matching the vector size API and avoiding a call-site cast that would merely conceal the mismatch. Its second parameter is named `nValues` rather than `numberOfCells`, because `Marray` is a generic equation-by-value container and also stores face data. This is a naming clarification, not a storage-layout change.

This is an interface-contract fix, not a claim that every index in the array subsystem has already been migrated. Index types and legacy loops can be reviewed separately when a concrete boundary or warning warrants it.

## Regression test

`tests/basic/MarrayTest.cpp` constructs an `Marray` using `HXSize_t` variables for both dimensions and checks the resulting equation count and per-equation value count. This makes the extent contract executable and puts the test in CTest through the `MarrayTest` target. It is intentionally a small contract test: it does not allocate extreme sizes or attempt to test allocator failure.

## Face-state buffer lifecycle

`LimField::BindFieldsAndAllocFaces` binds the cell-centered solution and gradient fields, then determines the equation count from the solution field. It retains `qf1` and `qf2` while both equation and face extents still match. If either buffer is absent or its extent differs, both are recreated with the required shape.

That reuse is important because the limiter is initialized repeatedly during Runge-Kutta stages. The face buffers are scratch storage, not persistent physical state: `GetQlQr` initializes face values from the adjacent left and right cells, reconstruction updates them, and boundary handling adjusts boundary-face values. Preserving the extent check avoids unnecessary allocation while keeping the shape tied to the current grid and equation system.

## Reconstruction algorithm overview

For each face, the method uses the left and right adjacent cells and the face-center displacement from each cell center. For a cell-side state, the linear reconstruction has the form

`q_face = q_cell + phi * (grad(q)_cell dot (x_face - x_cell))`

where `phi` is the cell limiter coefficient (or the minimum of the two side coefficients when the configured face reconstruction requires a shared limiter). A state-check callback accepts the candidate; when it rejects the candidate, the unmodified face value is retained. The weighted variant can blend the left and right states when a reconstructed candidate fails the check.

The limiter coefficient is computed from local neighbor extrema. The Barth limiter constrains extrapolation so the face prediction does not exceed the local allowable difference; the Venkatakrishnan variant applies a smooth limiter function with an epsilon related to the local solution range. Interface boundary faces are handled explicitly when building neighbor extrema, while ordinary interior faces update extrema for both adjacent cells.

## Architectural significance

The change clarifies a responsibility boundary:

1. Grid and field owners express counts using the shared size type.
2. The generic array container accepts those extents without narrowing.
3. Solver code remains responsible for deciding when scratch buffers can be reused.
4. Numerical reconstruction and limiter policy remain independent of the storage-type correction.

The goal is not to redesign the limiter or migrate every legacy index in one pass. It is to remove a real type mismatch at the reusable container boundary while preserving the algorithm and allocation behavior.

## Validation checklist

- Build the affected solver targets with MSVC and confirm C4267 is gone from the `Marray` construction path.
- Confirm the new `MarrayTest` is discovered by CTest and passes.
- Run the full configured test suite.
- Run the relevant CFD regression cases and compare numerical outputs if the build changes any template instantiations beyond the intended signature adjustment.
- Confirm `qf1` and `qf2` are still reused when equation and face extents match and reallocated when either extent changes.
