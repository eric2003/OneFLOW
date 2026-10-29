# Field Extents and Face-Limiter Reconstruction

## Scope

This note records the size-type contract behind the `MRField` face-state buffers and the workspace lifecycle in the finite-volume reconstruction path. The changes are intentionally narrow: they align the `Marray` constructor extents with `HXSize_t` and reuse limiter extrema storage without changing the limiter formulas or their initialization semantics.

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

## Limiter extrema workspace lifecycle

`Limiter::CalcLimiter` computes extrema for each equation before evaluating Barth or Venkatakrishnan coefficients. Previously, it constructed `minvf` and `maxvf` at every call and destroyed them after processing all equations. The arrays are scratch workspace: their values are rebuilt for the next limiter calculation, so their storage can persist even though their contents cannot.

The implementation now keeps these two vectors in the existing RAII-owned `Lim` object. Before each calculation, `assign(ug.nTCell, 0.0)` resets every entry and adjusts the logical extent to the current grid. When the grid size is unchanged, the vector can retain its allocated capacity and avoid repeated heap allocation. If the extent changes, the vector interface still accommodates the new size.

The full reset is deliberate. `CalcMinMaxDiff` seeds extrema from physical cells and then updates them through interface and interior faces. Other slots, including ghost-cell slots, were previously zero-initialized by constructing a fresh vector. Reusing old contents without clearing them would change those slots and could alter boundary-face limiter inputs. This change preserves the existing zero-initialization behavior rather than attempting to reinterpret or repair the ghost-cell algorithm.

This is a workspace-lifetime optimization, not a numerical-method change. It removes repeated allocation/deallocation while retaining the same initial values for the extrema arrays on every call. A future numerical change to ghost extrema should be handled separately, with boundary/interface tests and solver regression evidence.

## Reconstruction algorithm overview

For each face, the method uses the left and right adjacent cells and the face-center displacement from each cell center. For a cell-side state, the linear reconstruction has the form

`q_face = q_cell + phi * (grad(q)_cell dot (x_face - x_cell))`

where `phi` is the cell limiter coefficient (or the minimum of the two side coefficients when the configured face reconstruction requires a shared limiter). A state-check callback accepts the candidate; when it rejects the candidate, the unmodified face value is retained. The weighted variant can blend the left and right states when a reconstructed candidate fails the check.

The limiter coefficient is computed from local neighbor extrema. The Barth limiter constrains extrapolation so the face prediction does not exceed the local allowable difference; the Venkatakrishnan variant applies a smooth limiter function with an epsilon related to the local solution range. Interface boundary faces are handled explicitly when building neighbor extrema, while ordinary interior faces update extrema for both adjacent cells.

## Architectural significance

The change clarifies a responsibility boundary:

1. Grid and field owners express counts using the shared size type.
2. The generic array container accepts those extents without narrowing.
3. The limiter owns its extrema workspace and controls when its capacity can be reused; contents are reset at every calculation.
4. Numerical reconstruction and limiter policy remain independent of the storage-type correction.

The goal is not to redesign the limiter or migrate every legacy index in one pass. It is to remove a real type mismatch at the reusable container boundary and reduce repeated workspace allocation while preserving the numerical formulas and per-call initialization semantics.

## Validation checklist

- Build the affected solver targets with MSVC and confirm C4267 is gone from the `Marray` construction path.
- Confirm the new `MarrayTest` is discovered by CTest and passes.
- Run the full configured test suite.
- Run the relevant CFD regression cases and compare numerical outputs if the build changes any template instantiations beyond the intended signature adjustment.
- Confirm `qf1` and `qf2` are still reused when equation and face extents match and reallocated when either extent changes.
- Confirm `minvf` and `maxvf` are reused across limiter calls, reset to zero over the full `nTCell` extent each time, and resize correctly when the grid extent changes.
