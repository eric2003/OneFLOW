# OneFLOW 1D native HIP port

This isolated port preserves the scalar one-dimensional first-order upwind
smoke test and adds a compressible one-dimensional Euler validation path. It
does not modify or replace the original Windows-oriented WENO3 example and is
not a complete CFD GPU implementation.

The Euler path advances all three conserved variables `[rho, rho*u, rho*E]`
with a first-order Rusanov flux and SSP-RK3. HIP is compared with CPU at every
RK stage for face states, flux, residual, conserved and primitive variables. For the Sod discontinuity case, residual pointwise comparison excludes cells whose CPU reference `rho` or `p` jump by more than 2% to a neighbor; those cells are reported separately and are still covered by physical-state and state/flux checks.
The report includes absolute error, a scaled relative error, ULP distance,
finite and positive-state checks. Acceptance is
`abs(error) <= 1e-15 + 1e-15*max(1, abs(reference))`; the relative-error
scale is also reported with `max(1, abs(reference))`. This is a numerical
tolerance, not “15 decimal places” or “15 significant digits”.


## Isolated Lax-WENO5 replay

`OneDWeno5.cpp` is an independent Linux adapter that follows the existing
Windows-oriented `EulerField::LaxFriedrichs` path: local Lax wave speed,
WENO5 positive/negative flux reconstruction, and SSP-RK3. The original
`example/1d-linear-convection/weno3/cpp/01` tree is not modified.

`oneflow_1d_weno5_dump` writes a versioned binary CPU reference containing the
initial state and, for every time step and RK stage, split positive/negative
fluxes, reconstructed fluxes, numerical flux, residual, and state.
`oneflow_1d_weno5_replay` reads that artifact and advances HIP independently;
it compares all three conserved variables at every recorded intermediate, and
reports absolute error, scaled relative error, ULP distance, finite/physical
checks, and the residual smooth-region mask for the Sod case. The HIP path uses
a persistent per-case device state: allocation and the initial H2D transfer
happen once, RK stages reuse device buffers, and a production-style no-trace
run performs only one final D2H state transfer. The one-step API remains as a
compatibility wrapper around this persistent state.

This is a one-dimensional numerical replay validation, not a complete
compressible-CFD GPU implementation.
