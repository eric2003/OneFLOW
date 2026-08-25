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
