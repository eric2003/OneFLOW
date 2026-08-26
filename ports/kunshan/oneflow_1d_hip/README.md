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

## Current status and architecture policy

The isolated directory is the current one-dimensional validation path. CPU remains
the default reference and execution path; the HIP targets are validation adapters
and reserved integration points, not a claim that the full Navier–Stokes solver is
GPU-ported. CUDA, Kokkos, and multi-GPU remain extension interfaces only, and
Kokkos is not introduced here.

HIP architecture selection is intentionally not hard-coded:

1. `CMAKE_HIP_ARCHITECTURES`, when supplied explicitly;
2. `ONEFLOW_HIP_ARCHITECTURES`, accepting a comma- or semicolon-separated list;
3. the architecture reported by `rocminfo` on the target compute node;
4. `rocm_agent_enumerator` only when it yields a unique candidate;
5. an explicit configuration error when no safe architecture can be determined.

This leaves room for multiple targets such as `gfx906;gfx936` without assuming
that Kunshan Z100 and Zhengzhou BW1000 are interchangeable. The current evidence
is recorded in `doc/reports/one-dimensional-euler-hip-validation.md`: Zhengzhou
`gfx936` is an archived reference only; the Kunshan Z100 `gfx906` automatic-detection
Euler closure passed on August 26, 2026. This remains an isolated one-dimensional
validation path, not a complete production Navier–Stokes backend.

The original Windows-oriented WENO3 example remains outside this isolated adapter
and is not overwritten.
