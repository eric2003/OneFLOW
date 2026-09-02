# OneFLOW 1D native HIP port

This isolated port preserves the scalar one-dimensional first-order upwind
smoke test and adds a compressible one-dimensional Euler validation path. It
does not modify or replace the original Windows-oriented WENO3 example and is
not a complete CFD GPU implementation.

The Euler path advances all three conserved variables `[rho, rho*u, rho*E]`
with a first-order Rusanov flux and SSP-RK3. CPU and HIP share the
backend-neutral `EulerBackend::Step` seam; execution-specific details remain
behind the concrete backends. HIP is compared with CPU at every RK stage for face states, flux, residual, conserved and primitive variables. For the Sod discontinuity case, residual pointwise comparison excludes cells whose CPU reference `rho` or `p` jump by more than 2% to a neighbor; those cells are reported separately and are still covered by physical-state and state/flux checks.
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
GPU-ported. CUDA and Kokkos remain extension interfaces only, and Kokkos is not introduced here.
The current multi-GPU evidence is limited to the single-node 4-rank/4-DCU HIP MPI path.

HIP architecture selection is intentionally not hard-coded:

1. `CMAKE_HIP_ARCHITECTURES`, when supplied explicitly;
2. `ONEFLOW_HIP_ARCHITECTURES`, accepting a comma- or semicolon-separated list;
3. the architecture reported by `rocminfo` on the target compute node;
4. `rocm_agent_enumerator` only when it yields a unique candidate;
5. an explicit configuration error when no safe architecture can be determined.

This leaves room for multiple targets such as `gfx906;gfx936` without assuming
that Kunshan Z100 and Zhengzhou BW1000 are interchangeable. The current evidence
is recorded in `doc/reports/regression/oneflow-euler-validation.md`: Zhengzhou
`gfx936` is an archived reference only; the Kunshan Z100 `gfx906` automatic-detection
Euler closure passed on August 26, 2026. This remains an isolated one-dimensional
validation path, not a complete production Navier–Stokes backend.

The original Windows-oriented WENO3 example remains outside this isolated adapter
and is not overwritten.

## Euler backend benchmark

`oneflow_1d_euler_benchmark` runs the same one-dimensional Euler workload through
CPU and HIP and reports end-to-end time, HIP kernel time, H2D/D2H time,
allocation time, their ratios, and final-state error. GPU activity should be
checked together with the timing report using `rocm-smi` on the compute node.
The timing is intended to guide persistent device-state and reduced-D2H work;
it is not a substitute for the stage-by-stage correctness regression. On Kunshan
Z100 (`gfx906`), the current baseline reaches about 1.97--2.47x end-to-end speedup
for `nx=65536` through `4194304`, with final CPU/HIP state error 0. HIP event timing
shows the kernels execute, but kernels account for only about 1.56--2.51% of HIP
wall time while D2H dominates. `rocm-smi` sampling confirms activity on the target
HCU; the low average utilization is consistent with short kernels and synchronization,
not evidence of CPU fallback. The next optimization target is persistent device
state and fewer validation/trace copies, with correctness rechecked after each change.

The current implementation priority is the Euler stateful backend and persistent
device execution described in
`doc/reports/architecture/oneflow-euler-optimization-plan.md`. WENO5 and
dimensional expansion remain deferred until the Euler correctness and performance
acceptance gates pass.

`oneflow_1d_euler_stateful_benchmark` exercises the lifecycle API. The Kunshan
baseline uses one serial CPU thread and one DCU; for 100 steps, HIP steady advance
is 102--155x faster than the serial CPU across `nx=65536` through `4194304`, with
final state error 0. Kernel events account for 99.56--99.99% of `Advance` wall time,
and only one upload, one final download, and one synchronization occur per repeat.
This is the execution-structure result; the current multi-core/multi-card comparison is
recorded in the current performance report, while the full Navier-Stokes solver remains outside.
