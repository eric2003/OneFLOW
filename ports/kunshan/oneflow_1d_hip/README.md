# OneFLOW 1D native HIP port

This isolated port validates the scalar one-dimensional first-order upwind
reconstruction and inviscid residual used by the active 1D development line.
It does not modify or replace the original Windows-oriented example and is not
a complete CFD GPU implementation.

The CPU executable checks the reference formula independently. The HIP
executable compares left/right reconstructed states and residuals against the
double-precision CPU reference with an absolute tolerance of `1e-15`.
