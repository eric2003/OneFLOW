# Module Interface Contract (draft; freeze after human confirmation)

> Items marked `[CONFIRM]` need a human decision. After freezing, changes require an RFC.

```cpp
namespace cfd {

class Field;   // layout in data_layout.md; surface in core_api.md
class Mesh;    // structured/unstructured mesh, read-only for modules
class Comm;    // halo exchange and global reductions

// Spatial operator. Sign convention [CONFIRM]: dU/dt = R(U, t), with R containing the sign.
class ISpatialOperator {
public:
  virtual ~ISpatialOperator() = default;

  // Compute residual R. Halo exchange is the spatial operator's responsibility. [CONFIRM]
  // U is read-only; R is the output.
  virtual void residual(const Field& U, double t, Field& R) const = 0;

  // Local stability limit (CFL=1). The time integrator reduces the global minimum via Comm.
  virtual double max_stable_dt(const Field& U, double t) const = 0;
};

// Time integrator. Explicit only in v0.1. Implicit needs an RFC for the Jacobian interface. [CONFIRM]
class ITimeIntegrator {
public:
  virtual ~ITimeIntegrator() = default;

  // Advance U in place; update t; dt is provided by the caller.
  // comm is used for the global reduction of max_stable_dt().
  virtual void advance(Field& U, double& t, double dt,
                       const ISpatialOperator& L, Comm& comm) = 0;

  virtual int order() const = 0;
};

} // namespace cfd
```

## Responsibility split

- Halo exchange: spatial, inside `residual()`, via core `Comm`.
- Global dt: temporal calls `max_stable_dt()` and reduces through `Comm`.
- Work arrays: temporal allocates extra `Field`s with `Field::like()`.
- Boundary conditions: spatial, expressed through `residual()`.
- `max_stable_dt` is a suggestion, not the final dt; temporal clips it.

## Decoupled tests

- temporal: mock `ISpatialOperator` ODE operators in `tests/mocks`, verify time order without spatial.
- spatial: MMS manufactured solution, compare `residual()` with the analytic time derivative, verify spatial order without temporal.
- integration: only on the integration branch, using real verification cases.

## Open decisions

- `[CONFIRM]` sign convention and whether R includes cell-volume normalization.
- `[CONFIRM]` halo exchange ownership and whether ghosts are written in place.
- `[CONFIRM]` boundary-condition ownership.
- `[CONFIRM]` whether `dt` is clipped by the integrator or caller.
