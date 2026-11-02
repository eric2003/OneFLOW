# src/time_integration/rk3.jl

mutable struct RK3Integrator
    base::TimeIntegratorBase
end

function RK3Integrator(cfd::Any)
    RK3Integrator(TimeIntegratorBase(cfd))
end

function step(integrator::RK3Integrator, dt::Float64)
    base = integrator.base
    # Stage 1
    compute_residual(base)
    u1 = copy(base.solution.u)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        u1[i] += dt * base.solution.res[j]
    end
    base.solution.u .= u1
    apply_boundary(base)
    # Stage 2
    compute_residual(base)
    u2 = copy(base.solution.u)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        u2[i] = 0.75 * base.solution.un[i] +
                0.25 * base.solution.u[i] +
                0.25 * dt * base.solution.res[j]
    end
    base.solution.u .= u2
    apply_boundary(base)
    # Stage 3
    compute_residual(base)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        base.solution.u[i] = c1 * base.solution.un[i] +
                             c2 * base.solution.u[i] +
                             c3 * dt * base.solution.res[j]
    end
    apply_boundary(base)
    update_old_field(base.solution)
end