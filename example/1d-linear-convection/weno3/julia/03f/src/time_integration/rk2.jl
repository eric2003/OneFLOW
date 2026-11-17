# src/time_integration/rk2.jl

mutable struct RK2Integrator
    base::TimeIntegratorBase
end

function RK2Integrator(cfd::Any)
    RK2Integrator(TimeIntegratorBase(cfd))
end

function step(integrator::RK2Integrator, dt::Float64)
    base = integrator.base
    compute_residual(base)
    u_pred = copy(base.solution.u)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        u_pred[i] += dt * base.solution.res[j]
    end
    base.solution.u .= u_pred
    apply_boundary(base)
    compute_residual(base)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        base.solution.u[i] = 0.5 * base.solution.un[i] +
                             0.5 * base.solution.u[i] +
                             0.5 * dt * base.solution.res[j]
    end
    apply_boundary(base)
    update_old_field(base.solution)
end