# src/time_integration/rk1.jl

mutable struct RK1Integrator
    base::TimeIntegratorBase
end

function RK1Integrator(cfd::Any)
    RK1Integrator(TimeIntegratorBase(cfd))
end

function step(integrator::RK1Integrator, dt::Float64)
    compute_residual(integrator.base)
    for i in integrator.base.domain.ist:(integrator.base.domain.ied - 1)
        j = map_idx(integrator.base, i)
        integrator.base.solution.u[i] += dt * integrator.base.solution.res[j]
    end
    apply_boundary(integrator.base)
    update_old_field(integrator.base.solution)
end