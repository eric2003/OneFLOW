# src/time_integration/base.jl

"""
公共基类逻辑（对应 Python TimeIntegratorBase）
"""

mutable struct TimeIntegratorBase
    cfd::Any
    config::Any
    domain::Domain
    solution::Solution
    residual_calculator::Any
end

function TimeIntegratorBase(cfd::Any)
    config = cfd.config
    domain = cfd.domain
    solution = cfd.solution
    residual_calculator = cfd.residual_calculator
    TimeIntegratorBase(cfd, config, domain, solution, residual_calculator)
end

function compute_residual(integrator::TimeIntegratorBase)
    compute!(integrator.residual_calculator)
end

function apply_boundary(integrator::TimeIntegratorBase)
    apply!(integrator.cfd.boundary_condition, integrator.solution.u)
end

function map_idx(integrator::TimeIntegratorBase, i::Int)
    return i - integrator.domain.ist + 1
end