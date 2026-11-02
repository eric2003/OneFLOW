# src/problems/linear_advection.jl

"""
线性对流问题：定义初始条件、边界条件、解析解
"""
mutable struct LinearAdvectionProblem <: Problem
    ic_type::String
    boundary_type::String
    left_value::Float64
    right_value::Float64
    domain_length::Float64
    wave_speed::Float64  # 从 config 获取，未来可从 equation 获取

    function LinearAdvectionProblem(config::Any)
        new(
            getfield_safe(config, :ic_type, "step"),
            getfield_safe(config, :boundary_type, "periodic"),
            getfield_safe(config, :left_boundary_value, 1.0),
            getfield_safe(config, :right_boundary_value, 2.0),
            getfield_safe(config, :domain_length, 2.0),
            getfield_safe(config, :wave_speed, 1.0)
        )
    end
end

function create_initial_condition(prob::LinearAdvectionProblem, config::Any)
    return Main.InitialConditionFactory.create(prob.ic_type, config)
end

function create_boundary_condition(prob::LinearAdvectionProblem, cfd::Any)
    # 复用原有 BC 工厂（cfd 需包含 config + domain）
    return Main.BoundaryConditionFactory.create(cfd)
end

function exact_solution(prob::LinearAdvectionProblem, x::Vector{Float64}, t::Float64)::Vector{Float64}
    L = prob.domain_length
    c = prob.wave_speed
    # 周期性平移
    x_shifted = @. (x - c * t + L) % L
    # 重用 IC 的 evaluate_at
    ic = Main.InitialConditionFactory.create(prob.ic_type, (ic_type=prob.ic_type,))
    return Main.evaluate_at(ic, x_shifted)
end