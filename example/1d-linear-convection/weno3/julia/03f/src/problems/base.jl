# src/problems/base.jl

"""
抽象问题基类
所有具体问题必须实现：
- create_initial_condition(prob, config)
- create_boundary_condition(prob, cfd)
- exact_solution(prob, x, t)
"""
abstract type Problem end

function create_initial_condition(prob::Problem, config::Any)
    error("Not implemented for $(typeof(prob))")
end

function create_boundary_condition(prob::Problem, cfd::Any)
    error("Not implemented for $(typeof(prob))")
end

function exact_solution(prob::Problem, x::Vector{Float64}, t::Float64)::Vector{Float64}
    error("Not implemented for $(typeof(prob))")
end