# julia/initial_condition.jl
"""
初始条件模块（Julia 版）
目标：与 Python initial_condition.py 行为完全一致
"""

# ---------------------- 辅助函数 ----------------------
"""
模拟 Python 的 getattr(obj, name, default)
"""
function getfield_safe(obj, name::Symbol, default)
    return isdefined(obj, name) ? getfield(obj, name) : default
end

# ---------------------- 抽象接口（Julia 无继承，用函数约定） ----------------------
# 所有初始条件必须实现：
#   evaluate_at(ic, x::AbstractVector{Float64}) -> Vector{Float64}
#   apply(ic, solution)

# ---------------------- StepFunctionIC ----------------------
struct StepFunctionIC
    config::Any
end

function evaluate_at(ic::StepFunctionIC, x::AbstractVector{Float64})
    u0 = ones(Float64, length(x))
    for i in eachindex(x)
        if 0.5 <= x[i] <= 1.0
            u0[i] = 2.0
        end
    end
    return u0
end

function apply(ic::StepFunctionIC, solution)
    x = solution.domain.mesh.xcc
    u0 = evaluate_at(ic, x)
    _apply_to_interior(solution, u0)
end

# ---------------------- SineWaveIC ----------------------
struct SineWaveIC
    config::Any
end

function evaluate_at(ic::SineWaveIC, x::AbstractVector{Float64})
    L = getfield_safe(ic.config, :domain_length, 2.0)
    return sin.(2π * x / L)
end

function apply(ic::SineWaveIC, solution)
    x = solution.domain.mesh.xcc
    u0 = evaluate_at(ic, x)
    _apply_to_interior(solution, u0)
end

# ---------------------- GaussianPulseIC ----------------------
struct GaussianPulseIC
    config::Any
end

function evaluate_at(ic::GaussianPulseIC, x::AbstractVector{Float64})
    center = getfield_safe(ic.config, :pulse_center, 0.5)
    width = getfield_safe(ic.config, :pulse_width, 0.1)
    return exp.(-((x .- center) ./ width).^2)
end

function apply(ic::GaussianPulseIC, solution)
    x = solution.domain.mesh.xcc
    u0 = evaluate_at(ic, x)
    _apply_to_interior(solution, u0)
end

# ---------------------- 公共辅助函数 ----------------------
"""
将初始场 values 应用到 solution.u 的物理区域（含 ghost 的数组）
"""
function _apply_to_interior(solution, values)
    domain = solution.domain
    # Python: for i in range(domain.ist, domain.ied)
    # Julia:  for i in domain.ist:domain.ied-1 （因为 ied 是结束索引，不包含）
    for i in domain.ist:(domain.ied - 1)
        j = i - domain.ist + 1  # values 是 1-based，i - ist 是 0-based → +1
        solution.u[i] = values[j]
    end
end