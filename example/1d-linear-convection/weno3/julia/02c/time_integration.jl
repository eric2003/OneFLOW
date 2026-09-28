# julia/time_integration.jl
"""
时间推进器模块（与 time_integration.py 完全同构）
"""

include("mesh.jl")
include("domain.jl")
include("solution.jl")
include("residual.jl")
include("boundary.jl")

# ---------------------- 抽象时间推进器基类 ----------------------
abstract type TimeIntegrator end

mutable struct TimeIntegratorBase <: TimeIntegrator
    cfd::Any
    config::Any
    domain::Domain
    solution::Solution
    residual_calculator::Any  # ResidualCalculator
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
    return i - integrator.domain.ist + 1  # ← +1 转为 1-based
end

# ---------------------- RK1Integrator ----------------------
mutable struct RK1Integrator <: TimeIntegrator
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

# ---------------------- RK2Integrator ----------------------
mutable struct RK2Integrator <: TimeIntegrator
    base::TimeIntegratorBase
end

function RK2Integrator(cfd::Any)
    RK2Integrator(TimeIntegratorBase(cfd))
end

function step(integrator::RK2Integrator, dt::Float64)
    base = integrator.base
    # 阶段1：预测步
    compute_residual(base)
    u_pred = copy(base.solution.u)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        u_pred[i] += dt * base.solution.res[j]
    end
    base.solution.u .= u_pred
    apply_boundary(base)
    # 阶段2：校正步
    compute_residual(base)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        base.solution.u[i] = 0.5 * base.solution.un[i] + 0.5 * base.solution.u[i] + 0.5 * dt * base.solution.res[j]
    end
    apply_boundary(base)
    update_old_field(base.solution)
end

# ---------------------- RK3Integrator ----------------------
mutable struct RK3Integrator <: TimeIntegrator
    base::TimeIntegratorBase
end

function RK3Integrator(cfd::Any)
    RK3Integrator(TimeIntegratorBase(cfd))
end

function step(integrator::RK3Integrator, dt::Float64)
    base = integrator.base
    # 阶段1
    compute_residual(base)
    u1 = copy(base.solution.u)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        u1[i] += dt * base.solution.res[j]
    end
    base.solution.u .= u1
    apply_boundary(base)
    # 阶段2
    compute_residual(base)
    u2 = copy(base.solution.u)
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        u2[i] = 0.75 * base.solution.un[i] + 0.25 * base.solution.u[i] + 0.25 * dt * base.solution.res[j]
    end
    base.solution.u .= u2
    apply_boundary(base)
    # 阶段3
    compute_residual(base)
    c1, c2, c3 = 1.0/3.0, 2.0/3.0, 2.0/3.0
    for i in base.domain.ist:(base.domain.ied - 1)
        j = map_idx(base, i)
        base.solution.u[i] = c1 * base.solution.un[i] + c2 * base.solution.u[i] + c3 * dt * base.solution.res[j]
    end
    apply_boundary(base)
    update_old_field(base.solution)
end

# ---------------------- TimeIntegratorFactory ----------------------
module TimeIntegratorFactory

const _REGISTRY = Dict{String, Function}()

function register(name::String, ctor::Function)
    if haskey(_REGISTRY, name)
        @warn "时间积分器 '$name' 已注册，将被覆盖"
    end
    _REGISTRY[name] = ctor
end

function create(cfd::Any)
    if !hasproperty(cfd, :config)
        error("cfd 缺少 config 字段")
    end
    rk_order = cfd.config.rk_order
    if !(rk_order isa Integer) || rk_order < 1
        error("rk_order 必须为正整数，当前值: $rk_order")
    end

    name = "rk$rk_order"
    if !haskey(_REGISTRY, name)
        available = sort(collect(keys(_REGISTRY)))
        error("未注册的时间积分器: '$name'。可用选项: $available")
    end

    return _REGISTRY[name](cfd)
end

# ✅ 使用 Main. 前缀引用顶层类型
register("rk1", cfd -> Main.RK1Integrator(cfd))
register("rk2", cfd -> Main.RK2Integrator(cfd))
register("rk3", cfd -> Main.RK3Integrator(cfd))

end  # module TimeIntegratorFactory

export TimeIntegratorFactory