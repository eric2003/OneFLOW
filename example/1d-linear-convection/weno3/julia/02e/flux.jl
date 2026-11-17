# julia/flux.jl
"""
通量计算器模块（与 flux.py 完全同构）
- 抽象基类 + 具体实现
- 字段：cfd, config, mesh, wave_speed
"""

include("mesh.jl")

# ---------------------- 抽象基类 ----------------------
"""
InviscidFluxCalculator 抽象类型
Julia 无 ABC，用文档约定
所有子类型必须实现 compute!
"""
abstract type InviscidFluxCalculator end

# ---------------------- RusanovFluxCalculator ----------------------
mutable struct RusanovFluxCalculator <: InviscidFluxCalculator
    cfd::Any
    config::Any
    mesh::Mesh
    wave_speed::Float64
    
    function RusanovFluxCalculator(cfd::Any)
        config = cfd.config
        mesh = cfd.domain.mesh
        wave_speed = config.wave_speed
        new(cfd, config, mesh, wave_speed)
    end
end

function compute!(calc::RusanovFluxCalculator, q_face_left::Vector{Float64}, q_face_right::Vector{Float64}, flux::Vector{Float64})
    for i in 1:calc.mesh.nnodes
        u_L = q_face_left[i]
        u_R = q_face_right[i]
        c_L = calc.wave_speed
        c_R = calc.wave_speed
        F_L = c_L * u_L
        F_R = c_R * u_R
        Smax = max(abs(c_L), abs(c_R))
        flux[i] = 0.5 * (F_L + F_R) - 0.5 * Smax * (u_R - u_L)
    end
end

# ---------------------- EngquistOsherFluxCalculator ----------------------
mutable struct EngquistOsherFluxCalculator <: InviscidFluxCalculator
    cfd::Any
    config::Any
    mesh::Mesh
    wave_speed::Float64
    
    function EngquistOsherFluxCalculator(cfd::Any)
        config = cfd.config
        mesh = cfd.domain.mesh
        wave_speed = config.wave_speed
        new(cfd, config, mesh, wave_speed)
    end
end

function compute!(calc::EngquistOsherFluxCalculator, q_face_left::Vector{Float64}, q_face_right::Vector{Float64}, flux::Vector{Float64})
    for i in 1:calc.mesh.nnodes
        c = calc.wave_speed
        cp = 0.5 * (c + abs(c))
        cm = 0.5 * (c - abs(c))
        u_L = q_face_left[i]
        u_R = q_face_right[i]
        flux[i] = cp * u_L + cm * u_R
    end
end

# ---------------------- FluxCalculatorFactory ----------------------
module FluxCalculatorFactory

const _REGISTRY = Dict{String, Function}()

function register(name::String, ctor::Function)
    if haskey(_REGISTRY, name)
        @warn "通量计算器 '$name' 已注册，将被覆盖"
    end
    _REGISTRY[name] = ctor
end

function create(cfd::Any)
    if !hasproperty(cfd, :config)
        error("cfd 缺少 config 字段")
    end
    config = cfd.config
    if !hasproperty(config, :flux_type)
        error("cfd.config 缺少 flux_type 字段")
    end

    flux_type = config.flux_type
    if !(flux_type isa AbstractString)
        error("flux_type 必须为字符串，当前值: $flux_type")
    end

    if !haskey(_REGISTRY, flux_type)
        available = sort(collect(keys(_REGISTRY)))
        error("未注册的通量计算器: '$flux_type'。可用选项: $available")
    end

    return _REGISTRY[flux_type](cfd)
end

# ✅ 修正：使用 Main. 前缀
register("rusanov", cfd -> Main.RusanovFluxCalculator(cfd))
register("engquist-osher", cfd -> Main.EngquistOsherFluxCalculator(cfd))

end  # module FluxCalculatorFactory

export FluxCalculatorFactory