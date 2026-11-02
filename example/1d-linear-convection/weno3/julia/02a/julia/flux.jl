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

# ---------------------- 通量计算器注册表 + 工厂（方案2）----------------------

"""
全局通量注册表：flux 名称 → 构造函数
"""
const FLUX_REGISTRY = Dict{String, Function}()

"""
注册通量计算器
"""
function register_flux(name::String, ctor::Function)
    if haskey(FLUX_REGISTRY, name)
        @warn "通量计算器 '$name' 已注册，将被覆盖"
    end
    FLUX_REGISTRY[name] = ctor
end

"""
工厂函数：根据 cfd.config.flux_type 创建通量计算器
"""
function create_flux_calculator(cfd::Any)
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

    if !haskey(FLUX_REGISTRY, flux_type)
        available = sort(collect(keys(FLUX_REGISTRY)))
        error("未注册的通量计算器: '$flux_type'。可用选项: $available")
    end

    return FLUX_REGISTRY[flux_type](cfd)
end

# ---------------------- 注册内置通量计算器 ----------------------
register_flux("rusanov", cfd -> RusanovFluxCalculator(cfd))
register_flux("engquist-osher", cfd -> EngquistOsherFluxCalculator(cfd))