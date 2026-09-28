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