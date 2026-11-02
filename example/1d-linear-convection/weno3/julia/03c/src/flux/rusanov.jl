# src/flux/rusanov.jl

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