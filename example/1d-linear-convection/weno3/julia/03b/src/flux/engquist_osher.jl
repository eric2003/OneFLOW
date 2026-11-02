# src/flux/engquist_osher.jl
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

function compute!(calc::EngquistOsherFluxCalculator, qL::Vector{Float64}, qR::Vector{Float64}, flux::Vector{Float64})
    for i in 1:calc.mesh.nnodes
        c = calc.wave_speed
        cp = 0.5 * (c + abs(c))
        cm = 0.5 * (c - abs(c))
        u_L = qL[i]
        u_R = qR[i]
        flux[i] = cp * u_L + cm * u_R
    end
end