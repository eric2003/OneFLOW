# src/equations/linear_advection.jl

"""
线性对流方程: ∂u/∂t + c ∂u/∂x = 0
"""
mutable struct LinearAdvectionEquation <: Equation
    wave_speed::Float64
    function LinearAdvectionEquation(config::Any)
        c = getfield_safe(config, :wave_speed, 1.0)
        new(c)
    end
end

eq_flux(eq::LinearAdvectionEquation, u::Float64) = eq.wave_speed * u

max_wave_speed(eq::LinearAdvectionEquation, u::Float64) = abs(eq.wave_speed)