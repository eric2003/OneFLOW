# src/equations/base.jl

"""
抽象方程基类
所有具体方程必须实现：
- eq_flux(eq, u)
- max_wave_speed(eq, u)
"""
abstract type Equation end

"""
通量函数 F(u)
"""
function eq_flux(eq::Equation, u::Float64)::Float64
    error("Not implemented for $(typeof(eq))")
end

"""
最大波速（用于 Rusanov 通量和 CFL）
"""
function max_wave_speed(eq::Equation, u::Float64)::Float64
    error("Not implemented for $(typeof(eq))")
end

"""
方程变量数（标量方程返回 1）
"""
num_equations(eq::Equation) = 1