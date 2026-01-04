# julia/mesh.jl
"""
Mesh 结构体（与提供的 mesh.py 完全同构）
- 字段名、初始化顺序、循环逻辑完全一致
- 不做任何泛化或优化
"""

mutable struct Mesh
    xmin::Float64
    xmax::Float64
    ncells::Int
    nnodes::Int
    nx::Int
    x::Vector{Float64}
    xcc::Vector{Float64}
    L::Float64      # 在 init_mesh 中赋值
    dx::Float64     # 在 init_mesh 中赋值
    
    function Mesh()
        # 与 Python __init__ 完全一致
        xmin = 0.0
        xmax = 2.0
        ncells = 40
        nnodes = ncells + 1
        nx = ncells
        x = zeros(Float64, nnodes)
        xcc = zeros(Float64, ncells)
        
        # 调用 init_mesh（模拟 Python）
        L = xmax - xmin
        dx = L / ncells
        
        # Generate node coordinates: for i in range(self.nnodes)
        for i in 1:nnodes
            x[i] = xmin + (i - 1) * dx   # Python i → Julia i-1
        end
        
        # Generate cell center coordinates
        for i in 1:ncells
            xcc[i] = 0.5 * (x[i] + x[i+1])
        end
        
        new(xmin, xmax, ncells, nnodes, nx, x, xcc, L, dx)
    end
end