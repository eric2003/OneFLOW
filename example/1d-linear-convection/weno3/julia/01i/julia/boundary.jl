# julia/boundary.jl
# 目标：与 Python boundary.py 行为完全一致（1:1 同构）
# 前提：ist/ied 在 Julia 中按 1-based 设置（ist = nghosts + 1）

using NPZ  # 仅用于测试，实际模块可不依赖

# ---------------------- 辅助函数 ----------------------
"""
模拟 Python 的 getattr(obj, name, default)
"""
function getfield_safe(obj, name::Symbol, default)
    return isdefined(obj, name) ? getfield(obj, name) : default
end

# ---------------------- 抽象基类（Julia 风格：用函数接口） ----------------------
# Julia 无 abstract class，用文档+约定
# 所有子类型必须实现 apply!(bc, u)

# ---------------------- PeriodicBoundary ----------------------
mutable struct PeriodicBoundary
    cfd::Any
end

function apply!(self::PeriodicBoundary, u::Vector{Float64})
    nghosts = self.cfd.domain.nghosts
    ist = self.cfd.domain.ist   # 1-based, e.g., 3
    ied = self.cfd.domain.ied   # 1-based, e.g., 43

    # 左 ghost: u[ist - 1 - ig] = u[ied - 1 - ig]
    for ig in 0:(nghosts-1)
        u[ist - 1 - ig] = u[ied - 1 - ig]
    end
    # 右 ghost: u[ied + ig] = u[ist + ig]
    for ig in 0:(nghosts-1)
        u[ied + ig] = u[ist + ig]
    end
end

# ---------------------- DirichletBoundary ----------------------
mutable struct DirichletBoundary
    cfd::Any
end

function apply!(self::DirichletBoundary, u::Vector{Float64})
    nghosts = self.cfd.domain.nghosts
    ist = self.cfd.domain.ist
    ied = self.cfd.domain.ied

    left_value = getfield_safe(self.cfd.config, :left_boundary_value, 1.0)
    right_value = getfield_safe(self.cfd.config, :right_boundary_value, 2.0)

    # 左边界
    for ig in 0:(nghosts-1)
        u[ist - 1 - ig] = left_value
    end
    # 右边界
    for ig in 0:(nghosts-1)
        u[ied + ig] = right_value
    end

    # 调试信息
    if getfield_safe(self.cfd.config, :debug, false)
        println("  应用Dirichlet边界: 左值=$left_value, 右值=$right_value")
    end
end

# ---------------------- NeumannBoundary ----------------------
mutable struct NeumannBoundary
    cfd::Any
end

function apply!(self::NeumannBoundary, u::Vector{Float64})
    nghosts = self.cfd.domain.nghosts
    ist = self.cfd.domain.ist
    ied = self.cfd.domain.ied

    # 左边界零梯度：u[ist - 1 - ig] = u[ist + ig]
    for ig in 0:(nghosts-1)
        u[ist - 1 - ig] = u[ist + ig]
    end
    # 右边界零梯度：u[ied + ig] = u[ied - 1 - ig]  ← 注意是 ied - 1 - ig
    for ig in 0:(nghosts-1)
        u[ied + ig] = u[ied - 1 - ig]
    end

    # 调试信息
    if getfield_safe(self.cfd.config, :debug, false)
        println("  应用Neumann边界: 零梯度")
    end
end