# julia/boundary.jl
# 目标：与 Python boundary.py 行为完全一致（1:1 同构）
# 前提：ist/ied 在 Julia 中按 1-based 设置（ist = nghosts + 1）

#using NPZ  # 仅用于测试，实际模块可不依赖

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


# ---------------------- BoundaryConditionFactory ----------------------
module BoundaryConditionFactory

const _REGISTRY = Dict{String, Function}()

function register(name::String, ctor::Function)
    if haskey(_REGISTRY, name)
        @warn "边界条件 '$name' 已注册，将被覆盖"
    end
    _REGISTRY[name] = ctor
end

function create(cfd::Any)
    if !hasproperty(cfd, :config)
        error("cfd 缺少 config 字段")
    end
    bc_type = cfd.config.boundary_type
    if !(bc_type isa AbstractString)
        error("boundary_type 必须为字符串，当前值: $bc_type")
    end

    if !haskey(_REGISTRY, bc_type)
        available = sort(collect(keys(_REGISTRY)))
        error("未注册的边界条件: '$bc_type'。可用选项: $available")
    end

    return _REGISTRY[bc_type](cfd)
end

# ✅ 使用 Main. 前缀引用顶层类型
register("periodic", cfd -> Main.PeriodicBoundary(cfd))
register("dirichlet", cfd -> Main.DirichletBoundary(cfd))
register("neumann", cfd -> Main.NeumannBoundary(cfd))

end  # module BoundaryConditionFactory

export BoundaryConditionFactory