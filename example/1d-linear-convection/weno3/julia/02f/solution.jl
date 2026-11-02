# julia/solution.jl
"""
Solution 结构体（与真实 solution.py 完全同构）
字段顺序、方法、逻辑 1:1 对齐
"""

include("mesh.jl")
include("domain.jl")
include("initial_condition.jl")

mutable struct Solution
    domain::Domain
    q_face_left::Vector{Float64}
    q_face_right::Vector{Float64}
    flux::Vector{Float64}
    res::Vector{Float64}
    u::Vector{Float64}
    un::Vector{Float64}
    
    function Solution(config::Any, domain::Domain)
        mesh = domain.mesh
        
        q_face_left = zeros(Float64, mesh.nnodes)
        q_face_right = zeros(Float64, mesh.nnodes)
        flux = zeros(Float64, mesh.nnodes)
        res = zeros(Float64, mesh.ncells)
        u = zeros(Float64, domain.ntcells)
        un = zeros(Float64, domain.ntcells)
        
        sol = new(domain, q_face_left, q_face_right, flux, res, u, un)
        initialize_from_config(sol, config)
        return sol
    end
end

"""
重置解数组为初始状态
"""
function reset_solution(sol::Solution)
    fill!(sol.u, 0.0)
    fill!(sol.un, 0.0)
end

"""
根据配置初始化场
"""

function initialize_from_config(sol::Solution, config::Any)
    ic_type = getfield_safe(config, :ic_type, "step")
    ic = InitialConditionFactory.create(ic_type, config)
    apply(ic, sol)
end

"""
更新旧场：un = u
"""
function update_old_field(sol::Solution)
    sol.un .= sol.u  # 等价于 Python 的 un[:] = u[:]
end

# ---------------------- 辅助函数 ----------------------
function getfield_safe(obj, name::Symbol, default)
    return isdefined(obj, name) ? getfield(obj, name) : default
end