# julia/test/test_boundary.jl
using NPZ
include("../boundary.jl")

struct MockConfig
    boundary_type::String
    left_boundary_value::Float64
    right_boundary_value::Float64
    debug::Bool
end

struct MockDomain
    nghosts::Int
    ist::Int   # Julia 本地索引（1-based）
    ied::Int
    ntcells::Int
end

struct MockCfd
    config::MockConfig
    domain::MockDomain
end

# ===== 关键：使用 Julia 本地索引规则 =====
nghosts = 2
ncells = 40
ist = nghosts + 1      # = 3
ied = ist + ncells     # = 43
ntcells = ncells + 2 * nghosts  # = 44

config = MockConfig("dirichlet", 0.5, 1.5, true)
domain = MockDomain(nghosts, ist, ied, ntcells)
cfd_mock = MockCfd(config, domain)

# 加载 Python 生成的 u_input.npy
# 注意：Python u[0] → Julia u[1]，所以内容完全对应
u_input = npzread("../../python/u_input.npy")
@assert length(u_input) == ntcells "数组长度不匹配！"

# 测试三种边界
test_cases = [
    ("periodic", PeriodicBoundary(cfd_mock)),
    ("dirichlet", DirichletBoundary(cfd_mock)),
    ("neumann", NeumannBoundary(cfd_mock))
]

for (name, bc) in test_cases
    u = copy(u_input)
    apply!(bc, u)

    u_py = npzread("../../python/u_$(name)_py.npy")
    err = maximum(abs.(u .- u_py))
    println("边界: $(name) | 最大误差: $(err)")
    @assert err < 1e-12 "❌ $(name) 不一致！"
end

println("✅ 所有边界条件测试通过！")