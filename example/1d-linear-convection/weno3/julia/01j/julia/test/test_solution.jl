# julia/test/test_solution.jl
include("../mesh.jl")
include("../domain.jl")
include("../solution.jl")

# MockConfig
struct MockConfig
    recon_scheme::String
    spatial_order::Int
    ic_type::String
    domain_length::Float64
    pulse_center::Float64
    pulse_width::Float64
end

# 创建 solution
config = MockConfig("eno", 2, "step", 2.0, 0.5, 0.1)
mesh = Mesh()
domain = Domain(config, mesh)
sol = Solution(config, domain)

# 检查字段尺寸
@assert length(sol.q_face_left) == mesh.nnodes  # 41
@assert length(sol.flux) == mesh.nnodes         # 41
@assert length(sol.res) == mesh.ncells          # 40
@assert length(sol.u) == domain.ntcells         # 44

# 检查初始场
println("u[3] (物理起始): ", sol.u[3])   # 应为 1.0
println("u[23] (x=1.0): ", sol.u[23])    # 应为 2.0

# 测试 update_old_field
sol.u[3] = 999.0
update_old_field(sol)
@assert sol.un[3] == 999.0

# 测试 reset_solution
reset_solution(sol)
@assert sol.u[3] == 0.0
@assert sol.un[3] == 0.0

println("✅ Solution 测试通过！")