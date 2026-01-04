# julia/test/test_time_integration.jl
include("../time_integration.jl")
include("../flux.jl")
include("../boundary.jl")

# Mock CFD 对象
struct MockCfd
    config::Any
    domain::Domain
    solution::Solution
    residual_calculator::Any
    boundary_condition::Any
end

# 创建完整链路
mesh = Mesh()
domain = Domain((recon_scheme="eno", spatial_order=2), mesh)
solution = Solution((ic_type="step",), domain)

# 手动设置界面值（跳过 reconstructor）
for i in 1:mesh.nnodes
    solution.q_face_left[i] = Float64(i) * 0.1
    solution.q_face_right[i] = Float64(i) * 0.1
end

# Mock components
flux_calc = RusanovFluxCalculator((config=(wave_speed=1.0,), domain=domain))
res_calc = ResidualCalculator((config=nothing, domain=domain, solution=solution, reconstructor=nothing), flux_calc)
bc = PeriodicBoundary((domain=domain, config=(debug=false,)))
cfd = MockCfd((rk_order=1,), domain, solution, res_calc, bc)

# 测试 RK1
rk1 = RK1Integrator(cfd)
step(rk1, 0.025)

println("RK1 step completed")
@assert solution.u[3] != 0.0  # 应已更新

# 测试 RK2
cfd2 = MockCfd((rk_order=2,), domain, solution, res_calc, bc)
rk2 = RK2Integrator(cfd2)
step(rk2, 0.025)

println("RK2 step completed")

# 测试 RK3
cfd3 = MockCfd((rk_order=3,), domain, solution, res_calc, bc)
rk3 = RK3Integrator(cfd3)
step(rk3, 0.025)

println("RK3 step completed")

println("✅ Time Integration 逻辑测试通过！")