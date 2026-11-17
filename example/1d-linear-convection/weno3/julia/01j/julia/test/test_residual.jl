# julia/test/test_residual.jl
include("../mesh.jl")
include("../domain.jl")
include("../solution.jl")
include("../flux.jl")
include("../residual.jl")

# ===== Dummy Reconstructor =====
struct DummyReconstructor end

# ===== Mock Config =====
struct MockConfig
    flux_type::String
    wave_speed::Float64
end

# ===== Mock CFD (必须包含 reconstructor 字段) =====
struct MockCfd
    config::MockConfig
    domain::Domain
    solution::Solution
    reconstructor::DummyReconstructor  # ← 关键：添加此字段
end

# ===== 主测试 =====
config = MockConfig("rusanov", 1.0)
mesh = Mesh()
domain = Domain((recon_scheme="eno", spatial_order=2), mesh)
solution = Solution((ic_type="step",), domain)

# 手动设置界面值（跳过重建）
for i in 1:mesh.nnodes
    solution.q_face_left[i] = Float64(i) * 0.1
    solution.q_face_right[i] = Float64(i) * 0.1
end

flux_calc = RusanovFluxCalculator((config=config, domain=domain))
dummy_recon = DummyReconstructor()
cfd = MockCfd(config, domain, solution, dummy_recon)

# 创建残差计算器
res_calc = ResidualCalculator(cfd, flux_calc)

# 计算残差（注意：_reconstruct 仍被注释）
compute!(res_calc)

println("flux[1] = ", solution.flux[1])
println("res[1] = ", solution.res[1])

@assert abs(solution.flux[1] - 0.1) < 1e-12
@assert abs(solution.res[1] - (-2.0)) < 1e-12

println("✅ Residual 测试通过！")