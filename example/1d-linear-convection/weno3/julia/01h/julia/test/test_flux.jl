# julia/test/test_flux.jl
include("../mesh.jl")
include("../flux.jl")

# Mock 对象
struct MockConfig
    flux_type::String
    wave_speed::Float64
end

struct MockDomain
    mesh::Mesh
end

struct MockCfd
    config::MockConfig
    domain::MockDomain
end

# 创建 mock
mesh = Mesh()
config = MockConfig("rusanov", 1.0)
domain = MockDomain(mesh)
cfd = MockCfd(config, domain)

# 测试 Rusanov
rusanov = RusanovFluxCalculator(cfd)
N = mesh.nnodes
qL = [1.0, 2.0, 3.0, 4.0, 5.0, zeros(Float64, N-5)...]
qR = [2.0, 3.0, 4.0, 5.0, 6.0, zeros(Float64, N-5)...]
flux_rus = zeros(Float64, N)
compute!(rusanov, qL, qR, flux_rus)

println("Rusanov flux[1] = ", flux_rus[1])  # 应为 1.0
@assert abs(flux_rus[1] - 1.0) < 1e-12

# 测试 Engquist-Osher
eo = EngquistOsherFluxCalculator(cfd)
flux_eo = zeros(Float64, N)
compute!(eo, qL, qR, flux_eo)

println("EO flux[1] = ", flux_eo[1])  # c=1 → cp=1, cm=0 → flux=1*1 + 0*2 = 1.0
@assert abs(flux_eo[1] - 1.0) < 1e-12

# 测试 c = -1.0
config_neg = MockConfig("rusanov", -1.0)
cfd_neg = MockCfd(config_neg, domain)
rusanov_neg = RusanovFluxCalculator(cfd_neg)
compute!(rusanov_neg, qL, qR, flux_rus)
println("Rusanov (c=-1) flux[1] = ", flux_rus[1])  # F_L=-1, F_R=-2, Smax=1 → flux = -1.5 -0.5*(-1) = -1.0
@assert abs(flux_rus[1] - (-2.0)) < 1e-12  # ← 修正：-2.0 而非 -1.0

println("✅ Flux 测试通过！")