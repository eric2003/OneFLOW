# julia/test/test_weno3.jl
include("../mesh.jl")
include("../domain.jl")
include("../solution.jl")
include("../reconstructor/weno3.jl")

# Mock CFD
struct MockCfd
    domain::Domain
    solution::Solution
end

# 创建数据
mesh = Mesh()
# 注意：Domain 现在使用 1-based 索引（ist = nghosts + 1）
domain = Domain((recon_scheme="weno3", spatial_order=2), mesh)
solution = Solution((ic_type="step",), domain)

# 设置 u（物理区域 + ghost）
u = solution.u
for i in 1:domain.ntcells
    u[i] = Float64(i-1) * 0.1  # u = [0.0, 0.1, 0.2, ..., 4.3]
end

# 创建 reconstructor
weno3 = Weno3Reconstructor()
cfd = MockCfd(domain, solution)

# 重建
reconstruct(weno3, u, cfd)

println("q_face_left length = ", length(solution.q_face_left))  # 应为 41
println("q_face_right length = ", length(solution.q_face_right)) # 应为 41
println("q_face_left[1] = ", solution.q_face_left[1])
println("q_face_right[1] = ", solution.q_face_right[1])

# 验证左界面值（i = ist-1 = 2, j = 1）
# v1 = u[1] = 0.0, v2 = u[2] = 0.1, v3 = u[3] = 0.2
# qL = w0*(-0.5*0.0 + 1.5*0.1) + w1*(0.5*0.1 + 0.5*0.2) = w0*0.15 + w1*0.15 = 0.15
@assert abs(solution.q_face_left[1] - 0.15) < 1e-12

# 验证右界面值（i = ist = 3, j = 1）
# v1 = u[2] = 0.1, v2 = u[3] = 0.2, v3 = u[4] = 0.3
# qR = w0*(0.5*0.1 + 0.5*0.2) + w1*(1.5*0.2 - 0.5*0.3) = w0*0.15 + w1*0.15 = 0.15
@assert abs(solution.q_face_right[1] - 0.15) < 1e-12

println("✅ WENO3 测试通过！")