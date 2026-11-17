# julia/test/test_eno.jl
include("../mesh.jl")
include("../domain.jl")
include("../solution.jl")
include("../reconstructor/eno.jl")

struct MockCfd
    domain::Domain
    solution::Solution
end

mesh = Mesh()
# 注意：Domain 使用 1-based 索引（ist = nghosts + 1）
domain = Domain((recon_scheme="eno", spatial_order=2), mesh)
solution = Solution((ic_type="step",), domain)

# 线性初始场: u[i] = (i-1) * 0.1
u = solution.u
for i in 1:domain.ntcells
    u[i] = Float64(i-1) * 0.1
end

eno = EnoReconstructor(2, domain.ntcells)
cfd = MockCfd(domain, solution)
reconstruct(eno, u, cfd)

println("q_face_left[1] = ", solution.q_face_left[1])
println("q_face_right[1] = ", solution.q_face_right[1])

# ENO2 对线性函数应精确重构 0.15
@assert abs(solution.q_face_left[1] - 0.15) < 1e-12
@assert abs(solution.q_face_right[1] - 0.15) < 1e-12

println("✅ ENO 测试通过！")