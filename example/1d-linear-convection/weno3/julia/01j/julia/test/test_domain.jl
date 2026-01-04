# julia/test/test_domain.jl
include("../mesh.jl")
include("../domain.jl")

# MockConfig：模拟 Python CfdConfig
struct MockConfig
    recon_scheme::String
    spatial_order::Int
end

# 测试 ENO
config_eno = MockConfig("eno", 2)
mesh = Mesh()
domain_eno = Domain(config_eno, mesh)

println("ENO: nghosts = ", domain_eno.nghosts)  # 2
println("ENO: ist = ", domain_eno.ist)          # 2
println("ENO: ied = ", domain_eno.ied)          # 42
println("物理索引范围: ", collect(get_physical_indices(domain_eno))[1:3], " ... ", collect(get_physical_indices(domain_eno))[end-2:end])

# 测试 WENO（字符串 "weno"）
config_weno = MockConfig("weno", 2)
domain_weno = Domain(config_weno, mesh)
println("WENO: nghosts = ", domain_weno.nghosts)  # 2

# 测试 WENO3（字符串 "weno3"）
config_weno3 = MockConfig("weno3", 2)
domain_weno3 = Domain(config_weno3, mesh)
println("WENO3: nghosts = ", domain_weno3.nghosts)  # 2

# 测试 is_physical_cell
@assert is_physical_cell(domain_eno, 2) == true   # ist=2
@assert is_physical_cell(domain_eno, 41) == true  # ied-1=41
@assert is_physical_cell(domain_eno, 42) == false # ied=42

# ✅ 断言
@assert domain_eno.nghosts == 2
@assert domain_eno.ist == 2
@assert domain_eno.ied == 42
@assert domain_eno.ntcells == 44
@assert domain_weno.nghosts == 2
@assert domain_weno3.nghosts == 2

println("✅ Domain 测试通过！")