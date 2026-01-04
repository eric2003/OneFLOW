# julia/test/test_mesh.jl
include("../mesh.jl")

# 创建 mesh（与 Python 完全相同）
mesh = Mesh()

# 打印关键值（与 Python 对比）
println("xmin = ", mesh.xmin)      # 0.0
println("xmax = ", mesh.xmax)      # 2.0
println("ncells = ", mesh.ncells)  # 40
println("nnodes = ", mesh.nnodes)  # 41
println("nx = ", mesh.nx)          # 40
println("L = ", mesh.L)            # 2.0
println("dx = ", mesh.dx)          # 0.05

# 检查 x[1] (Python x[0]) 和 x[41] (Python x[40])
println("x[1] = ", mesh.x[1])      # 0.0
println("x[41] = ", mesh.x[41])    # 2.0

# 检查 xcc[1] (Python xcc[0]) 和 xcc[40] (Python xcc[39])
println("xcc[1] = ", mesh.xcc[1])   # 0.025
println("xcc[40] = ", mesh.xcc[40]) # 1.975

# ✅ 严格断言
@assert mesh.xmin == 0.0
@assert mesh.xmax == 2.0
@assert mesh.ncells == 40
@assert mesh.nnodes == 41
@assert mesh.nx == 40
@assert mesh.L == 2.0
@assert mesh.dx == 0.05
@assert mesh.x[1] == 0.0
@assert mesh.x[41] == 2.0
@assert abs(mesh.xcc[1] - 0.025) < 1e-12
@assert abs(mesh.xcc[40] - 1.975) < 1e-12

println("✅ Mesh 测试通过！")