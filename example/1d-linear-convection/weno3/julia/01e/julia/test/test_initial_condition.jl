# julia/test/test_initial_condition.jl
using NPZ

# 包含初始条件模块
include("../initial_condition.jl")

# 构造 mock config（与 Python 完全一致）
struct MockConfig
    ic_type::String
    domain_length::Float64
    pulse_center::Float64
    pulse_width::Float64
end

# 生成与 Python Mesh.xcc 完全相同的 x 坐标
function generate_xcc()
    xmin, xmax = 0.0, 2.0
    ncells = 40
    dx = (xmax - xmin) / ncells
    xcc = Vector{Float64}(undef, ncells)
    for i in 1:ncells
        xcc[i] = xmin + (i - 0.5) * dx  # i-1 + 0.5 → i-0.5
    end
    return xcc
end

# 主测试
xcc = generate_xcc()

test_cases = [
    ("step", StepFunctionIC(MockConfig("step", 2.0, 0.5, 0.1))),
    ("sin", SineWaveIC(MockConfig("sin", 2.0, 0.5, 0.1))),
    ("gaussian", GaussianPulseIC(MockConfig("gaussian", 2.0, 0.5, 0.1)))
]

for (name, ic) in test_cases
    u_jl = evaluate_at(ic, xcc)
    u_py = npzread("../../python/u_$(name)_interior_py.npy")
    
    err = maximum(abs.(u_jl .- u_py))
    println("IC: $(name) | 最大误差: $(err)")
    @assert err < 1e-12 "❌ $(name) 不一致！"
end

println("✅ 所有初始条件测试通过！")