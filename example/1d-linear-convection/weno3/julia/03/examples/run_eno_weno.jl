# examples/run_eno_weno.jl
"""
1:1 复刻 run_eno_weno.py 的 Julia 版本
"""

include("../src/config.jl")
include("../src/mesh.jl")
include("../src/solver.jl")
include("../src/plotter.jl")

function performEnoWenoAnalysis()
    # 1. 初始化网格
    mesh = Mesh()
    
    # 2. 配置并运行 ENO3 求解
    println("Running ENO3 solver...")
    config_eno3 = CfdConfig()                    # 初始化默认配置
    with_reconstruction(config_eno3, "eno", 3)   # 显式指定 3 阶
    config_eno3.dt = 0.0025                      # 覆盖默认值
    config_eno3.rk_order = 2
    
    cfd_eno3 = Cfd(config_eno3, mesh)
    run!(cfd_eno3)  # 求解并生成 result 字典

    # 3. 配置并运行 WENO3 求解
    println("Running WENO3 solver...")
    config_weno3 = CfdConfig()
    with_reconstruction(config_weno3, "weno", 3)  # 显式指定 3 阶（WENO 默认 5 阶）
    config_weno3.dt = 0.0025
    config_weno3.rk_order = 2

    cfd_weno3 = Cfd(config_weno3, mesh)
    run!(cfd_weno3)

    # 5. 绘制 ENO/WENO 对比图
    println("Plotting comparison results...")
    plot_eno_weno_comparison(
        cfd_eno3.result,
        cfd_weno3.result;
        save_path="eno_weno_comparison.png"
    )

    return cfd_eno3, cfd_weno3
end

# 主程序入口
if abspath(PROGRAM_FILE) == @__FILE__
    performEnoWenoAnalysis()
    println("Analysis completed!")
end