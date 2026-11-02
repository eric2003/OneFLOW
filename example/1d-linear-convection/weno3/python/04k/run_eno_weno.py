from solver import Cfd
from config import CfdConfig
from mesh import Mesh
from plotter import plot_eno_weno_comparison, CFDPlotter

def performEnoWenoAnalysis():
    # 1. 初始化网格
    #mesh = Mesh(ncells=100, L=2.0)
    mesh = Mesh()
    plotter = CFDPlotter()

    # 2. 配置并运行ENO3求解（使用你的链式调用）
    print("Running ENO3 solver...")
    config_eno3 = CfdConfig()  # 初始化默认配置
    config_eno3.with_reconstruction("eno", 3)  # 显式指定3阶（也可省略，ENO默认3阶）
    # 可选：覆盖默认值（如dt）
    config_eno3.dt = 0.0025
    config_eno3.rk_order = 1

    cfd_eno3 = Cfd(config_eno3, mesh)
    cfd_eno3.run()  # 求解并生成result字典

    # 可选：快速验证ENO3结果
    # plotter.plot_quick(cfd_eno3.result, title="ENO3 Quick Check")

    # 3. 配置并运行WENO3求解（注意：WENO默认5阶，这里显式指定3阶）
    print("Running WENO3 solver...")
    config_weno3 = CfdConfig()
    config_weno3.with_reconstruction("weno", 3)  # 显式指定3阶（默认是5阶）
    # 可选：覆盖默认值
    config_weno3.dt = 0.0025
    config_weno3.rk_order = 1

    cfd_weno3 = Cfd(config_weno3, mesh)
    cfd_weno3.run()

    # 4. 可选：保存结果（供离线绘图）
    # cfd_eno3.save_result("eno3_result.npz")
    # cfd_weno3.save_result("weno3_result.npz")

    # 5. 绘制ENO/WENO对比图
    print("Plotting comparison results...")
    plot_eno_weno_comparison(
        eno_result=cfd_eno3.result,
        weno_result=cfd_weno3.result,
        save_path="eno_weno_comparison.png"  # 可选：保存图片
    )

if __name__ == "__main__":
    # 主程序入口
    performEnoWenoAnalysis()
    print("Analysis completed!")