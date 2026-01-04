#!/usr/bin/env python3
# example/run_eno_weno.py

import sys
import os

src_path = os.path.join(os.path.dirname(__file__), '..', 'src')
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from core.solver import Cfd
from infrastructure.config import CfdConfig
from infrastructure.mesh import Mesh
from visualization.plotter import plot_eno_weno_comparison, CFDPlotter
from physics.problems.linear_advection import LinearAdvectionProblem

def performEnoWenoAnalysisBAK():
    # 1. 初始化网格
    mesh = Mesh()

    # 2. 配置并运行 ENO3 求解
    print("Running ENO3 solver...")
    config_eno3 = CfdConfig()
    config_eno3.with_reconstruction("eno", 3)
    config_eno3.dt = 0.0025
    config_eno3.rk_order = 2

    # ✅ 创建 Problem 实例（替代直接传 config）
    problem_eno3 = LinearAdvectionProblem(config_eno3)
    cfd_eno3 = Cfd(problem_eno3, mesh)  # ← 注入 problem
    cfd_eno3.run()

    # 3. 配置并运行 WENO3 求解
    print("Running WENO3 solver...")
    config_weno3 = CfdConfig()
    config_weno3.with_reconstruction("weno", 3)
    config_weno3.dt = 0.0025
    config_weno3.rk_order = 2

    # ✅ 创建另一个 Problem 实例
    problem_weno3 = LinearAdvectionProblem(config_weno3)
    cfd_weno3 = Cfd(problem_weno3, mesh)  # ← 注入 problem
    cfd_weno3.run()

    # 4. 绘制对比图（完全不变）
    print("Plotting comparison results...")
    plot_eno_weno_comparison(
        eno_result=cfd_eno3.result,
        weno_result=cfd_weno3.result,
        save_path="eno_weno_comparison.png"
    )
    
def performEnoWenoAnalysis():
    # 1. 初始化网格
    mesh = Mesh()

    # 2. 配置并运行 ENO5 求解
    print("Running ENO5 solver...")
    config_eno5 = CfdConfig()
    config_eno5.with_reconstruction("eno", 5)
    config_eno5.dt = 0.0025
    config_eno5.rk_order = 2

    # ✅ 创建 Problem 实例（替代直接传 config）
    problem_eno5 = LinearAdvectionProblem(config_eno5)
    cfd_eno5 = Cfd(problem_eno5, mesh)  # ← 注入 problem
    cfd_eno5.run()

    # 3. 配置并运行 WENO5求解
    print("Running WENO5 solver...")
    config_weno5 = CfdConfig()
    config_weno5.with_reconstruction("weno", 5)
    config_weno5.dt = 0.0025
    config_weno5.rk_order = 2

    # ✅ 创建另一个 Problem 实例
    problem_weno5 = LinearAdvectionProblem(config_weno5)
    cfd_weno5 = Cfd(problem_weno5, mesh)  # ← 注入 problem
    cfd_weno5.run()

    # 4. 绘制对比图（完全不变）
    print("Plotting comparison results...")
    plot_eno_weno_comparison(
        eno_result=cfd_eno5.result,
        weno_result=cfd_weno5.result,
        save_path="eno_weno_comparison.png"
    )    


if __name__ == "__main__":
    performEnoWenoAnalysis()
    print("Analysis completed!")