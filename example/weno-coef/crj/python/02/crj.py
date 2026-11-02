import matplotlib.pyplot as plt
import numpy as np

# 一维差分示意图
def plot_1d_fdm():
    # 定义网格点
    x = np.linspace(0, 1, 6)  # 6个网格点
    dx = x[1] - x[0]  # 网格间距

    # 绘制网格点
    plt.figure(figsize=(8, 2))
    plt.scatter(x, np.zeros_like(x), color='black', s=100, zorder=5)
    for i, xi in enumerate(x):
        plt.text(xi, -0.05, f'$x_{i}$', ha='center', fontsize=12)

    # 绘制网格线
    for i in range(len(x) - 1):
        plt.plot([x[i], x[i + 1]], [0, 0], color='black', linestyle='--', linewidth=1)

    # 添加标题和标签
    plt.title('One-Dimensional Finite Difference Grid', fontsize=14)
    plt.xlabel('x', fontsize=12)
    plt.yticks([])
    plt.grid(False)
    plt.show()

plot_1d_fdm()