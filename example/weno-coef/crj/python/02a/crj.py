import matplotlib.pyplot as plt
import numpy as np

# 二维差分示意图
def plot_2d_fdm():
    # 定义网格点
    x = np.linspace(0, 1, 5)  # 5个网格点
    y = np.linspace(0, 1, 5)  # 5个网格点
    dx = x[1] - x[0]  # 网格间距
    dy = y[1] - y[0]  # 网格间距

    # 创建网格
    X, Y = np.meshgrid(x, y)

    # 绘制网格点
    plt.figure(figsize=(6, 6))
    plt.scatter(X, Y, color='black', s=100, zorder=5)
    for i in range(len(x)):
        for j in range(len(y)):
            plt.text(X[j, i], Y[j, i], f'({i},{j})', ha='center', fontsize=10)

    # 绘制网格线
    for i in range(len(x)):
        plt.plot(X[i, :], Y[i, :], color='black', linestyle='--', linewidth=1)
    for j in range(len(y)):
        plt.plot(X[:, j], Y[:, j], color='black', linestyle='--', linewidth=1)

    # 添加标题和标签
    plt.title('Two-Dimensional Finite Difference Grid', fontsize=14)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.grid(False)
    plt.show()

plot_2d_fdm()