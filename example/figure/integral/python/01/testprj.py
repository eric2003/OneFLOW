import numpy as np
import matplotlib.pyplot as plt

# 定义积分区间和函数
a = 0
b = 1
n = 10  # 网格点数量
x = np.linspace(a, b, n+1)
f = lambda x: np.sin(2 * np.pi * x)  # 定义函数 f(x) = sin(2πx)

# 计算函数值
y = f(x)

# 绘制函数曲线
plt.figure(figsize=(8, 4))
plt.plot(x, y, label='f(x) = sin(2πx)', color='blue', linewidth=2)

# 绘制网格点
plt.scatter(x, y, color='red', zorder=5)

# 绘制积分区域
for i in range(n):
    plt.fill_between([x[i], x[i+1]], [0, 0], [y[i], y[i+1]], color='gray', alpha=0.3)

# 添加标题和标签
plt.title('One-Dimensional Line Integral', fontsize=14)
plt.xlabel('x', fontsize=12)
plt.ylabel('f(x)', fontsize=12)
plt.legend(fontsize=12)
plt.grid(True)

# 显示图表
plt.show()