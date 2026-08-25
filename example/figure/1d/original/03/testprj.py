import matplotlib.pyplot as plt
import numpy as np

# 定义当前计算点i和扩展范围r、s
i = 0  # 当前计算点
r = 3  # 左侧扩展范围
s = 3  # 右侧扩展范围

# 生成网格点索引
x_indices = np.arange(i - r, i + s + 1)

# 生成对应的值（这里用随机值表示，实际应用中可以替换为具体值）
u_values = np.random.rand(len(x_indices))

# 绘制网格点和对应的值
plt.figure(figsize=(10, 2))
plt.plot(x_indices, u_values, 'o-', label='Grid Points and Values')
plt.xlabel('Grid Index')
plt.ylabel('Values')
plt.title('Grid Points and Corresponding Values')
plt.grid(True)
plt.legend()

# 标注每个网格点
for idx, value in zip(x_indices, u_values):
    plt.text(idx, value, f'*', fontsize=12, ha='center', va='bottom')

# 标注当前计算点
plt.plot(i, u_values[i - (i - r)], 'ro', label='Current Point (i)')
plt.text(i, u_values[i - (i - r)], f'({i}, {u_values[i - (i - r)]:.2f})', fontsize=12, ha='center', va='bottom')

plt.show()