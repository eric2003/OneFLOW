import matplotlib.pyplot as plt
import numpy as np

# 定义网格点
i = 0  # 当前计算点
r = 2  # 左侧扩展范围
s = 2  # 右侧扩展范围

# 创建图形和子图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4), sharey=True)

# 左侧重建
x_left = np.arange(i - r, i + s + 1)
ax1.plot(x_left, np.zeros_like(x_left) + 1, 'o-', color='gray', markersize=8)
ax1.plot([i, i], [0, 1], 'k--')  # 绘制垂直线
ax1.plot([i + 1, i + 1], [0, 1], 'k--')  # 绘制垂直线
ax1.text(i + 0.5, 1.05, r'$u_{i+\frac{1}{2},L}^L$', fontsize=12, ha='center')
ax1.set_title('(a) Left-side reconstruction')
ax1.set_xticks(x_left)
ax1.set_xticklabels([f'{i-2}', f'{i-1}', f'{i}', f'{i+1}', f'{i+2}'])
ax1.set_xlim([i - r - 1, i + s + 2])
ax1.set_ylim([0, 1.5])
ax1.grid(False)

# 右侧重建
x_right = np.arange(i - 1, i + 4)
print("x_right=",x_right)
ax2.plot(x_right, np.zeros_like(x_right) + 1, 'o-', color='gray', markersize=8)
ax2.plot([i, i], [0, 1], 'k--')  # 绘制垂直线
ax2.plot([i + 2, i + 2], [0, 1], 'k--')  # 绘制垂直线
ax2.text(i + 1.5, 1.05, r'$u_{i+\frac{1}{2},R}^R$', fontsize=12, ha='center')
ax2.set_title('(b) Right-side reconstruction')
ax2.set_xticks(x_right)
ax2.set_xticklabels([f'{i-1}', f'{i}', f'{i+1}', f'{i+2}', f'{i+3}'])
ax2.set_xlim([i - 2, i + 3 + 1])
ax2.set_ylim([0, 1.5])
ax2.grid(False)

# 绘制边界点
for ax in [ax1, ax2]:
    ax.plot([i + s + 1, i + s + 1], [0, 1], 'ko')  # 绘制边界点
    ax.plot([i + s + 2, i + s + 2], [0, 1], 'ro')  # 绘制边界点
    ax.plot([i + s + 3, i + s + 3], [0, 1], 'ro')  # 绘制边界点

plt.tight_layout()
plt.show()