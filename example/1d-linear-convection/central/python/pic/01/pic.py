import matplotlib.pyplot as plt
import numpy as np

# Matplotlib 中文支持 + 负号修复
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 参数（可调整N查看效果）
N = 9  # 节点数，取奇数便于对称显示

# 图1：一维网格示意图（节点1到N，单元标注）
fig1, ax1 = plt.subplots(figsize=(10, 3))
x_nodes = np.arange(1, N+1)
ax1.plot(x_nodes, np.zeros(N), 'o', markersize=12, color='blue', label='节点 (nodes)')
for i in range(1, N):
    ax1.text(i + 0.5, 0.05, f'cell {i}', ha='center', fontsize=12, color='red')
ax1.axhline(0, color='black', linewidth=1)
ax1.set_ylim(-0.5, 0.5)
ax1.set_xlim(0.5, N+0.5)
ax1.set_yticks([])
ax1.set_xlabel('空间方向 x')
for i, txt in enumerate(range(1, N+1)):
    ax1.text(x_nodes[i], -0.1, f'{txt}', ha='center', fontsize=12)
ax1.set_title('一维计算网格示意图（节点编号 1 到 N，单元在节点之间）')
ax1.legend()

# 图2：空间中心差分 stencil（内点 + 边界周期绕回）
fig2, ax2 = plt.subplots(figsize=(12, 4))
x_stencil = np.arange(1, N+1)
ax2.plot(x_stencil, np.zeros(N), 'o', markersize=15, color='blue')

# 内点示例（取中间点 i=5）
i_center = 5
ax2.plot(x_stencil[i_center-1], 0, 'o', markersize=20, color='red', label='计算点 i')
ax2.annotate('i-1', (x_stencil[i_center-2], 0.1), ha='center', fontsize=14, color='green')
ax2.annotate('i', (x_stencil[i_center-1], 0.1), ha='center', fontsize=14, color='red')
ax2.annotate('i+1', (x_stencil[i_center], 0.1), ha='center', fontsize=14, color='green')

# 边界点示例（i=1 和 i=N，周期绕回）
ax2.arrow(x_stencil[-1], -0.2, - (N-2), 0, head_width=0.05, color='purple', linewidth=2, label='周期绕回 (i=1 的 i-1 → N)')
ax2.arrow(x_stencil[0], -0.3, N-2, 0, head_width=0.05, color='purple', linewidth=2, label='周期绕回 (i=N 的 i+1 → 1)')

ax2.set_ylim(-0.5, 0.5)
ax2.set_yticks([])
ax2.set_title('空间中心差分计算模板（stencil）及周期边界绕回')
ax2.legend()

# 图3：Crank-Nicolson 时空 stencil
fig3, ax3 = plt.subplots(figsize=(8, 6))
# 时间向上，空间水平
t_levels = [0, 1]  # n 层和 n+1 层
for t in t_levels:
    ax3.plot(x_stencil, np.full(N, t), 'o', markersize=12, color='blue' if t==0 else 'red')
    ax3.text(0, t, f'时间层 {"n" if t==0 else "n+1"}', fontsize=14)

# 示例内点 stencil（菱形）
i = 5
ax3.plot([x_stencil[i-2], x_stencil[i-1], x_stencil[i]], [0,1,0], 'k--', linewidth=1)
ax3.plot([x_stencil[i-2], x_stencil[i-1], x_stencil[i]], [1,0,1], 'k--', linewidth=1)
ax3.text(x_stencil[i-1], 0.5, 'CN 平均\n(时间+空间中心)', ha='center', bbox=dict(facecolor='yellow', alpha=0.5), fontsize=12)

ax3.set_yticks([0,1])
ax3.set_yticklabels(['t^n', 't^{n+1}'])
ax3.set_xlabel('空间节点 i')
ax3.set_title('Crank-Nicolson 时空计算模板（菱形平均）')

# 图4：矩阵结构示意（循环三对角）
fig4, ax4 = plt.subplots(figsize=(8, 8))
matrix_example = np.zeros((N, N))
np.fill_diagonal(matrix_example, 1)
for i in range(N-1):
    matrix_example[i, i+1] = 0.3  # 示例 α
    matrix_example[i+1, i] = -0.3
matrix_example[0, -1] = -0.3
matrix_example[-1, 0] = 0.3

im = ax4.imshow(matrix_example, cmap='coolwarm', vmin=-0.4, vmax=1)
ax4.set_title('矩阵 A 结构（循环三对角，角落项为周期边界）')
ax4.set_xlabel('列 j')
ax4.set_ylabel('行 i')
for i in range(N):
    for j in range(N):
        val = matrix_example[i,j]
        if val != 0:
            ax4.text(j, i, f'{val:.1f}' if abs(val)>0 else '1', ha='center', va='center', fontsize=10)

plt.colorbar(im, ax=ax4, label='系数值 (1 主对角, ±α 旁对角, 角落 ±α)')

plt.tight_layout()
plt.show()