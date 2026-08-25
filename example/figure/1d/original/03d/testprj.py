import matplotlib.pyplot as plt
import numpy as np

# 设置图形大小和样式
plt.figure(figsize=(8, 1))  # 调整图形大小以匹配图片高度和宽度

# 绘制网格点 (红点)
x_points = np.array([-2, -1, 0, 1, 2, 3, 4, 5])  # 对应 i = -2 到 5
y_points = np.zeros_like(x_points)  # 所有点在 y=0 线上

plt.plot(x_points, y_points, 'ro', markersize=8)  # 红点

# 绘制黑点 (在 x=-1 和 x=3)
plt.plot([-1, 3], [0, 0], 'ko', markersize=8)  # 黑点

# 绘制 x=0 和 x=2 的虚线
plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
plt.axvline(x=2, color='k', linestyle='--', alpha=0.5)

# 绘制紫色矩形模板 (从 x=-1 到 x=2)
template = plt.Rectangle((-1, -0.5), 3, 1, color='purple', alpha=0.3)
plt.gca().add_patch(template)

# 设置轴
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 保存或显示图形
plt.savefig('reconstruction_template.png', bbox_inches='tight', dpi=300)
plt.show()