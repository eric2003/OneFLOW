import matplotlib.pyplot as plt
import numpy as np

# 设置图形样式
plt.rcParams['font.size'] = 12

# 绘制 (a) Left-side reconstruction
plt.figure(figsize=(8, 2))  # 调整图形大小以匹配图片宽度

# 绘制网格点 (红点)
x_points = np.array([-2, -1, 0, 1, 2, 3, 4, 5, 6])  # 对应 i = -2 到 N+3
y_points = np.zeros_like(x_points)  # 所有点在 y=0 线上

plt.plot(x_points, y_points, 'ro', markersize=8)  # 红点

# 标记特定点 (i 值)
for i, x in enumerate(x_points):
    plt.text(x, 0.1, f'$i={x}$', fontsize=12, ha='center')

# 绘制 x=0 和 x=L 的虚线
plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
plt.axvline(x=6, color='k', linestyle='--', alpha=0.5)
plt.text(0, -0.2, '$x=0$', fontsize=12, ha='center')
plt.text(6, -0.2, '$x=L$', fontsize=12, ha='center')

# 绘制紫色矩形模板 (从 i-1 到 i+2, 覆盖 i+1)
template_left = plt.Rectangle((-0.5, -0.5), 3.5, 1, color='purple', alpha=0.3)
plt.gca().add_patch(template_left)
plt.text(1, 0.3, '$u_{i+1}^L$', fontsize=12, ha='center')

# 设置标题和轴
plt.title('(a) Left-side reconstruction', fontsize=14, pad=10)
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 保存或显示图形
plt.savefig('left_side_reconstruction.png', bbox_inches='tight', dpi=300)
plt.show()

# 绘制 (b) Right-side reconstruction
plt.figure(figsize=(8, 2))  # 调整图形大小以匹配图片宽度

# 绘制网格点 (红点)
plt.plot(x_points, y_points, 'ro', markersize=8)  # 红点

# 标记特定点 (i 值)
for i, x in enumerate(x_points):
    plt.text(x, 0.1, f'$i={x}$', fontsize=12, ha='center')

# 绘制 x=0 和 x=L 的虚线
plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
plt.axvline(x=6, color='k', linestyle='--', alpha=0.5)
plt.text(0, -0.2, '$x=0$', fontsize=12, ha='center')
plt.text(6, -0.2, '$x=L$', fontsize=12, ha='center')

# 绘制紫色矩形模板 (从 i-1 到 i+2, 覆盖 i+1)
template_right = plt.Rectangle((-0.5, -0.5), 3.5, 1, color='purple', alpha=0.3)
plt.gca().add_patch(template_right)
plt.text(1, 0.3, '$u_{i+1}^R$', fontsize=12, ha='center')

# 设置标题和轴
plt.title('(b) Right-side reconstruction', fontsize=14, pad=10)
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 保存或显示图形
plt.savefig('right_side_reconstruction.png', bbox_inches='tight', dpi=300)
plt.show()