import matplotlib.pyplot as plt
import numpy as np

# 设置图形大小和样式
plt.figure(figsize=(12, 1))  # 调整图形大小以适应新的点分布

# 定义新的点坐标 (从 i=-6 到 i=6，基于 new_x_points)
new_x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 4, 5, 6])
y_points = np.zeros_like(new_x_points)  # 所有点在 y=0 线上

# 绘制除中间 5 个点外的其他点 (内部红色，边缘黑色)
# 中间 5 个点对应 new_x_points[3:8] (i=-2, -1, 0, 1, 2)
# 其他点对应 new_x_points[:3] 和 new_x_points[8:] (i=-6, -5, -4, 4, 5, 6)
edge_points = np.concatenate([new_x_points[:3], new_x_points[8:]])  # 对应 i=-6, -5, -4, 4, 5, 6
plt.scatter(edge_points, np.zeros_like(edge_points), s=100, facecolor='red', edgecolor='black', linewidth=1)  # 红内黑边点

# 绘制中间 5 个点 (i=-2, -1, 0, 1, 2 对应新位置 -2, -1, 0, 1, 2) (内部黑色，边缘黑色，即纯黑色点)
middle_points = new_x_points[3:8]  # 对应新位置 i=-2, -1, 0, 1, 2
plt.scatter(middle_points, np.zeros_like(middle_points), s=100, facecolor='black', edgecolor='black', linewidth=1)  # 纯黑色点

# 绘制中间 5 个点的黑实线连接 (i=-2, -1, 0, 1, 2 对应新位置 -2, -1, 0, 1, 2)
plt.plot(middle_points, np.zeros_like(middle_points), 'k-', linewidth=1)  # 黑实线

# 绘制边缘第 3 个点与第 4 个点之间的混合线
# 左边 (i=-5 到 i=-4，新位置 -5 到 -4)
left_edge_points = new_x_points[1:3]  # 对应新位置 i=-5, -4
# 绘制一半实线 (从 i=-5 到中间)，一半虚线 (从中间到 i=-4)
x_half = np.linspace(-5, -4.5, 10)  # 实线部分
plt.plot(x_half, np.zeros_like(x_half), 'k-', linewidth=1)  # 黑实线
x_half_dash = np.linspace(-4.5, -4, 10)  # 虚线部分
plt.plot(x_half_dash, np.zeros_like(x_half_dash), 'k--', linewidth=1)  # 黑虚线

# 右边 (i=5 到 i=4，新位置 5 到 4)
right_edge_points = new_x_points[8:10]  # 对应新位置 i=5, 4
# 绘制一半实线 (从 i=5 到中间)，一半虚线 (从中间到 i=4)
x_half_right = np.linspace(5, 4.5, 10)  # 实线部分
plt.plot(x_half_right, np.zeros_like(x_half_right), 'k-', linewidth=1)  # 黑实线
x_half_dash_right = np.linspace(4.5, 4, 10)  # 虚线部分
plt.plot(x_half_dash_right, np.zeros_like(x_half_dash_right), 'k--', linewidth=1)  # 黑虚线

# 添加左侧点 (-6, -5) 下方的标签 "-2" 和 "-1"，与点保持一定距离
plt.text(-6, -0.5, '-2', fontsize=12, ha='center')  # 在 i=-6 下方 0.5 单位
plt.text(-5, -0.5, '-1', fontsize=12, ha='center')  # 在 i=-5 下方 0.5 单位

# 添加右侧点 (5, 6) 下方的标签 "N+2" 和 "N+3"，与点保持一定距离，并与左侧标签高度一致
plt.text(5, -0.5, '$N+2$', fontsize=12, ha='center')  # 在 i=5 下方 0.5 单位，使用 LaTeX 格式
plt.text(6, -0.5, '$N+3$', fontsize=12, ha='center')  # 在 i=6 下方 0.5 单位，使用 LaTeX 格式

# 设置轴
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 保存或显示图形
plt.savefig('eleven_points_new_positions.png', bbox_inches='tight', dpi=300)
plt.show()