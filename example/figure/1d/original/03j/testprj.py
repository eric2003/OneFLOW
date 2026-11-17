import matplotlib.pyplot as plt
import numpy as np

# 设置图形大小和样式
plt.figure(figsize=(10, 1))  # 调整图形大小以适应 11 个点

# 定义 11 个点的坐标 (从 i=-5 到 i=5)
x_points = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
y_points = np.zeros_like(x_points)  # 所有点在 y=0 线上

# 绘制除中间 5 个点外的其他点 (内部红色，边缘黑色)
edge_points = np.concatenate([x_points[:3], x_points[8:]])  # 对应 i=-5, -4, -3, 3, 4, 5
plt.scatter(edge_points, np.zeros_like(edge_points), s=100, facecolor='red', edgecolor='black', linewidth=1)  # 红内黑边点

# 绘制中间 5 个点 (i=-2, -1, 0, 1, 2) (内部黑色，边缘黑色，即纯黑色点)
middle_points = x_points[3:8]  # 对应 i=-2, -1, 0, 1, 2
plt.scatter(middle_points, np.zeros_like(middle_points), s=100, facecolor='black', edgecolor='black', linewidth=1)  # 纯黑色点

# 绘制中间 5 个点的黑实线连接 (i=-2, -1, 0, 1, 2)
plt.plot(middle_points, np.zeros_like(middle_points), 'k-', linewidth=1)  # 黑实线

# 绘制边缘第 3 个点与第 4 个点之间的混合线 (i=-4 到 i=-3, 和 i=4 到 i=3)
# 左边 (i=-4 到 i=-3)
left_edge_points = x_points[1:3]  # 对应 i=-4, -3
# 绘制一半实线 (从 i=-4 到中间)，一半虚线 (从中间到 i=-3)
x_half = np.linspace(-4, -3.5, 10)  # 实线部分
plt.plot(x_half, np.zeros_like(x_half), 'k-', linewidth=1)  # 黑实线
x_half_dash = np.linspace(-3.5, -3, 10)  # 虚线部分
plt.plot(x_half_dash, np.zeros_like(x_half_dash), 'k--', linewidth=1)  # 黑虚线

# 右边 (i=4 到 i=3)
right_edge_points = x_points[8:10]  # 对应 i=4, 3
# 绘制一半实线 (从 i=4 到中间)，一半虚线 (从中间到 i=3)
x_half_right = np.linspace(4, 3.5, 10)  # 实线部分
plt.plot(x_half_right, np.zeros_like(x_half_right), 'k-', linewidth=1)  # 黑实线
x_half_dash_right = np.linspace(3.5, 3, 10)  # 虚线部分
plt.plot(x_half_dash_right, np.zeros_like(x_half_dash_right), 'k--', linewidth=1)  # 黑虚线

# 添加左侧点 (-5, -4) 下方的标签 "-2" 和 "-1"，与点保持一定距离
plt.text(-5, -0.5, '-2', fontsize=12, ha='center')  # 在 i=-5 下方 0.5 单位
plt.text(-4, -0.5, '-1', fontsize=12, ha='center')  # 在 i=-4 下方 0.5 单位

# 设置轴
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 保存或显示图形
plt.savefig('eleven_points_middle_black_edge_black.png', bbox_inches='tight', dpi=300)
plt.show()