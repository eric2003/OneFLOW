import matplotlib.pyplot as plt
import numpy as np

# 设置图形大小和样式
plt.figure(figsize=(12, 1))  # 调整图形大小以适应新的点分布

# 定义原始 11 个点的坐标 (从 i=-5 到 i=5)
original_x_points = np.array([-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])
y_points = np.zeros_like(original_x_points)  # 所有点在 y=0 线上

# 计算新的点坐标：仅左边的第 3 个点和第 4 个点 (i=-3, -2) 之间的距离变为原距离的两倍
# 原距离为 1 (从 -3 到 -2)，新距离为 2
# 仅右边的第 3 个点和第 4 个点 (i=3, 2) 之间的距离变为原距离的两倍
new_x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 4, 5, 6])
# 解释：
# - 左边：i=-5, -4 保持不变 (-5, -4)
# - 左边的第 3 个点 i=-3 保持在 -3
# - 左边的第 4 个点 i=-2 移动到 -1 (距离从 -3 到 -1 为 2)
# - i=-1, 0, 1 保持不变 (-1, 0, 1)
# - 右边的第 3 个点 i=3 保持在 3
# - 右边的第 4 个点 i=2 移动到 1 (距离从 3 到 1 为 2)
# - 右边：i=4, 5 保持不变 (4, 5)

# 绘制除中间 5 个点外的其他点 (内部红色，边缘黑色)
edge_points = np.concatenate([new_x_points[:3], new_x_points[8:]])  # 对应 i=-5, -4, -3, 3, 4, 5
plt.scatter(edge_points, np.zeros_like(edge_points), s=100, facecolor='red', edgecolor='black', linewidth=1)  # 红内黑边点

# 绘制中间 5 个点 (i=-2, -1, 0, 1, 2 对应新位置 -1, -1, 0, 1, 1) (内部黑色，边缘黑色，即纯黑色点)
middle_points = new_x_points[3:8]  # 对应新位置 i=-1, -1, 0, 1, 1
plt.scatter(middle_points, np.zeros_like(middle_points), s=100, facecolor='black', edgecolor='black', linewidth=1)  # 纯黑色点

# 绘制中间 5 个点的黑实线连接 (i=-2, -1, 0, 1, 2 对应新位置 -1, -1, 0, 1, 1)
plt.plot(middle_points, np.zeros_like(middle_points), 'k-', linewidth=1)  # 黑实线

# 绘制边缘第 3 个点与第 4 个点之间的混合线
# 左边 (i=-4 到 i=-3，新位置 -4 到 -3)
left_edge_points = new_x_points[1:3]  # 对应新位置 i=-4, -3
# 绘制一半实线 (从 i=-4 到中间)，一半虚线 (从中间到 i=-3)
x_half = np.linspace(-4, -3.5, 10)  # 实线部分
plt.plot(x_half, np.zeros_like(x_half), 'k-', linewidth=1)  # 黑实线
x_half_dash = np.linspace(-3.5, -3, 10)  # 虚线部分
plt.plot(x_half_dash, np.zeros_like(x_half_dash), 'k--', linewidth=1)  # 黑虚线

# 右边 (i=4 到 i=3，新位置 4 到 3)
right_edge_points = new_x_points[8:10]  # 对应新位置 i=4, 3
# 绘制一半实线 (从 i=4 到中间)，一半虚线 (从中间到 i=3)
x_half_right = np.linspace(4, 3.5, 10)  # 实线部分
plt.plot(x_half_right, np.zeros_like(x_half_right), 'k-', linewidth=1)  # 黑实线
x_half_dash_right = np.linspace(3.5, 3, 10)  # 虚线部分
plt.plot(x_half_dash_right, np.zeros_like(x_half_dash_right), 'k--', linewidth=1)  # 黑虚线

# 添加左侧点 (-5, -4) 下方的标签 "-2" 和 "-1"，与点保持一定距离
plt.text(-5, -0.5, '-2', fontsize=12, ha='center')  # 在 i=-5 下方 0.5 单位
plt.text(-4, -0.5, '-1', fontsize=12, ha='center')  # 在 i=-4 下方 0.5 单位

# 添加右侧点 (4, 5) 下方的标签 "N+2" 和 "N+3"，与点保持一定距离，并与左侧标签高度一致
plt.text(4, -0.5, '$N+2$', fontsize=12, ha='center')  # 在 i=4 下方 0.5 单位，使用 LaTeX 格式
plt.text(5, -0.5, '$N+3$', fontsize=12, ha='center')  # 在 i=5 下方 0.5 单位，使用 LaTeX 格式

# 设置轴
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 保存或显示图形
plt.savefig('eleven_points_adjusted_third_fourth_distance.png', bbox_inches='tight', dpi=300)
plt.show()