import matplotlib.pyplot as plt
import numpy as np

# 设置图形大小和样式
plt.figure(figsize=(12, 1))  # 调整图形大小以适应新的点分布

# 定义新的点坐标 (从 i=-6 到 i=6，基于 new_x_points)
new_x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 4, 5, 6])
y_points = np.zeros_like(new_x_points)  # 所有点在 y=0 线上

# 绘制除中间 5 个点和特定边缘点外的其他点 (内部红色，边缘黑色)
# 中间 5 个点对应 new_x_points[3:8] (i=-2, -1, 0, 1, 2)
# 其他点需手动指定，排除左侧 i=-4 和右侧 i=4
edge_points_red = np.concatenate([new_x_points[:2], new_x_points[2:3], new_x_points[8:9], new_x_points[9:]])  # 对应 i=-6, -5, -4, 4, 5, 6，但 i=-4 和 i=4 单独处理
plt.scatter(edge_points_red, np.zeros_like(edge_points_red), s=100, facecolor='red', edgecolor='black', linewidth=1)  # 红内黑边点

# 绘制左侧第三点 (i=-4) 和右侧第三点 (i=4) 为纯黑色点 (内部黑色，边缘黑色)
special_black_points = np.array([-4, 4])  # 对应新位置 i=-4, 4
plt.scatter(special_black_points, np.zeros_like(special_black_points), s=100, facecolor='black', edgecolor='black', linewidth=1)  # 纯黑色点

# 绘制中间 5 个点 (i=-2, -1, 0, 1, 2 对应新位置 -2, -1, 0, 1, 2) (内部黑色，边缘黑色，即纯黑色点)
middle_points = new_x_points[3:8]  # 对应新位置 i=-2, -1, 0, 1, 2
plt.scatter(middle_points, np.zeros_like(middle_points), s=100, facecolor='black', edgecolor='black', linewidth=1)  # 纯黑色点

# 绘制中间 5 个点的黑实线连接 (i=-2, -1, 0, 1, 2 对应新位置 -2, -1, 0, 1, 2)
plt.plot(middle_points, np.zeros_like(middle_points), 'k-', linewidth=1)  # 黑实线

# 去掉左侧第二点和第三点的连线 (i=-5 到 i=-4)，以及右侧第二点和第三点的连线 (i=5 到 i=4)
# 原代码中的混合线已移除，不再绘制：
# - 左边 (i=-5 到 i=-4)
# - 右边 (i=5 到 i=4)

# 添加所有点的标签，与点保持一定距离，高度一致
plt.text(-6, -0.5, '$-2$', fontsize=12, ha='center')  # 在 i=-6 下方 0.5 单位
plt.text(-5, -0.5, '$-1$', fontsize=12, ha='center')  # 在 i=-5 下方 0.5 单位
plt.text(-4, -0.5, '$i=1$', fontsize=12, ha='center')  # 在 i=-4 下方 0.5 单位
plt.text(-2, -0.5, '$i-2$', fontsize=12, ha='center')  # 在 i=-2 下方 0.5 单位
plt.text(-1, -0.5, '$i-1$', fontsize=12, ha='center')  # 在 i=-1 下方 0.5 单位
plt.text(0, -0.5, '$i$', fontsize=12, ha='center')  # 在 i=0 下方 0.5 单位
plt.text(1, -0.5, '$i+1$', fontsize=12, ha='center')  # 在 i=1 下方 0.5 单位
plt.text(2, -0.5, '$i+2$', fontsize=12, ha='center')  # 在 i=2 下方 0.5 单位
plt.text(4, -0.5, '$i=N+1$', fontsize=12, ha='center')  # 在 i=4 下方 0.5 单位
plt.text(5, -0.5, '$N+2$', fontsize=12, ha='center')  # 在 i=5 下方 0.5 单位
plt.text(6, -0.5, '$N+3$', fontsize=12, ha='center')  # 在 i=6 下方 0.5 单位

# 设置轴
plt.axis('equal')  # 保持比例
plt.axis('off')  # 隐藏轴

# 保存或显示图形
plt.savefig('eleven_points_adjusted_edges.png', bbox_inches='tight', dpi=300)
plt.show()