import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# 设置字体为 Times New Roman
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

# 设置图形大小和样式
plt.figure(figsize=(12, 3))

# 定义新的点坐标 (从 i=-6 到 i=6，基于 new_x_points)
new_x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 4, 5, 6])
y_points = np.zeros_like(new_x_points)  # 所有点在 y=0 线上

# 绘制除中间 5 个点和特定边缘点外的其他点 (内部红色，边缘黑色)
edge_points_red = np.concatenate([new_x_points[:2], new_x_points[2:3], new_x_points[8:9], new_x_points[9:]])
plt.scatter(edge_points_red, np.zeros_like(edge_points_red), s=100, facecolor='red', edgecolor='black', linewidth=1)

# 绘制左侧第三点 (i=-4) 和右侧第三点 (i=4) 为纯黑色点
special_black_points = np.array([-4, 4])
plt.scatter(special_black_points, np.zeros_like(special_black_points), s=100, facecolor='black', edgecolor='black', linewidth=1)

# 绘制中间 5 个点 (i=-2, -1, 0, 1, 2)
middle_points = new_x_points[3:8]
plt.scatter(middle_points, np.zeros_like(middle_points), s=100, facecolor='black', edgecolor='black', linewidth=1)

# 绘制中间 5 个点的黑实线连接
plt.plot(middle_points, np.zeros_like(middle_points), 'k-', linewidth=1)

# 添加左起第三点和第四点之间的分段连线（-4到-2）
x_left = np.array([-4, -3.5, -2.5, -2])  # 按照1/4, 1/2, 1/4的比例分割
plt.plot([x_left[0], x_left[1]], [0, 0], 'k-', linewidth=1)  # 第一段实线
plt.plot([x_left[1], x_left[2]], [0, 0], 'k--', linewidth=1) # 第二段虚线
plt.plot([x_left[2], x_left[3]], [0, 0], 'k-', linewidth=1)  # 第三段实线

# 添加右起第三点和第四点之间的分段连线（2到4）
x_right = np.array([2, 2.5, 3.5, 4])    # 按照1/4, 1/2, 1/4的比例分割
plt.plot([x_right[0], x_right[1]], [0, 0], 'k-', linewidth=1)  # 第一段实线
plt.plot([x_right[1], x_right[2]], [0, 0], 'k--', linewidth=1) # 第二段虚线
plt.plot([x_right[2], x_right[3]], [0, 0], 'k-', linewidth=1)  # 第三段实线

xv = 0.5*(new_x_points[5]+new_x_points[6])
plt.plot([xv, xv], [-0.5, 0.5], 'k--')  # 绘制垂直线

# 添加圆角矩形背景
rectangle_color = (150/255, 150/255, 200/255)  # RGB颜色
#rectangle_color = 'lightblue'

x = -2.2
y = -0.1
width = 4.4
height = 0.2

rect = patches.FancyBboxPatch((x, y), width, height, 
                                       boxstyle="round,pad=0.1,rounding_size=0.1",
                                       edgecolor='none',
                                       facecolor=rectangle_color,
                                       zorder=0)

# 将矩形添加到坐标轴上
ax = plt.gca()
ax.add_patch(rect)
ax.set_title('Left-side reconstruction')

# 添加标签和其他设置（保持不变）
plt.text(-6, -0.5, '$-2$', fontsize=12, ha='center')
plt.text(-5, -0.5, '$-1$', fontsize=12, ha='center')
plt.text(-4, -0.5, '$i=1$', fontsize=12, ha='center')
plt.text(-2, -0.5, '$i-2$', fontsize=12, ha='center')
plt.text(-1, -0.5, '$i-1$', fontsize=12, ha='center')
plt.text(0, -0.5, '$i$', fontsize=12, ha='center')
plt.text(1, -0.5, '$i+1$', fontsize=12, ha='center')
plt.text(2, -0.5, '$i+2$', fontsize=12, ha='center')
plt.text(4, -0.5, '$i=N+1$', fontsize=12, ha='center')
plt.text(5, -0.5, '$N+2$', fontsize=12, ha='center')
plt.text(6, -0.5, '$N+3$', fontsize=12, ha='center')

plt.text(xv-0.3, 0.5, r'$u_{i+\frac{1}{2}}^L$', fontsize=12, ha='center')

plt.axis('equal')
plt.axis('off')

plt.savefig('eleven_points_adjusted_edges_with_rectangle.png', bbox_inches='tight', dpi=300)
plt.show()