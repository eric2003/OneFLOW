import matplotlib.pyplot as plt
import matplotlib.patches as patches

# 创建一个新的图形
fig, ax = plt.subplots()

# 添加一条竖向线段
line = ax.plot([0, 0], [0, 1], color='blue')[0]

# 在线段末端添加箭头
arrow = patches.FancyArrowPatch(
    (0, 1),  # 箭头起点坐标
    (0, 0),  # 箭头终点坐标
    arrowstyle='->',  # 箭头样式
    mutation_scale=20  # 箭头大小
)

# 将箭头添加到图形中
ax.add_patch(arrow)

# 设置坐标轴的范围
ax.set_xlim(-0.1, 0.1)
ax.set_ylim(-0.1, 1.1)

# 隐藏坐标轴
ax.axis('off')

# 显示图形
plt.show()