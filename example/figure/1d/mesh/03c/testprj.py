import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(10, 4))

# 画十字线
plt.plot([-6, 6], [0, 0], 'r-', linewidth=2)
plt.plot([0, 0], [-2, 2], 'b-', linewidth=2)

plt.xlim(-6, 6)
plt.ylim(-2, 2)

# 计算并设置精确的宽高比
ax = plt.gca()

# 数据范围比例
data_ratio = (6 - (-6)) / (2 - (-2))  # x_range / y_range = 12/4 = 3

# 获取图形区域的宽高（归一化坐标）
pos = ax.get_position()
fig_width, fig_height = pos.width, pos.height

# 计算显示比例
display_ratio = fig_width / fig_height
print(f"display_ratio={display_ratio}")
print(f"data_ratio={data_ratio}")

# 设置宽高比：data_ratio / (display_ratio)
#ax.set_aspect(data_ratio / display_ratio)
ax.set_aspect(display_ratio / data_ratio)

plt.title(f'Manual aspect ratio: {display_ratio/data_ratio:.2f}')
plt.grid(True)
plt.show()