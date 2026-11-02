import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(12, 4))  # 关键：将画布设为长方形 (宽12, 高4)

# 创建数据点
x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6], dtype=np.float64)

# 画十字线
plt.plot([-6, 6], [0, 0], 'r-', linewidth=2, label='Horizontal line y=0')  # 水平线
plt.plot([0, 0], [-2, 2], 'b-', linewidth=2, label='Vertical line x=0')    # 垂直线

# 设置坐标范围
plt.xlim(-6, 6)
plt.ylim(-2, 2)

plt.title('Long strip grid: x-range [-6,6], y-range [-2,2]')
plt.legend()

# 添加网格线
plt.grid(True, linestyle='--', alpha=0.7)

# 显示坐标轴刻度
plt.xticks(np.arange(-6, 7, 1))
plt.yticks(np.arange(-2, 3, 0.5))

plt.show()