import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(9, 3))  # 宽:高 = 3:1，对应数据范围比例12:4=3:1

# 画你的数据点（如果需要）
x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6])
y_points = np.zeros_like(x_points)  # 假设y坐标都是0

plt.scatter(x_points, y_points, color='red', s=50, zorder=5, label='Data points')

# 画坐标轴
plt.axhline(y=0, color='blue', linewidth=2, label='x-axis')
plt.axvline(x=0, color='green', linewidth=2, label='y-axis')

# 设置范围
plt.xlim(-6.5, 6.5)  # 稍微扩大一点
plt.ylim(-2, 2)

plt.title('1D Grid Visualization')
plt.xlabel('X coordinate (wide range)')
plt.ylabel('Y coordinate (narrow range)')
plt.legend()
plt.grid(True, alpha=0.3)

# 关键：不要设置'equal'，保持默认的'auto'
# plt.axis('equal')  # ← 注释掉或删除这一行！

plt.tight_layout()
plt.show()