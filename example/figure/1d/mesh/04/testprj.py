import matplotlib.pyplot as plt
import numpy as np

# 创建一个简单的边界面
fig, ax = plt.subplots(figsize=(8, 6))

# 主竖线
x_main = 0
ax.axvline(x=x_main, color='black', linewidth=3, label='Main vertical line')

# 在主竖线旁边添加斜线
x_offset = 0.2  # 斜线离主线的距离
num_slashes = 20  # 斜线数量
y_positions = np.linspace(-4, 4, num_slashes)

for i, y in enumerate(y_positions):
    # 短斜线的起点和终点
    length = 0.5
    angle = 45  # 斜线角度
    
    # 计算斜线端点
    x_start = x_main
    y_start = y
    x_end = x_start + length * np.cos(np.radians(angle))
    y_end = y_start + length * np.sin(np.radians(angle))
    
    ax.plot([x_start, x_end], [y_start, y_end], 
            color='red', linewidth=1.5, alpha=0.7)

ax.set_xlim(-2, 2)
ax.set_ylim(-5, 5)
ax.set_aspect('equal')
ax.set_title('Basic vertical line with diagonal slashes')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()