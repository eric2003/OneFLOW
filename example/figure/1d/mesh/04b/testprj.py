import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(10, 8))

# 两条主竖线
left_line = -2
right_line = 2
ax.axvline(x=left_line, color='black', linewidth=3, label='Left boundary')
ax.axvline(x=right_line, color='black', linewidth=3, label='Right boundary')

# 在两条竖线之间填充内容
x_fill = np.linspace(left_line, right_line, 100)
y_fill = np.sin(x_fill * 2) * 2
ax.fill_between(x_fill, -2, y_fill, alpha=0.2, color='skyblue', label='Area')

# 在边界上添加对称的斜线
num_slashes = 25
y_positions = np.linspace(-3, 3, num_slashes)

for y in y_positions:
    # 左边界斜线（向左）
    length = 0.6
    angle = 135  # 指向左下方
    
    # 左边界
    x_start = left_line
    y_start = y
    x_end = x_start + length * np.cos(np.radians(angle))
    y_end = y_start + length * np.sin(np.radians(angle))
    ax.plot([x_start, x_end], [y_start, y_end], 'r-', linewidth=1.5, alpha=0.6)
    
    # 右边界斜线（向右）
    angle = 45  # 指向右下方
    
    x_start = right_line
    y_start = y
    x_end = x_start + length * np.cos(np.radians(angle))
    y_end = y_start + length * np.sin(np.radians(angle))
    ax.plot([x_start, x_end], [y_start, y_end], 'r-', linewidth=1.5, alpha=0.6)

ax.set_xlim(-3, 3)
ax.set_ylim(-4, 4)
ax.set_aspect('equal')
ax.set_title('Boundary lines with symmetrical slashes')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()