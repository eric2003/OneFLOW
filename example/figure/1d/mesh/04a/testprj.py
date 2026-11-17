import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(10, 6))

# 主边界竖线
main_line_x = -3
ax.axvline(x=main_line_x, color='darkblue', linewidth=4, alpha=0.8)

# 创建装饰性斜线模式
num_patterns = 15
y_range = np.linspace(-4, 4, num_patterns)

# 不同样式的斜线
patterns = [
    {'color': 'red', 'length': 0.8, 'angle': 30, 'style': '-'},
    {'color': 'green', 'length': 0.6, 'angle': 45, 'style': '--'},
    {'color': 'orange', 'length': 0.4, 'angle': 60, 'style': '-.'},
]

for i, y in enumerate(y_range):
    pattern = patterns[i % len(patterns)]  # 循环使用不同样式
    
    x_start = main_line_x
    y_start = y
    angle_rad = np.radians(pattern['angle'])
    
    x_end = x_start + pattern['length'] * np.cos(angle_rad)
    y_end = y_start + pattern['length'] * np.sin(angle_rad)
    
    ax.plot([x_start, x_end], [y_start, y_end],
            color=pattern['color'],
            linewidth=2,
            linestyle=pattern['style'],
            alpha=0.7)

# 添加一些点装饰
for y in np.linspace(-3.5, 3.5, 8):
    ax.scatter(main_line_x + 0.1, y, color='purple', s=50, alpha=0.6)

ax.set_xlim(-4, 2)
ax.set_ylim(-5, 5)
ax.set_title('Decorative boundary pattern')
ax.grid(True, alpha=0.2)
plt.show()