import matplotlib.pyplot as plt
import numpy as np

def draw_boundary_pattern(ax, x_position, y_range=(-5, 5), 
                         num_slashes=30, slash_length=0.5,
                         slash_angle=45, color='blue',
                         line_width=2, main_line_width=3):
    """
    绘制一个边界面图案
    
    参数:
    ax: matplotlib坐标轴对象
    x_position: 主竖线的x坐标
    y_range: 斜线的y坐标范围 (min, max)
    num_slashes: 斜线数量
    slash_length: 斜线长度
    slash_angle: 斜线角度（度）
    color: 颜色
    """
    
    # 画主竖线
    ax.axvline(x=x_position, color=color, linewidth=main_line_width, alpha=0.8)
    
    # 生成斜线
    y_positions = np.linspace(y_range[0], y_range[1], num_slashes)
    angle_rad = np.radians(slash_angle)
    
    for i, y in enumerate(y_positions):
        # 交替改变斜线方向
        direction = 1 if i % 2 == 0 else -1
        current_angle = slash_angle * direction
        
        # 计算斜线端点
        x_start = x_position
        y_start = y
        x_end = x_start + slash_length * np.cos(np.radians(current_angle))
        y_end = y_start + slash_length * np.sin(np.radians(current_angle))
        
        # 绘制斜线
        ax.plot([x_start, x_end], [y_start, y_end],
                color=color, linewidth=line_width, alpha=0.6)
        
        # 在斜线末端添加小圆点
        ax.scatter(x_end, y_end, color=color, s=20, alpha=0.8, zorder=5)

# 使用自定义函数
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

# 不同的边界样式
configs = [
    {'x_position': -2, 'color': 'red', 'slash_angle': 30, 'num_slashes': 15},
    {'x_position': 0, 'color': 'green', 'slash_angle': 60, 'num_slashes': 20},
    {'x_position': 2, 'color': 'blue', 'slash_angle': 45, 'num_slashes': 25},
    {'x_position': -1, 'color': 'purple', 'slash_angle': 135, 'num_slashes': 18},
]

for ax, config in zip(axes, configs):
    draw_boundary_pattern(ax, **config, y_range=(-4, 4))
    ax.set_xlim(-3, 3)
    ax.set_ylim(-5, 5)
    ax.set_aspect('equal')
    ax.set_title(f"Boundary at x={config['x_position']}")
    ax.grid(True, alpha=0.2)

plt.tight_layout()
plt.show()