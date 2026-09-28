import matplotlib.pyplot as plt
import numpy as np

def draw_arrow_on_line(ax, x_start, y_start, x_end, y_end, position=0.5, 
                       arrow_style='->', color='blue', linewidth=2, 
                       head_size=15, zorder=2, show_connection=False):
    """
    Draw an arrow on a line segment at a specified position.
    
    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes to draw on
    x_start, y_start : float
        Starting point coordinates
    x_end, y_end : float
        Ending point coordinates
    position : float (default=0.5)
        Position along the line where to place the arrow (0=start, 1=end)
    arrow_style : str (default='->')
        Arrow style (e.g., '->', '-|>', '<-', '<->', 'fancy')
    color : str or tuple (default='blue')
        Arrow color
    linewidth : float (default=2)
        Arrow line width
    head_size : float (default=15)
        Arrow head size (mutation_scale parameter)
    zorder : int (default=2)
        Drawing order (higher values are drawn on top)
    show_connection : bool (default=False)
        Whether to show the connection line (arrow shaft)
    """
    # Calculate the coordinates for the arrow position
    arrow_x = x_start + position * (x_end - x_start)
    arrow_y = y_start + position * (y_end - y_start)
    
    # Calculate direction vector
    dx = x_end - x_start
    dy = y_end - y_start
    
    # Normalize direction vector
    length = np.sqrt(dx**2 + dy**2)
    if length > 0:
        dx_norm = dx / length
        dy_norm = dy / length
    else:
        dx_norm, dy_norm = 0, 0
    
    # Define arrow length
    arrow_length = 0.2 * length
    
    # Calculate arrow start and end points
    arrow_start_x = arrow_x - 0.4 * arrow_length * dx_norm
    arrow_start_y = arrow_y - 0.4 * arrow_length * dy_norm
    arrow_end_x = arrow_x + 0.4 * arrow_length * dx_norm
    arrow_end_y = arrow_y + 0.4 * arrow_length * dy_norm
    
    # Create arrow properties
    arrow_props = dict(arrowstyle=arrow_style,
                       color=color,
                       linewidth=linewidth,
                       mutation_scale=head_size,
                       shrinkA=0,
                       shrinkB=0)
    
    # If we don't want the connection line, set linestyle to 'none'
    if not show_connection:
        arrow_props['linestyle'] = 'none'
    
    # Draw the arrow
    arrow = ax.annotate('',
                        xy=(arrow_end_x, arrow_end_y),
                        xytext=(arrow_start_x, arrow_start_y),
                        arrowprops=arrow_props,
                        zorder=zorder)
    
    return arrow

# 方法2：使用箭头的偏移参数（更直接的方法）
def draw_arrow_only(ax, x_start, y_start, x_end, y_end, position=0.5,
                    arrow_style='->', color='blue', linewidth=2,
                    head_size=15, zorder=2):
    """
    Draw only the arrow head without the connecting line.
    """
    # Calculate arrow position
    arrow_x = x_start + position * (x_end - x_start)
    arrow_y = y_start + position * (y_end - y_start)
    
    # Calculate direction
    dx = x_end - x_start
    dy = y_end - y_start
    length = np.sqrt(dx**2 + dy**2)
    
    if length > 0:
        dx_norm = dx / length
        dy_norm = dy / length
    else:
        dx_norm, dy_norm = 0, 0
    
    # Very small offset to create just the arrow head
    offset = 0.001 * length
    
    # Create arrow
    arrow = ax.annotate('',
                        xy=(arrow_x + offset * dx_norm, arrow_y + offset * dy_norm),
                        xytext=(arrow_x - offset * dx_norm, arrow_y - offset * dy_norm),
                        arrowprops=dict(arrowstyle=arrow_style,
                                       color=color,
                                       linewidth=linewidth,
                                       mutation_scale=head_size,
                                       linestyle='none',  # No line!
                                       shrinkA=0,
                                       shrinkB=0),
                        zorder=zorder)
    
    return arrow

# 方法3：使用FancyArrowPatch直接控制
def draw_arrow_head_only(ax, x_start, y_start, x_end, y_end, position=0.5,
                         color='blue', linewidth=2, head_size=0.3, zorder=2):
    """
    Draw only the arrow head using FancyArrowPatch.
    """
    from matplotlib.patches import FancyArrowPatch
    
    # Calculate arrow position
    arrow_x = x_start + position * (x_end - x_start)
    arrow_y = y_start + position * (y_end - y_start)
    
    # Calculate direction
    dx = x_end - x_start
    dy = y_end - y_start
    length = np.sqrt(dx**2 + dy**2)
    
    if length > 0:
        dx_norm = dx / length
        dy_norm = dy / length
    else:
        dx_norm, dy_norm = 0, 0
    
    # Arrow length
    arrow_length = 0.2 * length
    
    # Create arrow with specified head properties
    arrow = FancyArrowPatch((arrow_x, arrow_y),
                           (arrow_x + arrow_length * dx_norm, 
                            arrow_y + arrow_length * dy_norm),
                           arrowstyle='->',
                           color=color,
                           linewidth=0,  # No line
                           mutation_scale=head_size * 20,  # Scale factor
                           zorder=zorder)
    
    ax.add_patch(arrow)
    return arrow

# 示例：对比不同方法
def compare_methods():
    """对比显示连接线和只显示箭头的方法"""
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    # 定义线段
    x_start, y_start = 0, 0
    x_end, y_end = 1, 1
    
    methods = [
        ("With connection line", True, draw_arrow_on_line),
        ("Without connection line", False, draw_arrow_on_line),
        ("draw_arrow_only method", None, draw_arrow_only),
        ("draw_arrow_head_only method", None, draw_arrow_head_only),
        ("Different arrow styles", None, None),
        ("Multiple arrows", None, None),
    ]
    
    for i, (title, show_conn, method_func) in enumerate(methods):
        ax = axes[i]
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.set_title(title, fontsize=12)
        
        # 绘制基准线段
        ax.plot([x_start, x_end], [y_start, y_end], 
                'gray', linewidth=2, alpha=0.2)
        
        if i == 0 or i == 1:
            # 方法1：使用linestyle控制
            draw_arrow_on_line(ax, x_start, y_start, x_end, y_end,
                              position=0.5,
                              color='blue',
                              show_connection=show_conn,
                              head_size=20)
            
        elif i == 2:
            # 方法2：draw_arrow_only
            draw_arrow_only(ax, x_start, y_start, x_end, y_end,
                           position=0.5,
                           color='red',
                           head_size=20)
            
        elif i == 3:
            # 方法3：draw_arrow_head_only
            draw_arrow_head_only(ax, x_start, y_start, x_end, y_end,
                                position=0.5,
                                color='green',
                                head_size=0.4)
            
        elif i == 4:
            # 不同箭头样式
            arrow_styles = ['->', '-|>', '<-', '<->', 'fancy']
            colors = ['red', 'blue', 'green', 'purple', 'orange']
            
            for j, (style, color) in enumerate(zip(arrow_styles, colors)):
                position = 0.2 + j * 0.15
                draw_arrow_only(ax, x_start, y_start, x_end, y_end,
                               position=position,
                               arrow_style=style,
                               color=color,
                               head_size=15)
                ax.text(position, 0.95, style, ha='center', fontsize=8)
                
        elif i == 5:
            # 多个箭头
            positions = [0.2, 0.4, 0.6, 0.8]
            colors = plt.cm.viridis(np.linspace(0, 1, len(positions)))
            
            for pos, color in zip(positions, colors):
                draw_arrow_only(ax, x_start, y_start, x_end, y_end,
                               position=pos,
                               color=color,
                               head_size=15)
    
    plt.tight_layout()
    plt.show()

# 实用示例：在折线上绘制箭头
def polyline_with_arrows():
    """在折线上绘制只显示箭头的例子"""
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 定义折线点
    points = [(0, 0), (2, 3), (4, 1), (6, 4), (8, 2)]
    x_coords = [p[0] for p in points]
    y_coords = [p[1] for p in points]
    
    # 绘制折线
    ax.plot(x_coords, y_coords, 'gray', linewidth=3, alpha=0.3, 
            marker='o', markersize=8, label='Polyline')
    
    # 在每个线段上绘制箭头（只显示箭头头）
    for i in range(len(points)-1):
        x1, y1 = points[i]
        x2, y2 = points[i+1]
        
        # 在线段的不同位置绘制箭头
        for position in [0.25, 0.5, 0.75]:
            draw_arrow_only(ax, x1, y1, x2, y2,
                           position=position,
                           arrow_style='-|>',
                           color=f'C{i}',
                           linewidth=1.5,
                           head_size=12)
    
    # 设置图形属性
    ax.set_xlim(-0.5, 8.5)
    ax.set_ylim(-0.5, 4.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title("Arrow Heads on Polyline (No Connection Lines)")
    ax.legend()
    
    plt.tight_layout()
    plt.show()

# 向量场示例
def vector_field_example():
    """向量场风格的箭头示例"""
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 创建网格
    x = np.linspace(0, 4, 5)
    y = np.linspace(0, 4, 5)
    X, Y = np.meshgrid(x, y)
    
    # 向量场函数（旋转场）
    U = -Y + 2  # X分量
    V = X - 2   # Y分量
    
    # 绘制向量场（只显示箭头头）
    for i in range(len(x)):
        for j in range(len(y)):
            # 计算向量长度
            length = np.sqrt(U[j, i]**2 + V[j, i]**2)
            
            if length > 0:
                # 归一化
                u_norm = U[j, i] / length
                v_norm = V[j, i] / length
                
                # 绘制箭头（只显示箭头头）
                draw_arrow_only(ax, 
                               x_start=X[j, i] - 0.1 * u_norm,
                               y_start=Y[j, i] - 0.1 * v_norm,
                               x_end=X[j, i] + 0.1 * u_norm,
                               y_end=Y[j, i] + 0.1 * v_norm,
                               position=0.5,
                               arrow_style='-|>',
                               color='blue',
                               head_size=10 + 5 * length)  # 长度越大箭头越大
    
    # 设置图形属性
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-0.5, 4.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.set_title("Vector Field with Arrow Heads Only")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("对比不同方法：")
    compare_methods()
    
    print("\n折线箭头示例：")
    polyline_with_arrows()
    
    print("\n向量场示例：")
    vector_field_example()