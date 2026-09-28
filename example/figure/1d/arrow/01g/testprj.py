import matplotlib.pyplot as plt
import numpy as np

def draw_arrow_on_line(ax, x_start, y_start, x_end, y_end, position=0.5, 
                       arrow_style='->', color='blue', linewidth=2, 
                       head_size=15, label=None, zorder=2, show_marker=False):
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
    label : str or None (default=None)
        Label for the arrow (for legend)
    zorder : int (default=2)
        Drawing order (higher values are drawn on top)
    show_marker : bool (default=False)
        Whether to show a marker at the arrow position
    
    Returns:
    --------
    arrow_annotation : matplotlib.text.Annotation
        The arrow annotation object
    """
    # Validate position parameter
    if position < 0 or position > 1:
        raise ValueError("Position must be between 0 and 1")
    
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
    
    # Define arrow length (relative to line length)
    arrow_length = 0.2 * length  # 20% of line length
    
    # Calculate arrow start and end points
    # Arrow points in the direction of the line
    arrow_start_x = arrow_x - 0.4 * arrow_length * dx_norm
    arrow_start_y = arrow_y - 0.4 * arrow_length * dy_norm
    arrow_end_x = arrow_x + 0.4 * arrow_length * dx_norm
    arrow_end_y = arrow_y + 0.4 * arrow_length * dy_norm
    
    # Draw the arrow
    arrow = ax.annotate('',
                        xy=(arrow_end_x, arrow_end_y),
                        xytext=(arrow_start_x, arrow_start_y),
                        arrowprops=dict(arrowstyle=arrow_style,
                                       color=color,
                                       linewidth=linewidth,
                                       mutation_scale=head_size,
                                       shrinkA=0,
                                       shrinkB=0),
                        zorder=zorder)
    
    # Add a small marker at the arrow position only if requested
    marker = None
    if show_marker:
        marker = ax.scatter(arrow_x, arrow_y, color=color, s=30, 
                           zorder=zorder+1, alpha=0.6)
    
    # Add label if provided
    if label:
        ax.text(arrow_x, arrow_y + 0.05*length, label, 
                color=color, fontsize=9, ha='center', va='bottom',
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
    
    return arrow, marker

# 简洁版本：完全不绘制标记点
def draw_arrow_on_line_simple(ax, x_start, y_start, x_end, y_end, position=0.5, 
                              arrow_style='->', color='blue', linewidth=2, 
                              head_size=15, zorder=2):
    """
    Simplified version - draw arrow without marker or label.
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
    
    # Arrow length
    arrow_length = 0.2 * length
    
    # Arrow start and end
    arrow_start_x = arrow_x - 0.4 * arrow_length * dx_norm
    arrow_start_y = arrow_y - 0.4 * arrow_length * dy_norm
    arrow_end_x = arrow_x + 0.4 * arrow_length * dx_norm
    arrow_end_y = arrow_y + 0.4 * arrow_length * dy_norm
    
    # Draw arrow
    arrow = ax.annotate('',
                        xy=(arrow_end_x, arrow_end_y),
                        xytext=(arrow_start_x, arrow_start_y),
                        arrowprops=dict(arrowstyle=arrow_style,
                                       color=color,
                                       linewidth=linewidth,
                                       mutation_scale=head_size,
                                       shrinkA=0,
                                       shrinkB=0),
                        zorder=zorder)
    
    return arrow

# Example usage
def example_usage():
    """Example demonstrating different options."""
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # 定义线段
    x_start, y_start = 0, 0
    x_end, y_end = 1, 1
    
    # 示例1: 不显示标记点（默认）
    ax1.set_title("Without marker (default)", fontsize=12)
    ax1.plot([x_start, x_end], [y_start, y_end], 'gray', linewidth=2, alpha=0.3)
    
    # 绘制多个箭头，不显示标记点
    positions = [0.2, 0.5, 0.8]
    colors = ['red', 'green', 'blue']
    
    for pos, color in zip(positions, colors):
        draw_arrow_on_line(ax1, x_start, y_start, x_end, y_end, 
                          position=pos, color=color, show_marker=False)
    
    ax1.set_xlim(-0.1, 1.1)
    ax1.set_ylim(-0.1, 1.1)
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3)
    
    # 示例2: 显示标记点
    ax2.set_title("With markers", fontsize=12)
    ax2.plot([x_start, x_end], [y_start, y_end], 'gray', linewidth=2, alpha=0.3)
    
    for pos, color in zip(positions, colors):
        draw_arrow_on_line(ax2, x_start, y_start, x_end, y_end, 
                          position=pos, color=color, show_marker=True)
    
    ax2.set_xlim(-0.1, 1.1)
    ax2.set_ylim(-0.1, 1.1)
    ax2.set_aspect('equal')
    ax2.grid(True, alpha=0.3)
    
    # 示例3: 使用简洁版本
    ax3.set_title("Using simple version", fontsize=12)
    ax3.plot([x_start, x_end], [y_start, y_end], 'gray', linewidth=2, alpha=0.3)
    
    # 使用简洁版本
    for pos, color in zip(positions, colors):
        draw_arrow_on_line_simple(ax3, x_start, y_start, x_end, y_end, 
                                 position=pos, color=color)
    
    ax3.set_xlim(-0.1, 1.1)
    ax3.set_ylim(-0.1, 1.1)
    ax3.set_aspect('equal')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

# 修改你的原始示例
def your_example():
    """修改后的你的示例"""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # 定义线段
    line_x = [0, 1]
    line_y = [0, 1]
    
    # 绘制线段
    ax.plot(line_x, line_y, 'gray', linewidth=3, alpha=0.3, label='Base line')
    
    # 绘制箭头（不显示标记点）
    draw_arrow_on_line(ax, 
                       x_start=line_x[0], y_start=line_y[0],
                       x_end=line_x[1], y_end=line_y[1],
                       position=0.5,  # 中间位置
                       color='blue',
                       show_marker=False  # 不显示标记点
                      )
    
    # 设置图形属性
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(-0.5, 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    ax.set_title("Arrow without marker point")
    
    plt.tight_layout()
    plt.show()

# 更多箭头示例
def multiple_arrows_example():
    """展示多个箭头且不显示标记点的例子"""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 定义多个线段
    lines = [
        ((0, 0), (4, 1), 'Line 1'),
        ((1, 0), (3, 3), 'Line 2'),
        ((0, 2), (4, 0), 'Line 3'),
    ]
    
    colors = ['blue', 'green', 'red']
    
    for i, ((x1, y1), (x2, y2), label) in enumerate(lines):
        color = colors[i]
        
        # 绘制线段
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=2, alpha=0.3, label=label)
        
        # 在线段上绘制3个箭头（不显示标记点）
        for position in [0.25, 0.5, 0.75]:
            draw_arrow_on_line(ax, x1, y1, x2, y2,
                              position=position,
                              arrow_style='-|>',
                              color=color,
                              linewidth=1.5,
                              head_size=12,
                              show_marker=False)
    
    # 设置图形属性
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-0.5, 3.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    ax.set_title("Multiple Arrows Without Marker Points")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("示例1: 你的原始示例（修改版）")
    your_example()
    
    print("\n示例2: 不同选项对比")
    example_usage()
    
    print("\n示例3: 多个箭头示例")
    multiple_arrows_example()