import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection  # 用于批量绘制（优化性能）

def draw_border_with_slants(
    border_x: float = 5,          # 主竖线的x坐标
    y_start: float = 1,           # 斜线起始y坐标
    y_end: float = 7,             # 斜线结束y坐标
    num_lines: int = 20,          # 斜线数量
    line_length: float = 0.8,     # 斜线长度
    angle: float = 30,            # 斜线与水平方向夹角（度，负数反向）
    border_color: str = 'black',  # 主竖线颜色
    border_width: float = 2,      # 主竖线宽度
    slant_color: str = 'red',     # 斜线颜色
    slant_width: float = 1,       # 斜线宽度
    slant_alpha: float = 0.7,     # 斜线透明度
    x_lim: tuple = (0, 10),       # x轴范围
    y_lim: tuple = (0, 8),        # y轴范围
    title: str = '竖线+短斜线边界',  # 图表标题
    save_path: str = None         # 保存路径（None则不保存）
) -> None:
    """
    绘制带短斜线的竖线边界图（Matplotlib版）
    
    参数说明：
    - border_x: 主竖线的x坐标
    - y_start/y_end: 斜线在y轴的分布范围
    - num_lines: 斜线数量（越多越密集）
    - line_length: 斜线的长度
    - angle: 斜线角度（正数向右倾，负数向左倾）
    - border_color/border_width: 主竖线样式
    - slant_color/slant_width/slant_alpha: 斜线样式
    - x_lim/y_lim: 坐标轴范围
    - title: 图表标题
    - save_path: 保存路径（如'border.png'），None则仅显示
    """
    # 1. 初始化画布
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_xlabel('X轴')
    ax.set_ylabel('Y轴')
    ax.set_title(title)
    ax.grid(alpha=0.3)  # 网格线（可选）
    
    # 2. 绘制主竖线
    ax.axvline(x=border_x, color=border_color, linewidth=border_width, label='边界线')
    
    # 3. 计算斜线坐标（批量生成，优化性能）
    y_positions = np.linspace(y_start, y_end, num_lines)
    angle_rad = np.deg2rad(angle)
    dx = line_length * np.cos(angle_rad)
    dy = line_length * np.sin(angle_rad)
    
    # 生成所有斜线的坐标对（格式：[[(x1,y1), (x2,y2)], ...]）
    lines = []
    for y in y_positions:
        x1, y1 = border_x, y
        x2, y2 = x1 + dx, y1 + dy
        lines.append([(x1, y1), (x2, y2)])
    
    # 4. 批量绘制斜线（比循环plot更高效，尤其斜线数量多时）
    lc = LineCollection(
        lines,
        colors=slant_color,
        linewidths=slant_width,
        alpha=slant_alpha
    )
    ax.add_collection(lc)
    
    # 5. 等比例显示（保证角度准确）+ 图例
    plt.axis('equal')
    ax.legend()
    
    # 6. 保存/显示
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')  # 高分辨率保存
    plt.show()

# ------------------------------
# 调用示例（直接运行即可）
# ------------------------------
if __name__ == '__main__':
    # 示例1：默认参数（向右倾斜的红色斜线）
    draw_border_with_slants()
    
    # 示例2：自定义参数（向左倾斜的蓝色斜线，更密集）
    draw_border_with_slants(
        border_x=4,
        num_lines=30,          # 更多斜线
        angle=-45,             # 向左倾斜
        slant_color='blue',    # 蓝色斜线
        title='向左倾斜的密集斜线边界',
        save_path='custom_border.png'  # 保存到本地
    )
    
    # 示例3：更粗的斜线+黑色主竖线
    draw_border_with_slants(
        slant_width=2,
        border_width=3,
        line_length=1.2,
        angle=60,
        slant_color='green'
    )