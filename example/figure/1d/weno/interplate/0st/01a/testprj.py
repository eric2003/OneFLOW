import numpy as np
import matplotlib.pyplot as plt

# ========================== 可调参数（只改这里！） ==========================
fig_width  = 16.0    # 图形宽度（英寸），比如 16、20、24 都行
fig_height = 4.0     # 图形高度（英寸），单行建议 3~5 之间好看
# =========================================================================

# 以 16×8 为基准的缩放（即使改成其他比例也完全自适应）
scale_x = fig_width / 16.0
scale_y = fig_height / 8.0

# 水平布局：5 个点，视觉上等距，左右边距接近 Word “窄”
visual_cell_width = fig_width / 7.2      # 7.2 这个数让左右边距≈0.4~0.6cm（Word窄）
center_x = fig_width / 2                 # 完美水平居中

# 垂直方向只需要一行，高度微调让标签不拥挤
vertical_unit = 1.6 * scale_y            # 控制标签与点的相对距离

# 主绘图函数（只画一行）
def plot_single_row():
    # 5 个网格点的位置：i-2, i-1, i, i+1, i+2
    offsets = [-2, -1, 0, 1, 2]
    xs = [center_x + m * visual_cell_width for m in offsets]

    # 1. 画网格线（蓝色粗横线 + 黑色细竖线）
    y_line = 0.0
    # 横线（覆盖整个 5 个 cell）
    plt.hlines(y_line, xs[0] - 0.5*visual_cell_width,
                     xs[-1] + 0.5*visual_cell_width,
                     color='b', linewidth=2.8*scale_x)

    # 竖线（每个点左右各一条细黑线）
    for x in xs:
        plt.vlines(x, y_line - 0.25*vertical_unit,
                        y_line + 0.25*vertical_unit,
                        color='k', linewidth=1.8*scale_x)

    # 2. 画中心黑点
    plt.scatter(xs, [0]*5, s=140*scale_x**2,
                facecolor='black', edgecolor='black', linewidth=1.2*scale_x)

    # 3. 标签（cell 标签：i-2, i-1, i, i+1, i+2）
    for m, x in zip(offsets, xs):
        if m == 0:
            plt.text(x, 0.68*vertical_unit, r'$i$', fontsize=16*scale_y,
                     ha='center', va='bottom', color='red', weight='bold')
        elif m > 0:
            plt.text(x, 0.68*vertical_unit, rf'$i+{m}$', fontsize=14*scale_y,
                     ha='center', va='bottom')
        else:
            plt.text(x, 0.68*vertical_unit, rf'$i{-m}$', fontsize=14*scale_y,
                     ha='center', va='bottom')

# ========================== 绘图 ==========================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])
fig = plt.figure(figsize=(fig_width, fig_height))

plot_single_row()

# 极窄边距（真正像 Word “窄”）
margin_x = 0.12 * visual_cell_width
margin_y = 0.4 * vertical_unit
plt.xlim(center_x - 2.5*visual_cell_width - margin_x,
         center_x + 2.5*visual_cell_width + margin_x)
plt.ylim(-1.2*vertical_unit - margin_y, 1.6*vertical_unit + margin_y)

plt.axis('off')
plt.savefig('cfd_stencil_5point_single_row.png', bbox_inches='tight',
            pad_inches=0.02, dpi=400)
plt.show()