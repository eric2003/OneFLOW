import numpy as np
import matplotlib.pyplot as plt

# ========================== 可调参数（只改这里！） ==========================
fig_width  = 16.0    # 随意改，18、20、24 都完美自适应
fig_height = 4.0     # 单行建议 3.5~5，太高会空
# =========================================================================

# 缩放（保持和原来完全一致的自适应能力）
scale_x = fig_width / 16.0
scale_y = fig_height / 8.0

visual_cell_width = fig_width / 6.9          # 左右边距≈0.4~0.6cm（Word 窄）
center_x = 3.0 * visual_cell_width           # 原始代码的完美居中方式（左2.5 + 右3.5 的中点）

# 垂直方向单位（和原代码一致）
base_dyref = 1.85
dyref = base_dyref * scale_y * (4.0 / 4)      # 原来是按行数缩放，这里固定为单行
vertical_unit = dyref * 0.54
yref = 0.0                                    # 只画一行，固定在 y=0

# ========================== 绘图函数（严格复刻你原始逻辑） ==========================
def plot_cell_center_rs(yref, r=2, s=2):       # 固定 r=2, s=2 → 对应 i-2 ~ i+2
    ms = list(range(-r, s + 1))               # [-2, -1, 0, 1, 2]
    xs = [center_x + m * visual_cell_width for m in ms]
    plt.scatter(xs, np.full_like(xs, yref), s=140*scale_x**2,
                facecolor='black', edgecolor='black', linewidth=1.2*scale_x)

def plot_mesh_rs(yref, r=2, s=2):
    ms = list(range(-r, s + 1))
    if not ms:
        return
    dy = 0.25 * vertical_unit
    # 竖线：位于 face 位置（m-0.5 和 m+0.5）
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:
        xv = center_x + v * visual_cell_width
        plt.plot([xv, xv], [yref - dy, yref + dy], 'k-', linewidth=1.8*scale_x)
    
    # 蓝色粗横线：每个 cell 内部
    for m in ms:
        left  = center_x + (m - 0.5) * visual_cell_width
        right = center_x + (m + 0.5) * visual_cell_width
        plt.plot([left, right], [yref, yref], 'b-', linewidth=2.8*scale_x)

def plot_label_rs(yref, r=2, s=2):
    ms = list(range(-r, s + 1))
    y_vertex = yref + 0.25 * vertical_unit   # x_{i±1/2} 标签位置
    y_cell   = yref - 0.68 * vertical_unit   # i, i±1 标签位置

    # 1. face 标签：x_{i-1/2}, x_{i+1/2}, ...
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:
        frac = int(round(abs(v) * 2))        # |v| = 0.5,1.5,2.5 → 1,3,5 → 1/2, 3/2, 5/2
        sign = '+' if v > 0 else ('-' if v < 0 else '')
        if frac % 2 == 1:                     # 奇数 → x_{i±1/2}, x_{i±3/2}, ...
            label = rf'$x_{{i{sign}\frac{{{frac}}}{{2}}}}$'
        else:
            label = rf'$x_{{i{sign}{frac//2}}}$'   # 暂时不会出现整数，但保留
        xv = center_x + v * visual_cell_width
        plt.text(xv, y_vertex, label, fontsize=14*scale_y, ha='center', va='bottom')

    # 2. cell 标签：i-2, i-1, i, i+1, i+2（写在 cell 中心）
    for m in ms:
        xc = center_x + m * visual_cell_width
        if m == 0:
            plt.text(xc, y_cell, r'$i$', fontsize=16*scale_y, ha='center',
                     color='red', weight='bold')
        elif m > 0:
            plt.text(xc, y_cell, rf'$i+{m}$', fontsize=14*scale_y, ha='center')
        else:
            plt.text(xc, y_cell, rf'$i-{-m}$', fontsize=14*scale_y, ha='center')

# ========================== 主程序 ==========================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])
plt.figure(figsize=(fig_width, fig_height))

plot_cell_center_rs(yref)
plot_mesh_rs(yref)
plot_label_rs(yref)

# 极窄边距（和原来完全一致）
margin_x = 0.12 * visual_cell_width
margin_y = 0.25 * vertical_unit
min_x = center_x - 2.5 * visual_cell_width - margin_x
max_x = center_x + 2.5 * visual_cell_width + margin_x   # 原来右边是3.5，这里对称改为2.5（5点更紧凑美观）
min_y = -1.3*vertical_unit - margin_y
max_y =  0.4*vertical_unit + margin_y

plt.xlim(min_x, max_x)
plt.ylim(min_y, max_y)
plt.axis('off')
plt.savefig('cfd_5point_stencil_standard.png', bbox_inches='tight', pad_inches=0.02, dpi=400)
plt.show()