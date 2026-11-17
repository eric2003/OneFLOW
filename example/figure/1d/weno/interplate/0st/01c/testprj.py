import numpy as np
import matplotlib.pyplot as plt

# ========================== 可调参数（只改这里！） ==========================
fig_width  = 16.0    # 随意改宽度
fig_height = 4.0     # 单行建议 3.8~4.5，最好别低于 3.5
# =========================================================================

# 完全复刻你原始代码的缩放与布局参数（一点都没改！）
base_width = 16.0
base_height = 8.0
scale_x = fig_width / base_width
scale_y = fig_height / base_height

visual_cell_width = fig_width / 6.9     # 原始值，保证 Word 窄边距
center_x = 3.0 * visual_cell_width      # 原始完美居中（左2.5 + 右3.5 的中点）

ffsize = 30

base_dyref = 1.85
dyref = base_dyref * scale_y * (4.0 / 4)   # 原始公式，这里固定 4 使单行间距与原来一行完全一致
vertical_unit = dyref * 0.54               # 完全不变
yref = 0.0                                 # 只画一行

# ========================== 绘图函数（100% 复刻原始逻辑） ==========================
def plot_cell_center_rs(yref, r=2, s=2):
    ms = list(range(-r, s + 1))            # [-2,-1,0,1,2]
    xs = [center_x + m * visual_cell_width for m in ms]
    plt.scatter(xs, np.full_like(xs, yref), s=140*scale_x**2,
                facecolor='black', edgecolor='black', linewidth=1.2*scale_x)

def plot_mesh_rs(yref, r=2, s=2):
    ms = list(range(-r, s + 1))
    dy = 0.1 * vertical_unit
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:                               # 画黑色细竖线（face位置）
        xv = center_x + v * visual_cell_width
        plt.plot([xv, xv], [yref - dy, yref + dy], 'k-', linewidth=1.8*scale_x)
    
    for m in ms:                                   # 画蓝色粗横线（每个cell）
        left  = center_x + (m - 0.5) * visual_cell_width
        right = center_x + (m + 0.5) * visual_cell_width
        plt.plot([left, right], [yref, yref], 'b-', linewidth=2.8*scale_x)

def plot_label_rs(yref, r=2, s=2):
    ms = list(range(-r, s + 1))
    y_vertex = yref + 0.2 * vertical_unit   # x_{i±1/2} 在竖线上方，和原来完全一样
    y_cell   = yref - 0.25 * vertical_unit   # i, i±1 在点下方，和原来完全一样

    # 1. face 标签 x_{i-3/2}, x_{i-1/2}, x_{i+1/2}, x_{i+3/2}, x_{i+5/2}
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:
        frac = int(round(abs(v) * 2))        # 1,3,5 → 1/2,3/2,5/2
        sign = '+' if v > 0 else '-'
        label = rf'$x_{{i{sign}\frac{{{frac}}}{{2}}}}$'
        xv = center_x + v * visual_cell_width
        plt.text(xv, y_vertex, label, fontsize=ffsize*scale_y, ha='center', va='bottom')

    # 2. cell 标签 i-2, i-1, i, i+1, i+2
    for m in ms:
        xc = center_x + m * visual_cell_width
        if m == 0:
            plt.text(xc, y_cell, r'$i$', fontsize=ffsize*scale_y, ha='center',
                     color='black', weight='bold')
        elif m > 0:
            plt.text(xc, y_cell, rf'$i+{m}$', fontsize=ffsize*scale_y, ha='center')
        else:
            plt.text(xc, y_cell, rf'$i-{-m}$', fontsize=ffsize*scale_y, ha='center')

# ========================== 主程序（和原来一模一样） ==========================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])
plt.figure(figsize=(fig_width, fig_height))

plot_cell_center_rs(yref)
plot_mesh_rs(yref)
plot_label_rs(yref)

# 边距也完全复刻原始代码（极窄，插入Word完美）
margin_x = 0.12 * visual_cell_width
margin_y = 0.25 * vertical_unit
min_x = center_x - 2.5 * visual_cell_width - margin_x
max_x = center_x + 3.5 * visual_cell_width + margin_x   # 保持原始左右不对称（右边多留一点，和原来一致）
min_y = -1.3*vertical_unit - margin_y
max_y =  0.4*vertical_unit + margin_y

plt.xlim(min_x, max_x)
plt.ylim(min_y, max_y)
plt.axis('off')
plt.savefig('cfd_5point_perfect_same_as_original.png', bbox_inches='tight',
            pad_inches=0.02, dpi=400)
plt.show()