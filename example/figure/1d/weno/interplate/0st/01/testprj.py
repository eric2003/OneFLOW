import numpy as np
import matplotlib.pyplot as plt

# ========================== 可调参数（只改这里！） ==========================
fig_width  = 16.0          # x方向宽度（英寸）
fig_height = 8.0           # y方向高度（英寸），随意改实现任意长宽比
# 示例：20x8、18x6、24x9 等都完美自适应不交叉

k = 3
# =========================================================================

# === 缩放因子（以 16×8 为基准）===
base_width  = 16.0
base_height = 8.0
scale_x = fig_width / base_width
scale_y = fig_height / base_height

# === 水平尺寸（图形略微放大 → 边距自然极小）===
visual_cell_width = fig_width / 6.9      # 6.9 → 边距更接近 Word “窄”边距（≈0.4~0.6cm）
center_x = 3.0 * visual_cell_width       # 完美居中（左2.5 + 右3.5 的平均）

# === 垂直尺寸 ===
base_dyref = 1.85
rv, sv = [], []
kk = k - 1
for m in range(0, k + 1):
    s_val = m
    r_val = kk - s_val
    rv.append(r_val)
    sv.append(s_val)
num_rows = len(rv)

dyref = base_dyref * scale_y * (4.0 / max(num_rows, 1))
vertical_unit = dyref * 0.54

def plot_cell_center_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    xs = [center_x + m * visual_cell_width for m in ms]
    plt.scatter(xs, np.full_like(xs, yref), s=140*scale_x**2,
                facecolor='black', edgecolor='black', linewidth=1.2*scale_x)

def plot_mesh_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    if not ms:
        return
    dy = 0.25 * vertical_unit
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))

    for v in v_rels:
        xv = center_x + v * visual_cell_width
        plt.plot([xv, xv], [yref - dy, yref + dy], 'k-', linewidth=1.8*scale_x)

    for m in ms:
        left  = center_x + (m - 0.5) * visual_cell_width
        right = center_x + (m + 0.5) * visual_cell_width
        plt.plot([left, right], [yref, yref], 'b-', linewidth=2.8*scale_x)

def plot_label_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    if not ms:
        return

    y_vertex = yref + 0.25 * vertical_unit
    y_cell   = yref - 0.68 * vertical_unit
    y_rs     = yref - 1.05 * vertical_unit

    # vertex 标签
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:
        frac = int(round(abs(v) * 2))
        sign = '+' if v > 0 else '-'
        label = rf'$x_{{i{sign}\frac{{{frac}}}{{2}}}}$'
        xv = center_x + v * visual_cell_width
        plt.text(xv, y_vertex, label, fontsize=14*scale_y, ha='center', va='bottom')

    # cell 标签
    for m in ms:
        xc = center_x + m * visual_cell_width
        if m == 0:
            plt.text(xc, y_cell, r'$i$', fontsize=16*scale_y, ha='center',
                     color='red', weight='bold')
        elif m > 0:
            plt.text(xc, y_cell, rf'$i+{m}$', fontsize=14*scale_y, ha='center')
        else:
            plt.text(xc, y_cell, rf'$i-{-m}$', fontsize=14*scale_y, ha='center')

    # (r,s)
    shift = (-r + s) / 2.0
    xc = center_x + shift * visual_cell_width
    plt.text(xc, y_rs, rf'$(r={r},\;s={s})$', fontsize=15*scale_y, ha='center',
             color='darkblue', weight='bold')

# ========================== 主程序 ==========================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])
plt.figure(figsize=(fig_width, fig_height))

for i in range(num_rows):
    yref = -i * dyref
    r = rv[i]
    s = sv[i]
    plot_cell_center_rs(yref, r, s)
    plot_mesh_rs(yref, r, s)
    plot_label_rs(yref, r, s)

# === 极小边距（真正像 Word “窄”边距）===
margin_x = 0.12 * visual_cell_width   # 左右 ≈0.3~0.5cm
margin_y = 0.25 * vertical_unit       # 上下极小（保证不裁剪最后一行的 (r,s)）

min_x = center_x - 2.5 * visual_cell_width - margin_x
max_x = center_x + 3.5 * visual_cell_width + margin_x
min_y = -(num_rows-1)*dyref - 1.3*vertical_unit - margin_y
max_y = 0 + 0.4*vertical_unit + margin_y

plt.xlim(min_x, max_x)
plt.ylim(min_y, max_y)
plt.axis('off')

plt.savefig('cfd_stencil_word_narrow.png', bbox_inches='tight', pad_inches=0.02, dpi=400)
plt.show()