import numpy as np
import matplotlib.pyplot as plt

# ========================== 可调参数（只需调这几个） ==========================
visual_cell_width = 2.0    # 2.0 时左右几乎完美贴边（12英寸宽）
center_x = 5.0             # 精确居中（左2.5 + 右3.5 = 6个单元 → center=2.5×2.0=5.0）
dyref = 2.35               # 行间距，彻底消除垂直交叉（可继续调大更疏朗）
k = 3                      # 保持3 → 产生你想要的4种偏置
# =========================================================================

def plot_cell_center_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    xs = [center_x + m * visual_cell_width for m in ms]
    plt.scatter(xs, np.full_like(xs, yref), s=140, facecolor='black', edgecolor='black', linewidth=1)

def plot_mesh_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    if not ms:
        return
    dy = 0.12 * visual_cell_width

    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))

    # 垂直黑线
    for v in v_rels:
        xv = center_x + v * visual_cell_width
        plt.plot([xv, xv], [yref - dy, yref + dy], 'k-', linewidth=1.6)

    # 水平蓝线（cell）
    for m in ms:
        left  = center_x + (m - 0.5) * visual_cell_width
        right = center_x + (m + 0.5) * visual_cell_width
        plt.plot([left, right], [yref, yref], 'b-', linewidth=2.6)

def plot_label_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    if not ms:
        return

    dyb = 0.5 * visual_cell_width                  # 1.0
    y_vertex = yref + 0.68 * dyb                    # vertex 标签（最上方）
    y_cell   = yref - 0.68 * dyb                    # cell 标签（中间偏下）
    y_rs     = yref - 1.05 * dyb                    # (r,s) 标签（最下方，保证不交叉）

    # === vertex 标签（只画当前行实际存在的）===
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:
        frac = int(round(abs(v) * 2))
        sign = '+' if v > 0 else '-'
        label = rf'$x_{{i{sign}\frac{{{frac}}}{{2}}}}$'
        xv = center_x + v * visual_cell_width
        plt.text(xv, y_vertex, label, fontsize=13, ha='center', va='bottom')

    # === cell 标签 ===
    for m in ms:
        xc = center_x + m * visual_cell_width
        if m == 0:
            plt.text(xc, y_cell, r'$i$', fontsize=15, ha='center', color='red', weight='bold')
        elif m > 0:
            plt.text(xc, y_cell, rf'$i+{m}$', fontsize=13, ha='center')
        else:
            plt.text(xc, y_cell, rf'$i-{-m}$', fontsize=13, ha='center')

    # === (r,s) 标注（放在最下方，永远不会和下一行 vertex 交叉）===
    shift = (-r + s) / 2.0
    xc = center_x + shift * visual_cell_width
    plt.text(xc, y_rs, rf'$(r={r},\;s={s})$', fontsize=14, ha='center', color='darkblue', weight='bold')

def getrs(k, rv, sv):
    kk = k - 1
    for m in range(0, k + 1):
        s_val = m
        r_val = kk - s_val
        rv.append(r_val)
        sv.append(s_val)

# ========================== 主程序 ==========================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])
plt.figure(figsize=(12, 7.8))   # 稍微高一点，4行更舒适

rv, sv = [], []
getrs(k, rv, sv)
print(f'生成的 (r,s) 顺序: {list(zip(rv, sv))}')  # 你会看到 [(2, 0), (1, 1), (0, 2), (-1, 3)]

for i in range(len(rv)):
    yref = -i * dyref
    r = rv[i]
    s = sv[i]
    plot_cell_center_rs(yref, r, s)
    plot_mesh_rs(yref, r, s)
    plot_label_rs(yref, r, s)

plt.axis('off')
plt.tight_layout()
plt.savefig('cfd_stencil_final.png', bbox_inches='tight', dpi=400)
plt.show()