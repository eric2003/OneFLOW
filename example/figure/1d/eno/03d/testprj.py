import numpy as np
import matplotlib.pyplot as plt

# ========================== 可调参数 ==========================
visual_cell_width = 2.0
center_x = 5.0             # 精确居中（左2.5 + 右3.5 = 6个单元宽）
dyref = 2.35               # 行间距（保证上下完全不交叉）
k = 3
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

    for v in v_rels:
        xv = center_x + v * visual_cell_width
        plt.plot([xv, xv], [yref - dy, yref + dy], 'k-', linewidth=1.6)

    for m in ms:
        left  = center_x + (m - 0.5) * visual_cell_width
        right = center_x + (m + 0.5) * visual_cell_width
        plt.plot([left, right], [yref, yref], 'b-', linewidth=2.6)

def plot_label_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    if not ms:
        return

    dyb = 0.5 * visual_cell_width                  # = 1.0
    # 关键调整：vertex 标签大幅贴近网格线
    y_vertex = yref + 0.28 * dyb     # ← 原来0.68，现在0.28，明显更近
    y_cell   = yref - 0.68 * dyb
    y_rs     = yref - 1.05 * dyb     # (r,s) 仍保持在最下方，安全距离

    # === vertex 标签（只画当前行实际存在的）===
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:
        frac = int(round(abs(v) * 2))
        sign = '+' if v > 0 else '-'
        label = rf'$x_{{i{sign}\frac{{{frac}}}{{2}}}}$'
        xv = center_x + v * visual_cell_width
        plt.text(xv, y_vertex, label, fontsize=13.5, ha='center', va='bottom', color='black')

    # === cell 标签 ===
    for m in ms:
        xc = center_x + m * visual_cell_width
        if m == 0:
            plt.text(xc, y_cell, r'$i$', fontsize=15, ha='center', color='red', weight='bold')
        elif m > 0:
            plt.text(xc, y_cell, rf'$i+{m}$', fontsize=13, ha='center')
        else:
            plt.text(xc, y_cell, rf'$i-{-m}$', fontsize=13, ha='center')

    # === (r,s) ===
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
plt.figure(figsize=(12, 7.8))

rv, sv = [], []
getrs(k, rv, sv)

for i in range(len(rv)):
    yref = -i * dyref
    r = rv[i]
    s = sv[i]
    plot_cell_center_rs(yref, r, s)
    plot_mesh_rs(yref, r, s)
    plot_label_rs(yref, r, s)

plt.axis('off')
plt.tight_layout()
plt.savefig('cfd_stencil_final_tighter.png', bbox_inches='tight', dpi=400)
plt.show()