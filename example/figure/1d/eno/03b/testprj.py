import numpy as np
import matplotlib.pyplot as plt

# ========================== 可调参数（只需调这几个） ==========================
visual_cell_width = 2.0    # 调大一点让整体占满（2.0 在 12 英寸宽下几乎完美左右贴边）
center_x = 6.0            # 6.0 基本居中（因为左右延伸不对称 2.5 vs 3.5）
dyref = 1.5               # 行间垂直间距
k = 3                    # 保持 3 → 产生 r = 2,1,0,-1 四行
# =========================================================================

def plot_cell_center_rs(yref, r, s):
    ms = list(range(-r, s + 1))                     # 自动支持 r 负（Python range 能处理）
    xs = [center_x + m * visual_cell_width for m in ms]
    plt.scatter(xs, np.full_like(xs, yref), s=120, facecolor='black', edgecolor='black')

def plot_mesh_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    if not ms:
        return
    dy = 0.12 * visual_cell_width

    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))

    # 垂直黑线
    for v in v_rels:
        xv = center_x + v * visual_cell_width
        plt.plot([xv, xv], [yref - dy, yref + dy], 'k-', linewidth=1.5)

    # 水平蓝线
    for m in ms:
        left  = center_x + (m - 0.5) * visual_cell_width
        right = center_x + (m + 0.5) * visual_cell_width
        plt.plot([left, right], [yref, yref], 'b-', linewidth=2.4)

def plot_label_rs(yref, r, s):
    ms = list(range(-r, s + 1))
    if not ms:
        return

    dyb = 0.5 * visual_cell_width
    yb = yref - dyb - 0.1   # (r,s) 稍微低一点
    ybc = yref - 0.5 * dyb
    ytt = yref + 0.52 * dyb   # vertex 标签在最上方

    # === vertex 标签（只画当前行实际存在的） ===
    v_rels = sorted(set(m - 0.5 for m in ms) | set(m + 0.5 for m in ms))
    for v in v_rels:
        offset = v
        frac = int(round(abs(offset) * 2))  # 1,3,5,7
        sign = '+' if offset > 0 else '-'
        label = rf'$x_{{i{sign}\frac{{{frac}}}{{2}}}}$'
        xv = center_x + offset * visual_cell_width
        plt.text(xv, ytt, label, fontsize=12, ha='center', va='bottom')

    # === cell 标签 ===
    for m in ms:
        xc = center_x + m * visual_cell_width
        if m == 0:
            plt.text(xc, ybc, r'$i$', fontsize=14, ha='center', color='red')
        elif m > 0:
            plt.text(xc, ybc, rf'$i+{m}$', fontsize=12, ha='center')
        else:
            plt.text(xc, ybc, rf'$i-{ -m }$', fontsize=12, ha='center')

    # (r,s) 标注
    shift = (-r + s) / 2.0
    xc = center_x + shift * visual_cell_width
    plt.text(xc, yb, rf'$(r={r},\; s={s})$', fontsize=13, ha='center', color='darkblue')

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
plt.figure(figsize=(12, 6.5))  # 稍微高一点容纳 4 行

rv, sv = [], []
getrs(k, rv, sv)
print(f'生成的 (r,s) 顺序: {list(zip(rv, sv))}')  # [ (2,0), (1,1), (0,2), (-1,3) ]

for i in range(len(rv)):
    yref = -i * dyref
    r = rv[i]
    s = sv[i]
    plot_cell_center_rs(yref, r, s)
    plot_mesh_rs(yref, r, s)
    plot_label_rs(yref, r, s)

plt.axis('off')
plt.xlim(-1, 13)      # 基本贴边（2.0 时 -5→+7 相对单位）
plt.ylim(-len(rv)*dyref - 1, 1.5)
plt.tight_layout()
plt.savefig('cfd_stencil_bias.png', bbox_inches='tight', dpi=400)
plt.show()