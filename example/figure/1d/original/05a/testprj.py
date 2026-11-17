import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.patches import PathPatch

def draw_curly_brace(x1, x2, y_attach, text, height=0.6, control_h=0.35, upward=True):
    ax = plt.gca()
    mid = (x1 + x2) / 2
    direction = 1 if upward else -1
    
    tip_y = y_attach + direction * height
    curve_y = tip_y + direction * control_h
    
    # 左半部分
    verts_left = [(x1, y_attach), (x1, tip_y), (mid, curve_y), (mid, tip_y)]
    codes_left = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3]
    ax.add_patch(PathPatch(Path(verts_left, codes_left), facecolor='none', edgecolor='k', lw=1.8))
    
    # 右半部分
    verts_right = [(mid, tip_y), (mid, curve_y), (x2, tip_y), (x2, y_attach)]
    codes_right = [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO]
    ax.add_patch(PathPatch(Path(verts_right, codes_right), facecolor='none', edgecolor='k', lw=1.8))
    
    # 文字统一放在同一高度（每个子图内部 r 标签在同一水平线）
    text_y_pos = y_attach + 2.0
    va = 'bottom' if upward else 'top'
    ax.text(mid, text_y_pos, text, ha='center', va=va, fontsize=12)

def plot_mixed_line(xst, xed, y0):
    ds = xed - xst
    ds1 = ds / 4
    x1 = xst + ds1
    x2 = xed - ds1
    plt.plot([xst, x1], [y0, y0], 'k-', linewidth=1)
    plt.plot([x1, x2], [y0, y0], 'k--', linewidth=1)
    plt.plot([x2, xed], [y0, y0], 'k-', linewidth=1)

def plot_cfd_line(x_points, y0):
    # 边缘红点
    edge_points_red = np.concatenate([x_points[:2], x_points[10:]])
    plt.scatter(edge_points_red, np.full_like(edge_points_red, y0), s=100,
                facecolor='red', edgecolor='black', linewidth=1)
    
    # 特殊黑点 i=-4 和 i=4
    special_black_points = np.array([-4, 4], dtype=np.float64)
    plt.scatter(special_black_points, np.full_like(special_black_points, y0), s=100,
                facecolor='black', edgecolor='black', linewidth=1)
    
    # 中间黑点 i=-2 到 i=3
    middle_points = x_points[3:9]  # [-2, -1, 0, 1, 2, 3]
    plt.scatter(middle_points, np.full_like(middle_points, y0), s=100,
                facecolor='black', edgecolor='black', linewidth=1)
    plt.plot(middle_points, np.full_like(middle_points, y0), 'k-', linewidth=1)
    
    # 左右混合虚实线
    plot_mixed_line(-4, -2, y0)
    plot_mixed_line(2, 4, y0)

def plot_rect(x, y, width, height, rectangle_color):
    rect = patches.FancyBboxPatch((x, y), width, height,
                                  boxstyle="round,pad=0.1,rounding_size=0.1",
                                  edgecolor='none', facecolor=rectangle_color, zorder=0)
    plt.gca().add_patch(rect)

def plot_label(y0, xv, lr, txt_name):
    plt.text(-6, y0-0.5, '$-2$', fontsize=12, ha='center')
    plt.text(-5, y0-0.5, '$-1$', fontsize=12, ha='center')
    plt.text(-4, y0-0.5, '$i=1$', fontsize=12, ha='center')   # 你原代码就是这样，保留
    plt.text(-2, y0-0.5, '$i-2$', fontsize=12, ha='center')
    plt.text(-1, y0-0.5, '$i-1$', fontsize=12, ha='center')
    plt.text(0,  y0-0.5, '$i$', fontsize=12, ha='center')
    plt.text(1,  y0-0.5, '$i+1$', fontsize=12, ha='center')
    plt.text(2,  y0-0.5, '$i+2$', fontsize=12, ha='center')
    plt.text(4,  y0-0.5, '$i=N+1$', fontsize=12, ha='center')
    plt.text(5,  y0-0.5, '$N+2$', fontsize=12, ha='center')
    plt.text(6,  y0-0.5, '$N+3$', fontsize=12, ha='center')

    lrname = r'$u_{i+\frac{1}{2}}^' + lr + r'$'
    plt.text(xv, y0+0.5, lrname, fontsize=12, ha='center')
    plt.text(-4, y0+0.3, '$x=0$', fontsize=12, ha='center')
    plt.text(+4, y0+0.3, '$x=L$', fontsize=12, ha='center')
    plt.text(0, y0-1.5, txt_name, fontsize=12, ha='center')

# ====================== 主程序 ======================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

plt.figure(figsize=(12, 10))   # 高度加大，留出大括号空间

x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6], dtype=np.float64)
xv = 0.5
rectangle_color = (150/255, 150/255, 200/255)

# ------------------- 上图 (a) -------------------
y0 = 3.0
plot_cfd_line(x_points, y0)
plt.plot([xv, xv], [y0-1.8, y0+2.2], 'k--', linewidth=1.5)   # 稍加长
plot_rect(-2.2, y0-0.1, 4.4, 0.2, rectangle_color)
plot_label(y0, xv-0.3, 'L', '(a) Left-side reconstruction')

# 大括号（向上，r 标签高度完全一致）
draw_curly_brace(-2, 0, y0, r'$r=2$', height=1.3, control_h=0.35)
draw_curly_brace(-1, 1, y0, r'$r=1$', height=0.9, control_h=0.35)
draw_curly_brace(0, 2, y0, r'$r=0$', height=0.55, control_h=0.35)

# ------------------- 下图 (b) -------------------
y1 = -4.0
plot_cfd_line(x_points, y1)
plt.plot([xv, xv], [y1-1.8, y1+2.2], 'k--', linewidth=1.5)
plot_rect(-1.2, y1-0.1, 4.4, 0.2, rectangle_color)
plot_label(y1, xv+0.3, 'R', '(b) Right-side reconstruction')

# 大括号（同样参数）
draw_curly_brace(-1, 1, y1, r'$r=2$', height=1.3, control_h=0.35)
draw_curly_brace(0,  2, y1, r'$r=1$', height=0.9, control_h=0.35)
draw_curly_brace(1,  3, y1, r'$r=0$', height=0.55, control_h=0.35)

plt.axis('equal')
plt.axis('off')
plt.savefig('cfd_final_with_braces.png', bbox_inches='tight', dpi=300)
plt.show()