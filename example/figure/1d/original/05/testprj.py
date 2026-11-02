import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches
from matplotlib.path import Path          # 新增导入
from matplotlib.patches import PathPatch   # 新增导入

def draw_curly_brace(x1, x2, y_attach, text, height=0.6, control_h=0.4, text_offset=0.5):
    ax = plt.gca()
    mid = (x1 + x2) / 2

    # 左半边
    verts_left = [(x1, y_attach),
                  (x1, y_attach - height),
                  (mid, y_attach - height - control_h),
                  (mid, y_attach - height)]
    codes_left = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3]
    ax.add_patch(PathPatch(Path(verts_left, codes_left), facecolor='none', edgecolor='k', lw=1.5))

    # 右半边
    verts_right = [(mid, y_attach - height),
                   (mid, y_attach - height - control_h),
                   (x2, y_attach - height),
                   (x2, y_attach)]
    codes_right = [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO]
    ax.add_patch(PathPatch(Path(verts_right, codes_right), facecolor='none', edgecolor='k', lw=1.5))

    # 文字
    ax.text(mid, y_attach - height - control_h - text_offset, text, ha='center', va='top', fontsize=12)

def plot_mixed_line(xst, xed, y0):
    ds = xed - xst
    ds1 = ds / 4
    x1 = xst + ds1
    x2 = xed - ds1
    x = np.array([xst, x1, x2, xed])
    plt.plot([x[0], x[1]], [y0, y0], 'k-', linewidth=1)
    plt.plot([x[1], x[2]], [y0, y0], 'k--', linewidth=1)
    plt.plot([x[2], x[3]], [y0, y0], 'k-', linewidth=1)

def plot_cfd_line(x_points, y0):
    edge_points_red = np.concatenate([x_points[:2], x_points[10:]])
    plt.scatter(edge_points_red, np.full_like(edge_points_red, y0), s=100, facecolor='red', edgecolor='black', linewidth=1)

    special_black_points = np.array([-4, 4], dtype=np.float64)
    plt.scatter(special_black_points, np.full_like(special_black_points, y0), s=100, facecolor='black', edgecolor='black', linewidth=1)

    middle_points = x_points[3:9]
    plt.scatter(middle_points, np.full_like(middle_points, y0), s=100, facecolor='black', edgecolor='black', linewidth=1)
    plt.plot(middle_points, np.full_like(middle_points, y0), 'k-', linewidth=1)

    plot_mixed_line(-4, -2, y0)
    plot_mixed_line(2, 4, y0)

def plot_rect(x, y, width, height, rectangle_color):
    rect = patches.FancyBboxPatch((x, y), width, height,
                                  boxstyle="round,pad=0.1,rounding_size=0.1",
                                  edgecolor='none', facecolor=rectangle_color, zorder=0)
    ax = plt.gca()
    ax.add_patch(rect)

def plot_label(y0, xv, lr, txt_name):
    plt.text(-6, y0-0.5, '$-2$', fontsize=12, ha='center')
    plt.text(-5, y0-0.5, '$-1$', fontsize=12, ha='center')
    plt.text(-4, y0-0.5, '$i=1$', fontsize=12, ha='center')   # 你原代码这里是i=1，可能打错，我保留
    plt.text(-2, y0-0.5, '$i-2$', fontsize=12, ha='center')
    plt.text(-1, y0-0.5, '$i-1$', fontsize=12, ha='center')
    plt.text(0, y0-0.5, '$i$', fontsize=12, ha='center')
    plt.text(1, y0-0.5, '$i+1$', fontsize=12, ha='center')
    plt.text(2, y0-0.5, '$i+2$', fontsize=12, ha='center')
    plt.text(4, y0-0.5, '$i=N+1$', fontsize=12, ha='center')
    plt.text(5, y0-0.5, '$N+2$', fontsize=12, ha='center')
    plt.text(6, y0-0.5, '$N+3$', fontsize=12, ha='center')

    lrname = r'$u_{i+\frac{1}{2}}^' + lr + r'$'
    plt.text(xv, y0+0.5, lrname, fontsize=12, ha='center')
    plt.text(-4, y0+0.3, '$x=0$', fontsize=12, ha='center')
    plt.text(+4, y0+0.3, '$x=L$', fontsize=12, ha='center')
    plt.text(0, y0-1.5, txt_name, fontsize=12, ha='center')

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

plt.figure(figsize=(12, 5))

x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6], dtype=np.float64)

y0 = 0.5
plot_cfd_line(x_points, y0)
xv = 0.5 * (x_points[5] + x_points[6])   # 0.5
plt.plot([xv, xv], [y0-0.5, y0+0.5], 'k--')
rectangle_color = (150/255, 150/255, 200/255)
plot_rect(-2.2, y0-0.1, 4.4, 0.2, rectangle_color)
plot_label(y0, xv-0.3, 'L', '(a) Left-side reconstruction')

# 上图添加三层大括号（从外到内）
draw_curly_brace(-2, 0, y0-0.5, '$r=2$', height=0.9, control_h=0.7, text_offset=0.6)
draw_curly_brace(-1, 1, y0-0.5, '$r=1$', height=0.7, control_h=0.5, text_offset=0.5)
draw_curly_brace(0, 2, y0-0.5, '$r=0$', height=0.5, control_h=0.3, text_offset=0.4)

y1 = -2.0
plot_cfd_line(x_points, y1)
plt.plot([xv, xv], [y1-0.5, y1+0.5], 'k--')
plot_rect(-1.2, y1-0.1, 4.4, 0.2, rectangle_color)
plot_label(y1, xv+0.3, 'R', '(b) Right-side reconstruction')

# 下图添加三层大括号（从外到内）
draw_curly_brace(-1, 1, y1-0.5, '$r=2$', height=0.9, control_h=0.7, text_offset=0.6)
draw_curly_brace(0, 2, y1-0.5, '$r=1$', height=0.7, control_h=0.5, text_offset=0.5)
draw_curly_brace(1, 3, y1-0.5, '$r=0$', height=0.5, control_h=0.3, text_offset=0.4)

plt.axis('equal')
plt.axis('off')
plt.savefig('cfd_with_braces.png', bbox_inches='tight', dpi=300)
plt.show()