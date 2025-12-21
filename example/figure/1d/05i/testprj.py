import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

# 折线括号：横线统一高度 + 颜色区分
def draw_simple_bracket(x1, x2, y_base, text, vert_bottom_y, color='k', linestyle='-', linewidth=1.6):
    top_y = y_base + 0.8                     # 横线完全统一高度
    mid = (x1 + x2) / 2
    
    plt.plot([x1, x1], [vert_bottom_y, top_y], color=color, ls=linestyle, lw=linewidth)
    plt.plot([x1, x2], [top_y, top_y], color=color, ls=linestyle, lw=linewidth)
    plt.plot([x2, x2], [top_y, vert_bottom_y], color=color, ls=linestyle, lw=linewidth)
    
    plt.text(mid, top_y + 0.15, text, ha='center', va='bottom', fontsize=12)

# 其余函数完全不变
def plot_mixed_line(xst, xed, y0):
    ds = xed - xst
    ds1 = ds / 4
    x1 = xst + ds1
    x2 = xed - ds1
    plt.plot([xst, x1], [y0, y0], 'k-', linewidth=1)
    plt.plot([x1, x2], [y0, y0], 'k--', linewidth=1)
    plt.plot([x2, xed], [y0, y0], 'k-', linewidth=1)

def plot_cfd_line(x_points, y0):
    edge_points_red = np.concatenate([x_points[:2], x_points[10:]])
    plt.scatter(edge_points_red, np.full_like(edge_points_red, y0), s=100,
                facecolor='red', edgecolor='black', linewidth=1)
    
    special_black_points = np.array([-4, 4], dtype=np.float64)
    plt.scatter(special_black_points, np.full_like(special_black_points, y0), s=100,
                facecolor='black', edgecolor='black', linewidth=1)
    
    middle_points = x_points[3:9]
    plt.scatter(middle_points, np.full_like(middle_points, y0), s=100,
                facecolor='black', edgecolor='black', linewidth=1)
    plt.plot(middle_points, np.full_like(middle_points, y0), 'k-', linewidth=1)
    
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
    plt.text(-4, y0-0.5, '$i=1$', fontsize=12, ha='center')
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
    plt.text(4, y0+0.3, '$x=L$', fontsize=12, ha='center')
    plt.text(0, y0-1.5, txt_name, fontsize=12, ha='center')

# ====================== 主程序 ======================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

plt.figure(figsize=(12, 7))  # 紧凑布局

x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6], dtype=np.float64)
xv = 0.5
rectangle_color = (150/255, 150/255, 200/255)

# 上图 (a)
y0 = 2.0
plot_cfd_line(x_points, y0)

# 关键修改：垂直虚线向上穿过括号横线 0.25 单位（到 y0+1.05），但不碰到 r 标签
plt.plot([xv, xv], [y0 - 0.8, y0 + 1.05], 'k--', linewidth=1.5)

plot_rect(-2.2, y0-0.15, 4.4, 0.3, rectangle_color)
plot_label(y0, xv-0.3, 'L', '(a) Left-side reconstruction')

# 灰度配色（最推荐）
draw_simple_bracket(-2, 0, y0, r'$r=2$', vert_bottom_y=y0 + 0.05, color='black',      linestyle=':',  linewidth=1.6)
draw_simple_bracket(-1, 1, y0, r'$r=1$', vert_bottom_y=y0 + 0.27, color=(0.3,0.3,0.3), linestyle='--', linewidth=1.6)
draw_simple_bracket( 0, 2, y0, r'$r=0$', vert_bottom_y=y0 + 0.49, color=(0.6,0.6,0.6), linestyle='-',  linewidth=1.6)

# 下图 (b)
y1 = -1.8
plot_cfd_line(x_points, y1)

# 同上处理
plt.plot([xv, xv], [y1 - 0.8, y1 + 1.05], 'k--', linewidth=1.5)

plot_rect(-1.2, y1-0.15, 4.4, 0.3, rectangle_color)
plot_label(y1, xv+0.3, 'R', '(b) Right-side reconstruction')

draw_simple_bracket(-1, 1, y1, r'$r=2$', vert_bottom_y=y1 + 0.05, color='black',      linestyle=':',  linewidth=1.6)
draw_simple_bracket( 0, 2, y1, r'$r=1$', vert_bottom_y=y1 + 0.27, color=(0.3,0.3,0.3), linestyle='--', linewidth=1.6)
draw_simple_bracket( 1, 3, y1, r'$r=0$', vert_bottom_y=y1 + 0.49, color=(0.6,0.6,0.6), linestyle='-',  linewidth=1.6)

plt.axis('equal')
plt.axis('off')
plt.savefig('cfd_final_with_piercing_line.png', bbox_inches='tight', dpi=300)
plt.show()