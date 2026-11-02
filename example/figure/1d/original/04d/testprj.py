import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_mixed_line( xst, xed, y0 ):
    ds = xed - xst
    ds1 = ds / 4
    x1 = xst + ds1
    x2 = xed - ds1
    x = np.array([xst, x1, x2, xed])  # 按照1/4, 1/2, 1/4的比例分割
    plt.plot([x[0], x[1]], [y0, y0], 'k-', linewidth=1)  # 第一段实线
    plt.plot([x[1], x[2]], [y0, y0], 'k--', linewidth=1) # 第二段虚线
    plt.plot([x[2], x[3]], [y0, y0], 'k-', linewidth=1)  # 第三段实线
    return
   
def plot_cfd_line( x_points, y0, lr ):
    # 绘制除中间 5 个点和特定边缘点外的其他点 (内部红色，边缘黑色)
    edge_points_red = np.concatenate([x_points[:2], x_points[10:]])
    plt.scatter(edge_points_red, np.full_like(edge_points_red, y0), s=100, facecolor='red', edgecolor='black', linewidth=1)
    
    # 绘制左侧第三点 (i=-4) 和右侧第三点 (i=4) 为纯黑色点
    special_black_points = np.array([-4, 4],dtype=np.float64)
    plt.scatter(special_black_points, np.full_like(special_black_points, y0), s=100, facecolor='black', edgecolor='black', linewidth=1)
    
    # 绘制中间 6 个点 (i=-2, -1, 0, 1, 2, 3)
    iend = 8
    if lr == 'R':
        iend = 9
    middle_points = x_points[3:iend]
    plt.scatter(middle_points, np.full_like(middle_points, y0), s=100, facecolor='black', edgecolor='black', linewidth=1)
    
    # 绘制中间 6 个点的黑实线连接
    plt.plot(middle_points, np.full_like(middle_points, y0), 'k-', linewidth=1)
    
    # 添加左起第三点和第四点之间的分段连线（-4到-2）
    plot_mixed_line(-4,-2, y0)
    
    # 添加右起第三点和第四点之间的分段连线（2到4）
    plot_mixed_line(2,4, y0)
    
    return
    
def plot_rect(x, y, width, height, rectangle_color):
    rect = patches.FancyBboxPatch((x, y), width, height, 
                                        boxstyle="round,pad=0.1,rounding_size=0.1",
                                        edgecolor='none',
                                        facecolor=rectangle_color,
                                        zorder=0)
    
    # 将矩形添加到坐标轴上
    ax = plt.gca()
    ax.add_patch(rect)
    return
    
def plot_label(y0, xv, lr, txt_name):
    # 添加标签和其他设置（保持不变）
    plt.text(-6, y0-0.5, '$-2$', fontsize=12, ha='center')
    plt.text(-5, y0-0.5, '$-1$', fontsize=12, ha='center')
    plt.text(-4, y0-0.5, '$i=1$', fontsize=12, ha='center')
    if lr == 'L':
        plt.text(-2, y0-0.5, '$i-2$', fontsize=12, ha='center')
    plt.text(-1, y0-0.5, '$i-1$', fontsize=12, ha='center')
    plt.text(0, y0-0.5, '$i$', fontsize=12, ha='center')
    plt.text(1, y0-0.5, '$i+1$', fontsize=12, ha='center')
    plt.text(2, y0-0.5, '$i+2$', fontsize=12, ha='center')
    if lr == 'R':
        plt.text(3, y0-0.5, '$i+3$', fontsize=12, ha='center')
    
    plt.text(4, y0-0.5, '$i=N+1$', fontsize=12, ha='center')
    plt.text(5, y0-0.5, '$N+2$', fontsize=12, ha='center')
    plt.text(6, y0-0.5, '$N+3$', fontsize=12, ha='center')
    
    lrname = r'$u_{i+\frac{1}{2}}^' + lr + r'$'
    plt.text(xv, y0+0.5, lrname, fontsize=12, ha='center')
    plt.text(-4, y0+0.3, '$x=0$', fontsize=12, ha='center')
    plt.text(+4, y0+0.3, '$x=L$', fontsize=12, ha='center')
    
    plt.text(0, y0-1.5, txt_name, fontsize=12, ha='center')
    return

# 设置字体为 Times New Roman
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif=['Times New Roman'])

# 设置图形大小和样式
plt.figure(figsize=(12, 5))

# 定义新的点坐标 (从 i=-6 到 i=6，基于 x_points)
x_points = np.array([-6, -5, -4, -2, -1, 0, 1, 2, 3, 4, 5, 6],dtype=np.float64)

y0 = 0.5
plot_cfd_line( x_points, y0, 'L' )

xv = 0.5*(x_points[5]+x_points[6])
plt.plot([xv, xv], [y0-0.5, y0+0.5], 'k--')  # 绘制垂直线

# 添加圆角矩形背景
rectangle_color = (150/255, 150/255, 200/255)  # RGB颜色

width = 4.4
height = 0.2

plot_rect(-2.2,y0-0.1,width,height,rectangle_color)
plot_label(y0,xv-0.3,'L','(a) Left-side reconstruction')

y1 = -2.0
plot_cfd_line( x_points, y1, 'R' )
plt.plot([xv, xv], [y1-0.5, y1+0.5], 'k--')  # 绘制垂直线

plot_rect(-1.2,y1-0.1,width,height,rectangle_color)
plot_label(y1,xv+0.3,'R','(b) Right-side reconstruction')


plt.axis('equal')
plt.axis('off')

plt.savefig('cfd.png', bbox_inches='tight', dpi=300)
plt.show()