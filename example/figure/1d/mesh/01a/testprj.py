import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

def plot_cfd_grid_storage():
    """绘制一维CFD网格和变量存储方式对比"""
    fig = plt.figure(figsize=(14, 8))
    
    # 创建网格
    n_cells = 5
    dx = 1.0
    x_vertices = np.linspace(0, n_cells*dx, n_cells + 1)
    x_centers = (x_vertices[:-1] + x_vertices[1:]) / 2
    
    # 1. 顶点中心存储
    ax1 = plt.subplot(2, 2, 1)
    ax1.set_title("顶点中心存储 (Vertex-centered)", fontsize=12, fontweight='bold')
    
    # 绘制网格线
    for x in x_vertices:
        ax1.axvline(x, color='gray', linestyle='-', alpha=0.5, linewidth=0.8)
    
    # 绘制顶点
    ax1.scatter(x_vertices, np.zeros_like(x_vertices), 
                s=100, color='red', zorder=5, label='变量存储点')
    
    # 标记顶点
    for i, x in enumerate(x_vertices):
        ax1.text(x, -0.15, f'$u_{i}$', ha='center', fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # 绘制控制体
    for i in range(len(x_vertices)-1):
        center = (x_vertices[i] + x_vertices[i+1]) / 2
        ax1.add_patch(Rectangle((x_vertices[i], -0.05), dx, 0.1,
                              alpha=0.2, color='blue', label='控制体' if i==0 else None))
        ax1.text(center, 0.1, f'Cell {i}', ha='center', fontsize=9)
    
    ax1.set_xlim(-0.5, n_cells*dx + 0.5)
    ax1.set_ylim(-0.3, 0.3)
    ax1.set_xlabel('x')
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # 2. 单元中心存储
    ax2 = plt.subplot(2, 2, 2)
    ax2.set_title("单元中心存储 (Cell-centered)", fontsize=12, fontweight='bold')
    
    # 绘制网格线
    for x in x_vertices:
        ax2.axvline(x, color='gray', linestyle='-', alpha=0.5, linewidth=0.8)
    
    # 绘制单元中心
    ax2.scatter(x_centers, np.zeros_like(x_centers), 
                s=100, color='green', zorder=5, label='变量存储点')
    
    # 标记变量
    for i, x in enumerate(x_centers):
        ax2.text(x, -0.15, f'$u_{i}$', ha='center', fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.7))
    
    # 绘制控制体
    for i, x in enumerate(x_vertices[:-1]):
        ax2.add_patch(Rectangle((x, -0.05), dx, 0.1,
                              alpha=0.2, color='orange', label='控制体' if i==0 else None))
        ax2.text(x + dx/2, 0.1, f'Cell {i}', ha='center', fontsize=9)
    
    ax2.set_xlim(-0.5, n_cells*dx + 0.5)
    ax2.set_ylim(-0.3, 0.3)
    ax2.set_xlabel('x')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # 3. 边界条件处理示意
    ax3 = plt.subplot(2, 2, 3)
    ax3.set_title("边界条件与虚拟点", fontsize=12, fontweight='bold')
    
    # 扩展网格显示虚拟点
    x_extended = np.linspace(-dx, (n_cells+1)*dx, n_cells + 4)
    x_real = x_extended[1:-2]  # 真实计算区域
    
    # 绘制所有点
    ax3.scatter(x_extended, np.zeros_like(x_extended), s=50, color='gray', alpha=0.5)
    ax3.scatter(x_real, np.zeros_like(x_real), s=100, color='blue', label='计算点')
    
    # 标记区域
    ax3.axvline(0, color='red', linestyle='--', linewidth=2, label='左边界')
    ax3.axvline(n_cells*dx, color='red', linestyle='--', linewidth=2, label='右边界')
    ax3.axvspan(-dx, 0, alpha=0.1, color='red', label='虚拟点区域')
    ax3.axvspan(n_cells*dx, (n_cells+1)*dx, alpha=0.1, color='red')
    
    # 标记点类型
    ax3.text(-dx/2, 0.1, '虚拟点', ha='center', fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="pink", alpha=0.7))
    ax3.text(n_cells*dx + dx/2, 0.1, '虚拟点', ha='center', fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="pink", alpha=0.7))
    
    ax3.set_xlim(-1.5*dx, (n_cells+1.5)*dx)
    ax3.set_ylim(-0.2, 0.3)
    ax3.set_xlabel('x')
    ax3.legend(loc='upper right')
    ax3.grid(True, alpha=0.3)
    
    # 4. 数值格式示意
    ax4 = plt.subplot(2, 2, 4)
    ax4.set_title("有限差分格式示意", fontsize=12, fontweight='bold')
    
    # 创建示例数据
    x = np.linspace(0, 10, 50)
    u = np.sin(x * 0.8) * np.exp(-0.1*x)
    
    # 选取几个点
    i = 25
    ax4.plot(x, u, 'b-', linewidth=2, label='真实解')
    ax4.scatter(x[i], u[i], s=150, color='red', zorder=5, label='计算点 $u_i$')
    ax4.scatter(x[i-1], u[i-1], s=100, color='green', zorder=4, label='$u_{i-1}$')
    ax4.scatter(x[i+1], u[i+1], s=100, color='orange', zorder=4, label='$u_{i+1}$')
    
    # 绘制差分示意
    ax4.plot([x[i-1], x[i+1]], [u[i-1], u[i+1]], 'k--', alpha=0.5)
    ax4.annotate('中心差分', xy=(x[i], u[i]), xytext=(x[i], u[i]+0.3),
                arrowprops=dict(arrowstyle="->", color='black'),
                ha='center', fontsize=10)
    
    ax4.annotate('迎风格式\n使用上游点', xy=(x[i], u[i]), xytext=(x[i]-2, u[i]-0.3),
                arrowprops=dict(arrowstyle="->", color='red'),
                ha='center', fontsize=9)
    
    ax4.set_xlabel('x')
    ax4.set_ylabel('u')
    ax4.legend(loc='upper right')
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('cfd_grid_illustration.png', dpi=300, bbox_inches='tight')
    plt.show()

# 运行绘图
plot_cfd_grid_storage()