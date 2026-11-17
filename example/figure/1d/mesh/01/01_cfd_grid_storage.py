"""
文件名: 01_cfd_grid_storage.py
功能: 绘制一维CFD基础网格与变量存储示意图
包含: 顶点中心存储、单元中心存储、边界条件处理、有限差分格式示意
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def setup_plot_style():
    """设置绘图样式"""
    plt.rcParams.update({
        'font.size': 11,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'figure.titlesize': 14,
        'figure.dpi': 100
    })

def plot_cfd_grid_storage():
    """绘制一维CFD网格和变量存储方式对比"""
    setup_plot_style()
    fig = plt.figure(figsize=(15, 10))
    
    # 创建网格
    n_cells = 5
    dx = 1.0
    x_vertices = np.linspace(0, n_cells*dx, n_cells + 1)
    x_centers = (x_vertices[:-1] + x_vertices[1:]) / 2
    
    # 1. 顶点中心存储
    ax1 = plt.subplot(2, 2, 1)
    ax1.set_title("Vertex-centered Storage", fontsize=14, fontweight='bold', pad=20)
    
    # 绘制网格线
    for x in x_vertices:
        ax1.axvline(x, color='gray', linestyle='-', alpha=0.5, linewidth=0.8)
    
    # 绘制顶点
    ax1.scatter(x_vertices, np.zeros_like(x_vertices), 
                s=120, color='red', zorder=5, edgecolors='black', linewidth=1.5)
    
    # 标记顶点
    for i, x in enumerate(x_vertices):
        if i == 0:
            label = f'Boundary\n$u_0$'
            color = "orange"
        elif i == len(x_vertices) - 1:
            label = f'Boundary\n$u_{i}$'
            color = "orange"
        else:
            label = f'Storage\n$u_{i}$'
            color = "yellow"
        
        ax1.text(x, -0.2, label, ha='center', fontsize=10, va='top',
                bbox=dict(boxstyle="round,pad=0.4", facecolor=color, alpha=0.8, edgecolor='black'))
    
    # 绘制控制体
    for i in range(len(x_vertices)-1):
        center = (x_vertices[i] + x_vertices[i+1]) / 2
        rect = ax1.add_patch(Rectangle((x_vertices[i], -0.08), dx, 0.16,
                              alpha=0.2, color='blue', 
                              edgecolor='blue', linewidth=1.5))
        ax1.text(center, 0.15, f'Cell {i}', ha='center', fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8))
    
    ax1.set_xlim(-0.5, n_cells*dx + 0.5)
    ax1.set_ylim(-0.35, 0.35)
    ax1.set_xlabel('Position x')
    ax1.set_ylabel('Variable Storage')
    ax1.text(0.5, 0.95, 'Variables stored at vertices', transform=ax1.transAxes,
            ha='center', fontsize=11, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow"))
    ax1.grid(True, alpha=0.3)
    
    # 2. 单元中心存储
    ax2 = plt.subplot(2, 2, 2)
    ax2.set_title("Cell-centered Storage", fontsize=14, fontweight='bold', pad=20)
    
    # 绘制网格线
    for x in x_vertices:
        ax2.axvline(x, color='gray', linestyle='-', alpha=0.5, linewidth=0.8)
    
    # 绘制单元中心
    ax2.scatter(x_centers, np.zeros_like(x_centers), 
                s=120, color='green', zorder=5, edgecolors='black', linewidth=1.5)
    
    # 标记变量
    for i, x in enumerate(x_centers):
        label = f'Storage\n$u_{i}$'
        ax2.text(x, -0.2, label, ha='center', fontsize=10, va='top',
                bbox=dict(boxstyle="round,pad=0.4", facecolor="lightgreen", alpha=0.8, edgecolor='black'))
    
    # 绘制控制体
    for i, x in enumerate(x_vertices[:-1]):
        rect = ax2.add_patch(Rectangle((x, -0.08), dx, 0.16,
                              alpha=0.2, color='orange', 
                              edgecolor='orange', linewidth=1.5))
        ax2.text(x + dx/2, 0.15, f'Cell {i}', ha='center', fontsize=11,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="peachpuff", alpha=0.8))
    
    # 标记边界
    ax2.axvline(x_vertices[0], color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax2.axvline(x_vertices[-1], color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax2.text(x_vertices[0], 0.25, 'Boundary', ha='center', color='red', fontsize=11, fontweight='bold')
    ax2.text(x_vertices[-1], 0.25, 'Boundary', ha='center', color='red', fontsize=11, fontweight='bold')
    
    ax2.set_xlim(-0.5, n_cells*dx + 0.5)
    ax2.set_ylim(-0.35, 0.35)
    ax2.set_xlabel('Position x')
    ax2.set_ylabel('Variable Storage')
    ax2.text(0.5, 0.95, 'Variables stored at cell centers', transform=ax2.transAxes,
            ha='center', fontsize=11, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow"))
    ax2.grid(True, alpha=0.3)
    
    # 3. 边界条件处理示意
    ax3 = plt.subplot(2, 2, 3)
    ax3.set_title("Boundary Conditions and Ghost Cells", fontsize=14, fontweight='bold', pad=20)
    
    # 扩展网格显示虚拟点
    x_extended = np.linspace(-dx, (n_cells+1)*dx, n_cells + 4)
    x_real = x_extended[1:-2]  # 真实计算区域
    
    # 绘制所有点
    ax3.scatter(x_extended, np.zeros_like(x_extended), s=80, color='gray', alpha=0.5)
    ax3.scatter(x_real, np.zeros_like(x_real), s=120, color='blue', 
                edgecolors='black', linewidth=1.5, zorder=5)
    
    # 标记边界
    ax3.axvline(0, color='red', linestyle='-', linewidth=3, alpha=0.8)
    ax3.axvline(n_cells*dx, color='red', linestyle='-', linewidth=3, alpha=0.8)
    
    # 填充虚拟点区域
    ax3.axvspan(-dx, 0, alpha=0.15, color='red', hatch='//')
    ax3.axvspan(n_cells*dx, (n_cells+1)*dx, alpha=0.15, color='red', hatch='//')
    
    # 标记点类型
    ax3.text(-dx/2, 0.15, 'Ghost Cell', ha='center', fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="pink", alpha=0.9, edgecolor='red'))
    ax3.text(n_cells*dx + dx/2, 0.15, 'Ghost Cell', ha='center', fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="pink", alpha=0.9, edgecolor='red'))
    
    # 标记计算区域
    ax3.axvspan(0, n_cells*dx, alpha=0.1, color='green')
    ax3.text(n_cells*dx/2, -0.2, 'Computational Domain', ha='center', fontsize=12,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="lightgreen", alpha=0.8))
    
    # 添加箭头示意边界条件
    ax3.annotate('BC: u=0', xy=(0, 0), xytext=(-1.2*dx, 0.25),
                arrowprops=dict(arrowstyle="->", color='darkred', lw=2),
                ha='center', fontsize=10, color='darkred',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))
    
    ax3.annotate('BC: ∂u/∂x=0', xy=(n_cells*dx, 0), 
                xytext=(n_cells*dx + 1.2*dx, 0.25),
                arrowprops=dict(arrowstyle="->", color='darkred', lw=2),
                ha='center', fontsize=10, color='darkred',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))
    
    ax3.set_xlim(-1.8*dx, (n_cells+1.8)*dx)
    ax3.set_ylim(-0.3, 0.4)
    ax3.set_xlabel('Position x')
    ax3.set_ylabel('Domain')
    ax3.grid(True, alpha=0.3)
    
    # 4. 数值格式示意
    ax4 = plt.subplot(2, 2, 4)
    ax4.set_title("Finite Difference Schemes", fontsize=14, fontweight='bold', pad=20)
    
    # 创建示例数据
    x = np.linspace(0, 10, 50)
    u = np.sin(x * 0.8) * np.exp(-0.1*x)
    
    # 选取几个点
    i = 25
    ax4.plot(x, u, 'b-', linewidth=3, alpha=0.5, label='Exact Solution')
    ax4.scatter(x[i], u[i], s=200, color='red', zorder=5, 
                edgecolors='black', linewidth=1.5, label='Current point $u_i$')
    ax4.scatter(x[i-1], u[i-1], s=150, color='green', zorder=4, 
                edgecolors='black', linewidth=1, label='Upstream $u_{i-1}$')
    ax4.scatter(x[i+1], u[i+1], s=150, color='orange', zorder=4, 
                edgecolors='black', linewidth=1, label='Downstream $u_{i+1}$')
    
    # 绘制差分示意线
    ax4.plot([x[i-1], x[i+1]], [u[i-1], u[i+1]], 'k--', alpha=0.5, linewidth=1.5)
    
    # 标注差分格式
    ax4.annotate('Central Difference\n(2nd order)', xy=(x[i], u[i]), xytext=(x[i], u[i]+0.4),
                arrowprops=dict(arrowstyle="->", color='blue', lw=2),
                ha='center', fontsize=11, color='blue',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))
    
    ax4.annotate('Upwind Scheme\n(1st order)', xy=(x[i], u[i]), xytext=(x[i]-2.5, u[i]-0.3),
                arrowprops=dict(arrowstyle="->", color='red', lw=2),
                ha='center', fontsize=11, color='red',
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.9))
    
    ax4.set_xlabel('Position x')
    ax4.set_ylabel('Variable u')
    ax4.legend(loc='upper right', fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    # 添加整体标题
    plt.suptitle('1D CFD Grid and Variable Storage Illustration', fontsize=16, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('01_cfd_grid_storage.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()

if __name__ == "__main__":
    plot_cfd_grid_storage()