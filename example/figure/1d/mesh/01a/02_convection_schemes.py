"""
文件名: 02_convection_schemes.py
功能: 绘制一维对流方程不同数值格式对比
包含: 迎风格式、中心差分格式、QUICK格式的比较
"""

import numpy as np
import matplotlib.pyplot as plt

def setup_plot_style():
    """设置绘图样式"""
    plt.rcParams.update({
        'font.size': 11,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'figure.titlesize': 14
    })

def plot_convection_schemes():
    """绘制不同对流格式示意图"""
    setup_plot_style()
    
    # 创建数据
    x = np.linspace(0, 2*np.pi, 100)
    u_exact = np.sin(x)
    
    # 添加数值噪声模拟数值解
    np.random.seed(42)
    u_upwind = u_exact + 0.1*np.random.randn(len(x))  # 迎风格式（有耗散）
    u_central = u_exact + 0.05*np.random.randn(len(x))  # 中心差分（有振荡）
    u_quick = u_exact + 0.02*np.random.randn(len(x))   # QUICK格式
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. 迎风格式
    ax = axes[0, 0]
    ax.plot(x, u_exact, 'k-', linewidth=3, alpha=0.7, label='Exact Solution')
    ax.plot(x, u_upwind, 'r--', linewidth=2, marker='o', markersize=4, 
            markevery=5, label='Upwind Scheme')
    ax.fill_between(x, u_exact-0.15, u_exact+0.15, alpha=0.1, color='gray')
    ax.set_title('Upwind Scheme', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.annotate('Numerical\nDissipation', xy=(3, 0), xytext=(4, -0.8),
                arrowprops=dict(arrowstyle="->", color='red'),
                fontsize=10, color='red')
    
    # 2. 中心差分格式
    ax = axes[0, 1]
    ax.plot(x, u_exact, 'k-', linewidth=3, alpha=0.7, label='Exact Solution')
    ax.plot(x, u_central, 'b--', linewidth=2, marker='s', markersize=4,
            markevery=5, label='Central Difference')
    ax.set_title('Central Difference Scheme', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.annotate('Numerical\nOscillation', xy=(1.5, 0.5), xytext=(0.5, 0.8),
                arrowprops=dict(arrowstyle="->", color='blue'),
                fontsize=10, color='blue')
    
    # 3. QUICK格式
    ax = axes[1, 0]
    ax.plot(x, u_exact, 'k-', linewidth=3, alpha=0.7, label='Exact Solution')
    ax.plot(x, u_quick, 'g--', linewidth=2, marker='^', markersize=4,
            markevery=5, label='QUICK Scheme')
    ax.set_title('QUICK Scheme (3rd Order)', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. 格式示意
    ax = axes[1, 1]
    ax.set_title('Grid Point Dependencies', fontsize=12, fontweight='bold')
    
    # 绘制网格点
    points_x = np.array([0, 1, 2, 3, 4])
    points_y = np.zeros_like(points_x)
    
    ax.scatter(points_x, points_y, s=200, color='gray')
    
    # 标记点
    labels = ['$u_{i-2}$', '$u_{i-1}$', '$u_i$', '$u_{i+1}$', '$u_{i+2}$']
    for i, (x_pos, label) in enumerate(zip(points_x, labels)):
        ax.text(x_pos, 0.1, label, ha='center', fontsize=12, fontweight='bold')
    
    # 迎风格式依赖
    ax.annotate('', xy=(2, 0), xytext=(1, 0),
                arrowprops=dict(arrowstyle='<-', color='red', lw=2))
    ax.text(1.5, -0.15, 'Upwind', ha='center', color='red', fontweight='bold')
    
    # 中心差分依赖
    ax.annotate('', xy=(2, 0), xytext=(1, -0.05),
                arrowprops=dict(arrowstyle='<-', color='blue', lw=2))
    ax.annotate('', xy=(2, 0), xytext=(3, -0.05),
                arrowprops=dict(arrowstyle='<-', color='blue', lw=2))
    ax.text(2, -0.2, 'Central', ha='center', color='blue', fontweight='bold')
    
    # QUICK格式依赖
    ax.annotate('', xy=(2, 0), xytext=(0, 0.05),
                arrowprops=dict(arrowstyle='<-', color='green', lw=2, alpha=0.7))
    ax.annotate('', xy=(2, 0), xytext=(1, 0.05),
                arrowprops=dict(arrowstyle='<-', color='green', lw=2, alpha=0.7))
    ax.annotate('', xy=(2, 0), xytext=(3, 0.05),
                arrowprops=dict(arrowstyle='<-', color='green', lw=2, alpha=0.7))
    ax.text(2, 0.2, 'QUICK', ha='center', color='green', fontweight='bold')
    
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-0.3, 0.3)
    ax.set_xlabel('Grid Point Index')
    ax.grid(True, alpha=0.3)
    ax.set_yticks([])
    
    plt.suptitle('Comparison of Convection Schemes for 1D CFD', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('02_convection_schemes.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    plot_convection_schemes()