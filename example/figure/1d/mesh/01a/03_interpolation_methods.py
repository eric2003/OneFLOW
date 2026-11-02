"""
文件名: 03_interpolation_methods.py
功能: 绘制不同插值方法对比
包含: 线性插值、二次插值、三次样条插值的比较
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

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

def plot_interpolation_methods():
    """绘制不同插值方法对比"""
    setup_plot_style()
    
    # 创建粗网格和细网格
    x_coarse = np.linspace(0, 10, 6)
    u_coarse = np.sin(x_coarse * 0.8)
    
    x_fine = np.linspace(0, 10, 100)
    u_exact = np.sin(x_fine * 0.8)
    
    # 不同插值方法
    # 线性插值
    f_linear = interpolate.interp1d(x_coarse, u_coarse, kind='linear')
    u_linear = f_linear(x_fine)
    
    # 二次插值
    f_quadratic = interpolate.interp1d(x_coarse, u_coarse, kind='quadratic')
    u_quadratic = f_quadratic(x_fine)
    
    # 三次样条插值
    f_cubic = interpolate.CubicSpline(x_coarse, u_coarse)
    u_cubic = f_cubic(x_fine)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. 线性插值
    ax = axes[0, 0]
    ax.plot(x_fine, u_exact, 'k-', alpha=0.3, linewidth=3, label='Exact Solution')
    ax.plot(x_fine, u_linear, 'r--', linewidth=2, label='Linear Interpolation')
    ax.scatter(x_coarse, u_coarse, s=100, color='blue', zorder=5, label='Known Points')
    ax.set_title('Linear Interpolation', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. 二次插值
    ax = axes[0, 1]
    ax.plot(x_fine, u_exact, 'k-', alpha=0.3, linewidth=3, label='Exact Solution')
    ax.plot(x_fine, u_quadratic, 'g--', linewidth=2, label='Quadratic Interpolation')
    ax.scatter(x_coarse, u_coarse, s=100, color='blue', zorder=5, label='Known Points')
    ax.set_title('Quadratic Interpolation', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. 三次样条插值
    ax = axes[1, 0]
    ax.plot(x_fine, u_exact, 'k-', alpha=0.3, linewidth=3, label='Exact Solution')
    ax.plot(x_fine, u_cubic, 'b--', linewidth=2, label='Cubic Spline')
    ax.scatter(x_coarse, u_coarse, s=100, color='blue', zorder=5, label='Known Points')
    ax.set_title('Cubic Spline Interpolation', fontsize=12, fontweight='bold')
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. 误差比较
    ax = axes[1, 1]
    errors = {
        'Linear': np.abs(u_linear - u_exact),
        'Quadratic': np.abs(u_quadratic - u_exact),
        'Cubic Spline': np.abs(u_cubic - u_exact)
    }
    
    x_pos = np.arange(len(errors))
    mean_errors = [np.mean(err) for err in errors.values()]
    max_errors = [np.max(err) for err in errors.values()]
    
    width = 0.35
    ax.bar(x_pos - width/2, mean_errors, width, label='Mean Error', color='skyblue')
    ax.bar(x_pos + width/2, max_errors, width, label='Max Error', color='salmon')
    
    ax.set_xlabel('Interpolation Method')
    ax.set_ylabel('Error')
    ax.set_title('Interpolation Error Comparison', fontsize=12, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(list(errors.keys()))
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # 添加数值标签
    for i, (mean_err, max_err) in enumerate(zip(mean_errors, max_errors)):
        ax.text(i - width/2, mean_err + 0.001, f'{mean_err:.3f}', 
                ha='center', va='bottom', fontsize=9)
        ax.text(i + width/2, max_err + 0.001, f'{max_err:.3f}', 
                ha='center', va='bottom', fontsize=9)
    
    plt.suptitle('Comparison of Interpolation Methods for CFD', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('03_interpolation_methods.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    plot_interpolation_methods()