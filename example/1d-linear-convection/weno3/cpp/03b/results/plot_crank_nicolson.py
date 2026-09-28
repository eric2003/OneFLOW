# ==================== scripts/plot_crank_nicolson.py (修复中文乱码) ====================
#!/usr/bin/env python3
"""
Crank-Nicolson与解析解对比绘图脚本（修复中文乱码）
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# ===================================================================
# 1. 中文显示设置（解决乱码问题）
# ===================================================================
def setup_chinese_font():
    """设置中文字体支持"""
    import matplotlib
    
    # 尝试多种中文字体方案
    font_options = [
        # Windows系统字体
        "Microsoft YaHei",        # 微软雅黑
        "SimHei",                 # 黑体
        "SimSun",                 # 宋体
        "NSimSun",                # 新宋体
        "FangSong",               # 仿宋
        "KaiTi",                  # 楷体
        
        # 通用字体回退
        "DejaVu Sans",
        "Arial",
        "sans-serif"
    ]
    
    # 尝试设置中文字体
    for font_name in font_options:
        try:
            matplotlib.rcParams['font.sans-serif'] = [font_name]
            matplotlib.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题
            print(f"✓ 使用字体: {font_name}")
            return True
        except:
            continue
    
    print("⚠ 警告: 未找到合适的中文字体，将使用默认字体")
    return False

# 初始化中文字体
setup_chinese_font()

def read_solution_file(filename):
    """读取结果文件"""
    try:
        data = np.loadtxt(filename, comments='#')
        if data.shape[1] == 2:
            x = data[:, 0]
            u = data[:, 1]
            return x, u
        else:
            print(f"Warning: File {filename} has wrong format, expected 2 columns, got {data.shape[1]}")
            return np.array([]), np.array([])
    except Exception as e:
        print(f"Failed to read {filename}: {e}")
        return np.array([]), np.array([])

def plot_comparison(cn_file="crank_nicolson_results.txt", 
                   analytical_file="analytical_results.txt",
                   output_file="crank_nicolson_comparison.png"):
    """绘制对比图"""
    
    # 检查文件是否存在
    if not os.path.exists(cn_file):
        print(f"Error: Cannot find {cn_file}")
        print("Please make sure Crank-Nicolson simulation has been run")
        return False
    
    if not os.path.exists(analytical_file):
        print(f"Error: Cannot find {analytical_file}")
        print("Please make sure analytical solution file exists")
        return False
    
    print("Reading data...")
    x_cn, u_cn = read_solution_file(cn_file)
    x_ana, u_ana = read_solution_file(analytical_file)
    
    if len(x_cn) == 0 or len(u_cn) == 0:
        print("Cannot read Crank-Nicolson results")
        return False
    
    if len(x_ana) == 0 or len(u_ana) == 0:
        print("Cannot read analytical solution")
        return False
    
    # 检查数据长度是否匹配
    if len(x_cn) != len(x_ana):
        print(f"Warning: Data length mismatch (CN: {len(x_cn)}, Analytical: {len(x_ana)})")
        print("Using Crank-Nicolson x-coordinates")
    
    print(f"Number of data points: {len(x_cn)}")
    
    # 创建图形
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    
    # 第一个子图：数值解与解析解对比
    ax1 = axes[0]
    ax1.plot(x_cn, u_cn, 'bo-', linewidth=1, markersize=3, 
             markerfacecolor='none', label='Crank-Nicolson (Numerical)')
    ax1.plot(x_ana, u_ana, 'r-', linewidth=2, label='Theoretical (Analytical)')
    
    ax1.set_xlabel('x-coordinate', fontsize=12)
    ax1.set_ylabel('u(x)', fontsize=12)
    ax1.set_title('Crank-Nicolson Method: Numerical vs Analytical Solution', 
                  fontsize=14, fontweight='bold')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # 计算并显示误差统计
    error = u_cn - u_ana
    l1_error = np.mean(np.abs(error))
    l2_error = np.sqrt(np.mean(error**2))
    linf_error = np.max(np.abs(error))
    
    stats_text = f'L1 Error: {l1_error:.3e}\nL2 Error: {l2_error:.3e}\nL∞ Error: {linf_error:.3e}'
    ax1.text(0.02, 0.98, stats_text,
             transform=ax1.transAxes,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
             fontsize=10)
    
    # 第二个子图：误差分布
    ax2 = axes[1]
    ax2.plot(x_cn, error, 'g-', linewidth=1.5, label='Error (Numerical - Analytical)')
    ax2.axhline(y=0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    
    ax2.set_xlabel('x-coordinate', fontsize=12)
    ax2.set_ylabel('Error', fontsize=12)
    ax2.set_title('Error Distribution', fontsize=14, fontweight='bold')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)
    
    # 在误差图中添加误差范围
    error_mean = np.mean(error)
    error_std = np.std(error)
    ax2.axhline(y=error_mean, color='r', linestyle=':', linewidth=1, alpha=0.7, 
                label=f'Mean: {error_mean:.3e}')
    ax2.fill_between(x_cn, error_mean - error_std, error_mean + error_std, 
                     alpha=0.2, color='gray', label=f'±1 Std Dev')
    
    ax2.legend(loc='upper right', fontsize=9)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图像
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved as: {output_file}")
    
    # 显示图像
    plt.show()
    
    # 打印详细统计信息
    print("\n" + "="*50)
    print("Error Statistics Report")
    print("="*50)
    print(f"Number of data points: {len(x_cn)}")
    print(f"L1 Error (MAE): {l1_error:.6e}")
    print(f"L2 Error (RMSE): {l2_error:.6e}")
    print(f"L∞ Error (Max AE): {linf_error:.6e}")
    print(f"Error Mean: {error_mean:.6e}")
    print(f"Error Std Dev: {error_std:.6e}")
    print(f"Error Range: [{np.min(error):.6e}, {np.max(error):.6e}]")
    
    # 计算相对误差
    mask = np.abs(u_ana) > 1e-12  # Avoid division by zero
    if np.any(mask):
        relative_error = np.abs(error[mask] / u_ana[mask])
        max_rel_error = np.max(relative_error)
        avg_rel_error = np.mean(relative_error)
        print(f"Max Relative Error: {max_rel_error:.6e}")
        print(f"Avg Relative Error: {avg_rel_error:.6e}")
    
    return True

def save_error_data(cn_file="crank_nicolson_results.txt", 
                   analytical_file="analytical_results.txt"):
    """保存误差数据到文件"""
    x_cn, u_cn = read_solution_file(cn_file)
    x_ana, u_ana = read_solution_file(analytical_file)
    
    if len(x_cn) == 0 or len(x_ana) == 0:
        return False
    
    error = u_cn - u_ana
    
    # 保存详细的误差数据
    with open("error_analysis.txt", "w", encoding='utf-8') as f:
        f.write("# Crank-Nicolson Error Analysis\n")
        f.write("# ============================\n")
        f.write("# Columns: 1-x-coordinate, 2-CN-numerical, 3-Analytical, 4-Absolute-Error, 5-Relative-Error\n")
        f.write("#\n")
        
        for i in range(len(x_cn)):
            abs_error = abs(error[i])
            rel_error = abs_error / abs(u_ana[i]) if abs(u_ana[i]) > 1e-12 else 0.0
            f.write(f"{x_cn[i]:12.6f} {u_cn[i]:12.6f} {u_ana[i]:12.6f} "
                   f"{error[i]:12.6e} {rel_error:12.6e}\n")
    
    print("Detailed error analysis saved to: error_analysis.txt")
    return True

def main():
    """主函数"""
    print("="*60)
    print("Crank-Nicolson Results Visualization")
    print("="*60)
    
    # 使用命令行参数（如果提供）
    if len(sys.argv) >= 3:
        cn_file = sys.argv[1]
        analytical_file = sys.argv[2]
        output_file = sys.argv[3] if len(sys.argv) >= 4 else "crank_nicolson_comparison.png"
    else:
        cn_file = "crank_nicolson_results.txt"
        analytical_file = "analytical_results.txt"
        output_file = "crank_nicolson_comparison.png"
    
    print(f"Files:")
    print(f"  Crank-Nicolson results: {cn_file}")
    print(f"  Analytical solution: {analytical_file}")
    print(f"  Output image: {output_file}")
    print()
    
    # 绘制对比图
    success = plot_comparison(cn_file, analytical_file, output_file)
    
    if success:
        # 保存详细的误差数据
        save_error_data(cn_file, analytical_file)
        
        print("\n" + "="*60)
        print("✅ Visualization completed!")
        print("="*60)
        print("\nGenerated files:")
        print(f"  1. {output_file} - Comparison plot")
        print(f"  2. error_analysis.txt - Detailed error data")
        print("\nNext time you can run directly:")
        print(f"  python {sys.argv[0]} [CN-file] [analytical-file] [output-image]")
        print("\nExample:")
        print(f"  python {sys.argv[0]} my_cn_results.txt my_analytical.txt my_plot.png")
    else:
        print("\n❌ Visualization failed, please check input files")

if __name__ == "__main__":
    main()