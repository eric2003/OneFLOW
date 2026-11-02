import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 设置中文字体
rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
rcParams['axes.unicode_minus'] = False

def load_and_plot():
    """加载并绘制精确参数的结果"""
    data = np.loadtxt('precise_wave.dat', comments='#')
    x = data[:, 0]
    u = data[:, 1]
    
    # 计算参数
    delta = 0.005
    alpha = 10.0
    beta = np.log(2) / (36 * delta**2)
    
    # 创建图形
    fig = plt.figure(figsize=(14, 10))
    
    # 1. 整体视图
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(x, u, 'b-', linewidth=2, label='初始条件 u(x,0)')
    ax1.fill_between(x, 0, u, alpha=0.2, color='blue')
    
    # 标记区域
    regions = [(-0.8, -0.6, f'高斯脉冲\nβ={beta:.1f}'),
               (-0.4, -0.2, '方波'),
               (0.0, 0.2, '三角波'),
               (0.4, 0.6, f'半椭圆\nα={alpha}')]
    
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4']
    
    for (xmin, xmax, label), color in zip(regions, colors):
        ax1.axvspan(xmin, xmax, alpha=0.1, color=color)
        ax1.text((xmin+xmax)/2, 1.05, label, 
                ha='center', va='bottom', fontsize=10,
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('u(x,0)')
    ax1.set_title('复杂波形初始条件（精确参数）\n' +
                 f'a=0.5, z=-0.7, δ={delta}, α={alpha}, β=log2/(36δ²)={beta:.2f}')
    ax1.set_xlim(-1.0, 1.0)
    ax1.set_ylim(-0.05, 1.15)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    
    # 2. 区域1详细视图（高斯脉冲）
    ax2 = plt.subplot(2, 2, 2)
    mask1 = (x >= -0.8) & (x <= -0.6)
    x1 = x[mask1]
    u1 = u[mask1]
    
    ax2.plot(x1, u1, 'b-', linewidth=2, label='组合脉冲')
    
    # 计算三个高斯分量
    z = -0.7
    g1 = np.exp(-beta * (x1 - (z - delta))**2) / 6.0
    g2 = np.exp(-beta * (x1 - (z + delta))**2) / 6.0
    g3 = 4.0 * np.exp(-beta * (x1 - z)**2) / 6.0
    
    ax2.plot(x1, g1, 'r--', alpha=0.6, linewidth=1, label='G(x,β,z-δ)/6')
    ax2.plot(x1, g2, 'g--', alpha=0.6, linewidth=1, label='G(x,β,z+δ)/6')
    ax2.plot(x1, g3, 'm--', alpha=0.6, linewidth=1, label='4G(x,β,z)/6')
    
    # 标记关键点
    ax2.axvline(x=z, color='k', linestyle=':', alpha=0.5, label=f'中心 z={z}')
    ax2.axvline(x=z-delta, color='r', linestyle=':', alpha=0.3)
    ax2.axvline(x=z+delta, color='g', linestyle=':', alpha=0.3)
    
    ax2.set_xlabel('x')
    ax2.set_ylabel('u(x,0)')
    ax2.set_title('区域1：高斯脉冲组合（详细）')
    ax2.set_xlim(-0.8, -0.6)
    ax2.set_ylim(0.85, 1.05)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=9)
    
    # 3. 区域4详细视图（半椭圆）
    ax3 = plt.subplot(2, 2, 3)
    mask4 = (x >= 0.4) & (x <= 0.6)
    x4 = x[mask4]
    u4 = u[mask4]
    
    ax3.plot(x4, u4, 'b-', linewidth=2, label='组合椭圆')
    
    # 计算三个半椭圆分量
    a = 0.5
    # F(x,α,a) = √max(1 - α²(x-a)², 0)
    def ellipse_func(x_val, center):
        arg = 1.0 - (alpha * (x_val - center))**2
        return np.sqrt(np.maximum(arg, 0))
    
    f1 = ellipse_func(x4, a - delta) / 6.0
    f2 = ellipse_func(x4, a + delta) / 6.0
    f3 = 4.0 * ellipse_func(x4, a) / 6.0
    
    ax3.plot(x4, f1, 'r--', alpha=0.6, linewidth=1, label='F(x,α,a-δ)/6')
    ax3.plot(x4, f2, 'g--', alpha=0.6, linewidth=1, label='F(x,α,a+δ)/6')
    ax3.plot(x4, f3, 'm--', alpha=0.6, linewidth=1, label='4F(x,α,a)/6')
    
    ax3.axvline(x=a, color='k', linestyle=':', alpha=0.5, label=f'中心 a={a}')
    ax3.axvline(x=a-delta, color='r', linestyle=':', alpha=0.3)
    ax3.axvline(x=a+delta, color='g', linestyle=':', alpha=0.3)
    
    ax3.set_xlabel('x')
    ax3.set_ylabel('u(x,0)')
    ax3.set_title('区域4：半椭圆组合（详细）')
    ax3.set_xlim(0.4, 0.6)
    ax3.set_ylim(0.85, 1.05)
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc='upper left', fontsize=9)
    
    # 4. 所有区域的统计信息
    ax4 = plt.subplot(2, 2, 4)
    ax4.axis('off')
    
    # 计算各区域统计
    stats_text = "各区域统计信息:\n\n"
    
    for (xmin, xmax, label), color in zip(regions, colors):
        mask = (x >= xmin) & (x <= xmax)
        if np.any(mask):
            u_region = u[mask]
            max_val = np.max(u_region)
            mean_val = np.mean(u_region)
            area = np.trapz(u_region, x[mask])
            
            stats_text += f"{label}:\n"
            stats_text += f"  范围: [{xmin}, {xmax}]\n"
            stats_text += f"  最大值: {max_val:.4f}\n"
            stats_text += f"  平均值: {mean_val:.4f}\n"
            stats_text += f"  面积: {area:.4f}\n\n"
    
    # 总体统计
    total_area = np.trapz(u, x)
    total_variation = np.sum(np.abs(np.diff(u)))
    
    stats_text += f"总体统计:\n"
    stats_text += f"总定义域: [{x.min():.1f}, {x.max():.1f}]\n"
    stats_text += f"总点数: {len(x)}\n"
    stats_text += f"总变差: {total_variation:.4f}\n"
    stats_text += f"总面积: {total_area:.4f}\n"
    stats_text += f"全局最大值: {np.max(u):.4f}\n"
    stats_text += f"全局最小值: {np.min(u):.4f}"
    
    ax4.text(0.05, 0.95, stats_text, fontsize=10, 
             verticalalignment='top',
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    
    ax4.set_title('统计信息')
    
    plt.tight_layout()
    plt.savefig('precise_wave_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    return x, u, beta, delta, alpha

def compare_with_expected():
    """与预期值比较"""
    print("=== 关键点验证 ===")
    
    # 预期值（从Fortran输出获得）
    expected_values = {
        'Region 1 center (x=-0.7)': 0.993643362556305,
        'Region 4 center (x=0.5)': 0.999583072590636,
        'Square wave (x=-0.3)': 1.0,
        'Triangle peak (x=0.1)': 1.0
    }
    
    data = np.loadtxt('precise_wave.dat', comments='#')
    x = data[:, 0]
    u = data[:, 1]
    
    for desc, expected in expected_values.items():
        if 'x=' in desc:
            # 提取x坐标
            import re
            match = re.search(r'x=([-+]?\d*\.?\d+)', desc)
            if match:
                x_val = float(match.group(1))
                idx = np.argmin(np.abs(x - x_val))
                actual = u[idx]
                error = abs(actual - expected)
                print(f"{desc}:")
                print(f"  预期: {expected:.15f}")
                print(f"  实际: {actual:.15f}")
                print(f"  误差: {error:.2e}")
                if error < 1e-10:
                    print("  ✓ 通过")
                else:
                    print("  ✗ 未通过")
                print()

def analyze_gaussian_shape():
    """分析高斯脉冲形状"""
    print("\n=== 高斯脉冲形状分析 ===")
    
    delta = 0.005
    beta = np.log(2) / (36 * delta**2)
    
    print(f"参数: δ = {delta}, β = {beta:.4f}")
    print()
    
    # 计算半高宽
    fwhm = 2 * np.sqrt(np.log(2) / beta)
    print(f"1. 半高宽 (FWHM):")
    print(f"   FWHM = 2√(ln2/β) = {fwhm:.6f}")
    print(f"   在网格上的点数 (dx=0.00125): {fwhm/0.00125:.1f}")
    print()
    
    # 计算标准差
    sigma = 1 / np.sqrt(2 * beta)
    print(f"2. 标准差:")
    print(f"   σ = 1/√(2β) = {sigma:.6f}")
    print(f"   3σ宽度: {6*sigma:.6f}")
    print()
    
    # 计算在几个关键点的值
    z = -0.7
    points = [z, z + delta, z + 2*delta, z + fwhm/2]
    
    print(f"3. 关键点的值:")
    for point in points:
        value = np.exp(-beta * (point - z)**2)
        print(f"   x = {point:.4f}: G = {value:.6f}")
    
    # 检查是否像单个脉冲
    print(f"\n4. 脉冲重叠分析:")
    print(f"   在 x = z ± δ 处: G = {np.exp(-beta * delta**2):.6f}")
    print(f"   在 x = z ± 2δ 处: G = {np.exp(-beta * (2*delta)**2):.6f}")
    print(f"   由于 δ 很小 (0.005) 且 β 很大 (770.16)，")
    print(f"   三个高斯脉冲几乎完全重叠，看起来像单个脉冲。")

if __name__ == "__main__":
    print("=== 精确复杂波形分析 ===\n")
    x, u, beta, delta, alpha = load_and_plot()
    compare_with_expected()
    analyze_gaussian_shape()
    
    # 保存参数文件
    with open('wave_parameters.txt', 'w') as f:
        f.write("=== 复杂波形精确参数 ===\n\n")
        f.write(f"数学公式:\n")
        f.write(f"u₀(x) = (1/6)[G(x,β,z-δ) + G(x,β,z+δ) + 4G(x,β,z)], -0.8 ≤ x ≤ -0.6\n")
        f.write(f"        1,                                           -0.4 ≤ x ≤ -0.2\n")
        f.write(f"        1 - |10(x-0.1)|,                              0 ≤ x ≤ 0.2\n")
        f.write(f"        (1/6)[F(x,α,a-δ) + F(x,α,a+δ) + 4F(x,α,a)],   0.4 ≤ x ≤ 0.6\n")
        f.write(f"        0,                                           otherwise\n\n")
        f.write(f"其中:\n")
        f.write(f"  G(x,β,z) = exp(-β(x-z)²)\n")
        f.write(f"  F(x,α,a) = √max(1 - α²(x-a)², 0)\n\n")
        f.write(f"参数值:\n")
        f.write(f"  a = {0.5}        (区域4中心)\n")
        f.write(f"  z = {-0.7}       (区域1中心)\n")
        f.write(f"  δ = {delta}     (偏移量)\n")
        f.write(f"  α = {alpha}        (椭圆参数)\n")
        f.write(f"  β = {beta:.6f}  (高斯参数, = log2/(36δ²))\n\n")
        f.write(f"高斯脉冲性质:\n")
        fwhm = 2 * np.sqrt(np.log(2) / beta)
        sigma = 1 / np.sqrt(2 * beta)
        f.write(f"  标准差 σ = {sigma:.6f}\n")
        f.write(f"  半高宽 FWHM = {fwhm:.6f}\n")
        f.write(f"  在区域1宽度 (0.2) 中的比例: {fwhm/0.2:.1%}\n")
    
    print(f"\n参数信息已保存到 wave_parameters.txt")