import numpy as np
import matplotlib.pyplot as plt
import configparser

def parse_value(value_str):
    """解析配置值，移除可能的注释"""
    if '#' in value_str:
        value_str = value_str.split('#')[0]
    return value_str.strip()

def read_config(filename='config_diffusion_effect.txt'):
    """从配置文件读取参数"""
    config = configparser.ConfigParser()
    
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 移除行内注释
    lines = content.split('\n')
    cleaned_lines = []
    for line in lines:
        if '#' in line:
            line = line.split('#')[0]
        cleaned_lines.append(line)
    
    cleaned_content = '\n'.join(cleaned_lines)
    config.read_string(cleaned_content)
    
    # 读取参数
    params = {
        'xmin': float(parse_value(config['Domain']['xmin'])),
        'xmax': float(parse_value(config['Domain']['xmax'])),
        'total_time': float(parse_value(config['Time']['total_time'])),
        'Nx': int(parse_value(config['Mesh']['Nx'])),
        'Nt': int(parse_value(config['Mesh']['Nt'])),
        'a': float(parse_value(config['Physics']['a'])),  # 对流速度
        'nu_list': parse_value(config['Physics']['nu_list']),  # 扩散系数列表
        'case': parse_value(config['Case']['type']),  # 测试案例类型
        'include_zero_diffusion': parse_value(config.get('Options', 'include_zero_diffusion', 
                                                        fallback='True')).lower() in ['true', 'yes', '1', 't']
    }
    
    # 解析扩散系数列表
    nu_values = []
    for nu_str in params['nu_list'].split(','):
        nu_str = nu_str.strip()
        if nu_str:
            nu_values.append(float(nu_str))
    params['nu_values'] = nu_values
    
    return params

def initial_condition(x, case_type):
    """
    初始条件函数
    """
    if case_type == 'gaussian':
        # 高斯分布
        t0 = 0.1
        nu_small = 0.01  # 为了计算高斯宽度
        return np.exp(-x**2 / (4*nu_small*t0)) / np.sqrt(4*np.pi*nu_small*t0)
    
    elif case_type == 'combined_wave':
        # 组合波形
        u0 = np.zeros_like(x)
        
        # 区域1: 高斯波包 [-0.8, -0.6]
        mask1 = (-0.8 <= x) & (x <= -0.6)
        if np.any(mask1):
            beta = 100.0
            z = -0.7
            delta = 0.005
            def G(x, beta, z):
                return np.exp(-beta * (x - z)**2)
            u0[mask1] = (1/6) * (G(x[mask1], beta, z-delta) + 
                                G(x[mask1], beta, z+delta) + 
                                4*G(x[mask1], beta, z))
        
        # 区域2: 方形波 [-0.4, -0.2]
        mask2 = (-0.4 <= x) & (x <= -0.2)
        u0[mask2] = 1.0
        
        # 区域3: 三角形波 [0, 0.2]
        mask3 = (0 <= x) & (x <= 0.2)
        u0[mask3] = 1 - np.abs(10 * (x[mask3] - 0.1))
        
        # 区域4: 半椭圆波 [0.4, 0.6]
        mask4 = (0.4 <= x) & (x <= 0.6)
        if np.any(mask4):
            alpha = 10.0
            a_ellipse = 0.5
            delta = 0.005
            def F(x, alpha, a):
                return np.sqrt(np.maximum(1 - alpha**2 * (x - a)**2, 0))
            u0[mask4] = (1/6) * (F(x[mask4], alpha, a_ellipse-delta) + 
                                F(x[mask4], alpha, a_ellipse+delta) + 
                                4*F(x[mask4], alpha, a_ellipse))
        return u0
    
    elif case_type == 'square_wave':
        # 简单方波
        u0 = np.zeros_like(x)
        mask = (-0.5 <= x) & (x <= 0.5)
        u0[mask] = 1.0
        return u0
    
    else:
        raise ValueError(f"未知的案例类型: {case_type}")

def exact_solution(x, t, a, case_type):
    """
    纯对流方程的精确解: u(x,t) = u₀(x-at)
    注意：这是"理论"精确解，假设没有扩散
    """
    # 计算在x-at位置上的初始条件
    x_shifted = x - a * t
    return initial_condition(x_shifted, case_type)

def solve_convection_diffusion_cn(x, a, nu, total_time, Nt, u0):
    """
    使用Crank-Nicolson格式求解对流扩散方程
    方程: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²
    """
    Nx = len(x) - 1
    dx = x[1] - x[0]
    dt = total_time / Nt
    
    u = u0.copy()
    
    # CFL数和扩散数
    CFL = a * dt / dx
    diffusion_number = nu * dt / (dx**2)
    
    if CFL > 1.0:
        print(f"  警告: CFL数={CFL:.3f} > 1，可能不稳定")
    if diffusion_number > 0.5:
        print(f"  警告: 扩散数={diffusion_number:.3f} > 0.5，可能不稳定")
    
    # Crank-Nicolson格式时间推进
    for n in range(Nt):
        # 构造三对角矩阵
        N = Nx + 1
        main_diag = np.ones(N)
        lower_diag = np.zeros(N)
        upper_diag = np.zeros(N)
        b = np.zeros(N)
        
        # 系数
        conv_coef = a / (2*dx)  # 对流项系数（中心差分）
        diff_coef = nu / (dx**2)  # 扩散项系数
        
        # 内部点
        for i in range(1, Nx):
            # 隐式矩阵系数 (I + 0.5*dt*L)
            main_diag[i] = 1.0 + dt * diff_coef
            lower_diag[i] = -0.5 * dt * (conv_coef + diff_coef)
            upper_diag[i] = 0.5 * dt * (conv_coef - diff_coef)
            
            # 右端向量 ((I - 0.5*dt*L)u^n)
            b[i] = (1.0 - dt * diff_coef) * u[i] + \
                   0.5 * dt * (conv_coef + diff_coef) * u[i-1] + \
                   0.5 * dt * (-conv_coef + diff_coef) * u[i+1]
        
        # 边界条件：零梯度
        main_diag[0] = 1.0
        upper_diag[0] = -1.0  # u₀ = u₁
        b[0] = 0.0
        
        main_diag[Nx] = 1.0
        lower_diag[Nx] = -1.0  # u_N = u_{N-1}
        b[Nx] = 0.0
        
        # Thomas算法求解三对角系统
        # 前向消元
        c_prime = np.zeros(N-1)
        d_prime = np.zeros(N)
        
        c_prime[0] = upper_diag[0] / main_diag[0]
        d_prime[0] = b[0] / main_diag[0]
        
        for i in range(1, N):
            if i < N-1:
                denom = main_diag[i] - lower_diag[i] * c_prime[i-1]
                c_prime[i] = upper_diag[i] / denom
            denom = main_diag[i] - lower_diag[i] * c_prime[i-1]
            d_prime[i] = (b[i] - lower_diag[i] * d_prime[i-1]) / denom
        
        # 回代
        u_new = np.zeros(N)
        u_new[-1] = d_prime[-1]
        for i in range(N-2, -1, -1):
            u_new[i] = d_prime[i] - c_prime[i] * u_new[i+1]
        
        u = u_new
    
    return u

def calculate_errors(u_num, u_exact, x):
    """计算各种误差指标"""
    error = u_num - u_exact
    
    metrics = {
        'max_error': np.max(np.abs(error)),
        'l2_error': np.sqrt(np.trapz(error**2, x)),
        'l1_error': np.trapz(np.abs(error), x),
        'relative_max_error': np.max(np.abs(error)) / (np.max(np.abs(u_exact)) + 1e-12),
        'mass_error': abs(np.trapz(u_num, x) - np.trapz(u_exact, x)) / (abs(np.trapz(u_exact, x)) + 1e-12),
        'peak_height': np.max(u_num) if len(u_num) > 0 else 0,
        'peak_position': x[np.argmax(u_num)] if len(u_num) > 0 else 0
    }
    
    return metrics

def plot_results(x, u_initial, u_exact, solutions, errors_list, params):
    """
    绘制对比图
    """
    a = params['a']
    total_time = params['total_time']
    case_type = params['case']
    nu_values = params['nu_values']
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    
    # 1. 所有数值解与精确解对比
    ax = axes[0, 0]
    ax.plot(x, u_initial, 'k:', linewidth=1.5, alpha=0.5, label='Initial (t=0)')
    ax.plot(x, u_exact, 'r--', linewidth=2.5, label='Exact (pure convection)')
    
    colors = ['b', 'g', 'm', 'c', 'y', 'orange']
    for i, (nu, u_num) in enumerate(solutions):
        color = colors[i % len(colors)]
        label = f'ν={nu:.2e}'
        if nu == 0:
            label = 'ν=0 (pure convection)'
        ax.plot(x, u_num, color=color, linewidth=1.5, label=label, alpha=0.8)
    
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title(f'Effect of Diffusion Coefficient\nConvection a={a}, Time T={total_time}')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # 2. 误差随扩散系数的变化
    ax = axes[0, 1]
    nu_list = [nu for nu, _ in solutions]
    max_errors = [err['max_error'] for err in errors_list]
    l2_errors = [err['l2_error'] for err in errors_list]
    
    ax.loglog(nu_list, max_errors, 'bo-', linewidth=2, markersize=8, label='Max Error')
    ax.loglog(nu_list, l2_errors, 'rs--', linewidth=2, markersize=8, label='L2 Error')
    ax.set_xlabel('Diffusion Coefficient ν (log scale)')
    ax.set_ylabel('Error (log scale)')
    ax.set_title('Error vs Diffusion Coefficient')
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')
    
    # 3. 峰值高度随扩散系数的变化
    ax = axes[0, 2]
    peak_heights = [err['peak_height'] for err in errors_list]
    exact_peak = np.max(u_exact)
    
    ax.semilogx(nu_list, peak_heights, 'go-', linewidth=2, markersize=8, label='Numerical peak')
    ax.axhline(y=exact_peak, color='r', linestyle='--', linewidth=2, label='Exact peak')
    ax.set_xlabel('Diffusion Coefficient ν (log scale)')
    ax.set_ylabel('Peak Height')
    ax.set_title('Peak Height vs Diffusion Coefficient')
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')
    
    # 4. 质量误差随扩散系数的变化
    ax = axes[1, 0]
    mass_errors = [err['mass_error'] for err in errors_list]
    
    ax.loglog(nu_list, mass_errors, 'mo-', linewidth=2, markersize=8)
    ax.set_xlabel('Diffusion Coefficient ν (log scale)')
    ax.set_ylabel('Mass Relative Error (log scale)')
    ax.set_title('Mass Conservation Error')
    ax.grid(True, alpha=0.3, which='both')
    
    # 5. 局部放大：峰值区域
    ax = axes[1, 1]
    # 找到精确解的峰值位置
    idx_exact_peak = np.argmax(u_exact)
    x_center = x[idx_exact_peak]
    zoom_width = 0.3
    
    zoom_mask = (x >= x_center - zoom_width) & (x <= x_center + zoom_width)
    
    if np.any(zoom_mask):
        ax.plot(x[zoom_mask], u_exact[zoom_mask], 'r--', linewidth=2.5, label='Exact')
        
        for i, (nu, u_num) in enumerate(solutions[:3]):  # 只显示前3个
            color = colors[i % len(colors)]
            label = f'ν={nu:.2e}' if nu > 0 else 'ν=0'
            ax.plot(x[zoom_mask], u_num[zoom_mask], color=color, linewidth=1.5, label=label, alpha=0.8)
        
        ax.set_xlabel('x')
        ax.set_ylabel('u(x,t)')
        ax.set_title(f'Zoom: Peak Region (x≈{x_center:.2f})')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
    
    # 6. 数值解与精确解的绝对误差
    ax = axes[1, 2]
    # 选择一个小扩散系数和一个大扩散系数的误差
    if len(solutions) >= 2:
        # 最小的非零扩散系数
        small_nu_idx = 1 if nu_values[0] == 0 else 0
        # 最大的扩散系数
        large_nu_idx = len(solutions) - 1
        
        if small_nu_idx < len(solutions):
            nu_small, u_small = solutions[small_nu_idx]
            error_small = np.abs(u_small - u_exact)
            ax.plot(x, error_small, 'g-', linewidth=1.5, label=f'ν={nu_small:.2e}', alpha=0.8)
        
        if large_nu_idx != small_nu_idx and large_nu_idx < len(solutions):
            nu_large, u_large = solutions[large_nu_idx]
            error_large = np.abs(u_large - u_exact)
            ax.plot(x, error_large, 'b-', linewidth=1.5, label=f'ν={nu_large:.2e}', alpha=0.8)
    
    ax.set_xlabel('x')
    ax.set_ylabel('Absolute Error')
    ax.set_title('Error Distribution for Different ν')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.suptitle(f'Diffusion Effect Study: {case_type} Initial Condition', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()
    
    # 打印误差表格
    print("\n" + "="*80)
    print("误差分析汇总表:")
    print("="*80)
    print(f"{'ν':<12} {'Max Error':<15} {'L2 Error':<15} {'L1 Error':<15} {'Mass Error':<15} {'Peak Height':<15}")
    print("-"*80)
    
    for i, (nu, metrics) in enumerate(zip([nu for nu, _ in solutions], errors_list)):
        print(f"{nu:<12.2e} {metrics['max_error']:<15.2e} {metrics['l2_error']:<15.2e} "
              f"{metrics['l1_error']:<15.2e} {metrics['mass_error']:<15.2e} {metrics['peak_height']:<15.4f}")

def main():
    """
    主函数：研究扩散系数对计算结果的影响
    """
    print("="*70)
    print("扩散系数影响研究：精确解为纯对流，数值解包含扩散项")
    print("方程: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²")
    print("精确解（理论）: u(x,t) = u₀(x-at) (假设ν=0)")
    print("数值解: 包含扩散项ν的影响")
    print("="*70)
    
    try:
        # 读取配置文件
        print("正在读取配置文件 config_diffusion_effect.txt...")
        params = read_config('config_diffusion_effect.txt')
        
        print(f"计算域: x∈[{params['xmin']}, {params['xmax']}]")
        print(f"时间: T={params['total_time']}")
        print(f"网格: Nx={params['Nx']}, Nt={params['Nt']}")
        print(f"对流速度: a={params['a']}")
        print(f"测试案例: {params['case']}")
        print(f"扩散系数列表: {params['nu_values']}")
        
        # 创建空间网格
        x = np.linspace(params['xmin'], params['xmax'], params['Nx']+1)
        
        # 计算初始条件
        print("\n正在计算初始条件...")
        u_initial = initial_condition(x, params['case'])
        
        # 计算精确解（纯对流，ν=0）
        print("正在计算精确解（纯对流）...")
        u_exact = exact_solution(x, params['total_time'], params['a'], params['case'])
        
        # 存储所有解
        solutions = []
        errors_list = []
        
        # 为每个扩散系数计算数值解
        for nu in params['nu_values']:
            print(f"\n计算扩散系数 ν={nu:.2e} 的数值解...")
            
            if nu == 0:
                print("  注意: ν=0，这是纯对流情况，数值解应与精确解一致")
            
            # 计算数值解
            u_num = solve_convection_diffusion_cn(
                x, params['a'], nu, params['total_time'], params['Nt'], u_initial
            )
            
            # 计算误差
            metrics = calculate_errors(u_num, u_exact, x)
            
            solutions.append((nu, u_num))
            errors_list.append(metrics)
            
            print(f"  最大误差: {metrics['max_error']:.2e}")
            print(f"  L2误差: {metrics['l2_error']:.2e}")
            print(f"  峰值高度: {metrics['peak_height']:.4f} (精确解: {np.max(u_exact):.4f})")
            print(f"  峰值位置: x={metrics['peak_position']:.3f} (精确解: x={x[np.argmax(u_exact)]:.3f})")
        
        # 绘制结果
        print("\n正在生成图像...")
        plot_results(x, u_initial, u_exact, solutions, errors_list, params)
        
        # 分析扩散效应
        print("\n" + "="*70)
        print("扩散效应分析:")
        print("="*70)
        
        # 找到最小和最大的扩散系数（非零）
        non_zero_nu = [nu for nu in params['nu_values'] if nu > 0]
        if non_zero_nu:
            min_nu = min(non_zero_nu)
            max_nu = max(non_zero_nu)
            
            # 找到对应的索引
            min_idx = params['nu_values'].index(min_nu)
            max_idx = params['nu_values'].index(max_nu)
            
            min_metrics = errors_list[min_idx]
            max_metrics = errors_list[max_idx]
            
            print(f"最小扩散系数 ν={min_nu:.2e}:")
            print(f"  最大误差: {min_metrics['max_error']:.2e}")
            print(f"  峰值高度衰减: {1 - min_metrics['peak_height']/np.max(u_exact):.1%}")
            
            print(f"\n最大扩散系数 ν={max_nu:.2e}:")
            print(f"  最大误差: {max_metrics['max_error']:.2e}")
            print(f"  峰值高度衰减: {1 - max_metrics['peak_height']/np.max(u_exact):.1%}")
            
            # 计算误差增长比例
            error_growth = max_metrics['max_error'] / min_metrics['max_error']
            print(f"\n误差增长比例: {error_growth:.2f}倍")
            
            # 扩散特征长度
            diffusion_length = np.sqrt(2 * max_nu * params['total_time'])
            print(f"最大扩散的特征长度: √(2νT) = {diffusion_length:.4f}")
        
        print("\n" + "="*70)
        print("研究完成")
        print("="*70)
        
    except FileNotFoundError:
        print("错误: 找不到配置文件 config_diffusion_effect.txt")
        print("正在创建默认配置文件...")
        
        # 创建默认配置文件
        default_config = """# 扩散效应研究参数配置

[Domain]
xmin = -1.5
xmax = 1.5

[Time]
total_time = 1.0

[Mesh]
Nx = 300
Nt = 2000

[Physics]
a = 1.0
nu_list = 0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2

[Case]
type = combined_wave  # 可选: gaussian, combined_wave, square_wave

[Options]
include_zero_diffusion = True
"""
        
        with open('config_diffusion_effect.txt', 'w', encoding='utf-8') as f:
            f.write(default_config)
        
        print("已创建默认配置文件 config_diffusion_effect.txt")
        print("请重新运行程序")
        
    except Exception as e:
        print(f"运行时错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()