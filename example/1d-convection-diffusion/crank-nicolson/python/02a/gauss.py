import numpy as np
import matplotlib.pyplot as plt
import configparser
import re
from scipy import sparse
from scipy.sparse.linalg import spsolve

def parse_value(value_str):
    """解析配置值，移除可能的注释"""
    if '#' in value_str:
        value_str = value_str.split('#')[0]
    return value_str.strip()

def read_config(filename='config_advanced.txt'):
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
    
    # 读取全局参数
    global_params = {
        'xmin': float(parse_value(config['Global']['xmin'])),
        'xmax': float(parse_value(config['Global']['xmax'])),
        'total_time': float(parse_value(config['Global']['total_time'])),
        'Nx': int(parse_value(config['Global']['Nx'])),
        'Nt': int(parse_value(config['Global']['Nt'])),
        't0': float(parse_value(config['Global']['t0'])),
        'initial_condition': parse_value(config['Global']['initial_condition'])
    }
    
    # 读取初始条件参数
    init_params = {
        'beta': float(parse_value(config['InitialCondition']['beta'])),
        'z': float(parse_value(config['InitialCondition']['z'])),
        'delta': float(parse_value(config['InitialCondition']['delta'])),
        'alpha': float(parse_value(config['InitialCondition']['alpha'])),
        'a_ellipse': float(parse_value(config['InitialCondition']['a']))
    }
    
    # 读取物理参数
    physics_params = {
        'a': float(parse_value(config['Physics']['a'])),
        'nu': float(parse_value(config['Physics']['nu']))
    }
    
    # 读取输出设置
    output_params = {
        'plot_title': parse_value(config['Output']['plot_title']),
        'save_figure': parse_value(config['Output']['save_figure']).lower() in ['true', 'yes', '1', 't'],
        'figure_name': parse_value(config['Output']['figure_name'])
    }
    
    # 读取验证参数
    verification_params = {
        'use_reference_solution': parse_value(config.get('Verification', 'use_reference_solution', 
                                                        fallback='True')).lower() in ['true', 'yes', '1', 't'],
        'ref_Nx_multiplier': int(parse_value(config.get('Verification', 'ref_Nx_multiplier', fallback='4'))),
        'ref_Nt_multiplier': int(parse_value(config.get('Verification', 'ref_Nt_multiplier', fallback='4'))),
        'use_high_order': parse_value(config.get('Verification', 'use_high_order', 
                                                 fallback='False')).lower() in ['true', 'yes', '1', 't']
    }
    
    return global_params, init_params, physics_params, output_params, verification_params

def combined_wave_initial(x, params):
    """组合波形初始条件"""
    beta = params['beta']
    z = params['z']
    delta = params['delta']
    alpha = params['alpha']
    a_ellipse = params['a_ellipse']
    
    u0 = np.zeros_like(x)
    
    def G(x, beta, z):
        return np.exp(-beta * (x - z)**2)
    
    def F(x, alpha, a):
        return np.sqrt(np.maximum(1 - alpha**2 * (x - a)**2, 0))
    
    # 区域1: 高斯波包 [-0.8, -0.6]
    mask1 = (-0.8 <= x) & (x <= -0.6)
    if np.any(mask1):
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
        u0[mask4] = (1/6) * (F(x[mask4], alpha, a_ellipse-delta) + 
                            F(x[mask4], alpha, a_ellipse+delta) + 
                            4*F(x[mask4], alpha, a_ellipse))
    
    return u0

def solve_convection_diffusion_cn(x, a, nu, total_time, Nt, u0):
    """
    使用Crank-Nicolson格式求解（二阶精度，无条件稳定，可能振荡）
    """
    Nx = len(x) - 1
    dx = x[1] - x[0]
    dt = total_time / Nt
    
    u = u0.copy()
    
    for n in range(Nt):
        N = Nx + 1
        main_diag = np.ones(N)
        lower_diag = np.zeros(N)
        upper_diag = np.zeros(N)
        b = np.zeros(N)
        
        # 内部点
        for i in range(1, Nx):
            # 中心差分（可能振荡）
            conv_coef = a / (2*dx)
            diff_coef = nu / (dx**2)
            
            main_diag[i] = 1.0 + dt * diff_coef
            lower_diag[i] = -0.5 * dt * (conv_coef + diff_coef)
            upper_diag[i] = 0.5 * dt * (conv_coef - diff_coef)
            
            b[i] = (1.0 - dt * diff_coef) * u[i] + \
                   0.5 * dt * (conv_coef + diff_coef) * u[i-1] + \
                   0.5 * dt * (-conv_coef + diff_coef) * u[i+1]
        
        # 边界条件：零梯度
        main_diag[0] = 1.0
        upper_diag[0] = -1.0
        b[0] = 0.0
        
        main_diag[Nx] = 1.0
        lower_diag[Nx] = -1.0
        b[Nx] = 0.0
        
        A = sparse.diags([lower_diag[1:], main_diag, upper_diag[:-1]], 
                        [-1, 0, 1], format='csr')
        u = spsolve(A, b)
    
    return u

def solve_convection_diffusion_upwind(x, a, nu, total_time, Nt, u0):
    """
    使用迎风格式求解（一阶精度，无振荡，耗散较大）
    更稳定，适合作为参考解
    """
    Nx = len(x) - 1
    dx = x[1] - x[0]
    dt = total_time / Nt
    
    u = u0.copy()
    
    # 迎风格式参数
    if a >= 0:
        # 向后差分
        conv_scheme = 'backward'
    else:
        # 向前差分
        conv_scheme = 'forward'
    
    for n in range(Nt):
        u_old = u.copy()
        
        for i in range(1, Nx):
            # 对流项：迎风格式
            if conv_scheme == 'backward':
                conv = a * (u_old[i] - u_old[i-1]) / dx
            else:
                conv = a * (u_old[i+1] - u_old[i]) / dx
            
            # 扩散项：中心差分
            diff = nu * (u_old[i+1] - 2*u_old[i] + u_old[i-1]) / (dx**2)
            
            u[i] = u_old[i] - dt * conv + dt * diff
        
        # 边界条件：零梯度
        u[0] = u[1]
        u[Nx] = u[Nx-1]
    
    return u

def solve_convection_diffusion_beam_warming(x, a, nu, total_time, Nt, u0):
    """
    使用Beam-Warming格式求解（二阶精度，迎风格式，更稳定）
    比中心差分更少振荡
    """
    Nx = len(x) - 1
    dx = x[1] - x[0]
    dt = total_time / Nt
    
    u = u0.copy()
    
    for n in range(Nt):
        u_old = u.copy()
        
        for i in range(2, Nx-1):
            # Beam-Warming格式（二阶迎风）
            if a >= 0:
                # 向后差分，二阶精度
                conv = a * (3*u_old[i] - 4*u_old[i-1] + u_old[i-2]) / (2*dx)
            else:
                # 向前差分，二阶精度
                conv = a * (-3*u_old[i] + 4*u_old[i+1] - u_old[i+2]) / (2*dx)
            
            # 扩散项：中心差分
            diff = nu * (u_old[i+1] - 2*u_old[i] + u_old[i-1]) / (dx**2)
            
            u[i] = u_old[i] - dt * conv + dt * diff
        
        # 边界点：使用一阶格式
        for i in [1, Nx-1]:
            if a >= 0:
                conv = a * (u_old[i] - u_old[i-1]) / dx
            else:
                conv = a * (u_old[i+1] - u_old[i]) / dx
            
            diff = nu * (u_old[i+1] - 2*u_old[i] + u_old[i-1]) / (dx**2)
            u[i] = u_old[i] - dt * conv + dt * diff
        
        # 边界条件：零梯度
        u[0] = u[1]
        u[Nx] = u[Nx-1]
    
    return u

def compute_reference_solution(x_ref, a, nu, total_time, u0_func, params, 
                              use_high_order=False, ref_Nt_multiplier=4):
    """
    计算参考解（使用密网格和稳定格式）
    """
    print("计算参考解...")
    
    # 参考解网格更密
    Nx_ref = len(x_ref) - 1
    Nt_ref = ref_Nt_multiplier * Nx_ref  # 时间步也更密
    
    # 计算参考解的初始条件
    u0_ref = u0_func(x_ref, params)
    
    # 选择格式：如果要高阶且稳定，用Beam-Warming；否则用迎风
    if use_high_order:
        print(f"  使用Beam-Warming格式（二阶迎风）")
        print(f"  参考网格: Nx={Nx_ref}, Nt={Nt_ref}")
        u_ref = solve_convection_diffusion_beam_warming(
            x_ref, a, nu, total_time, Nt_ref, u0_ref
        )
    else:
        print(f"  使用迎风格式（一阶，无振荡）")
        print(f"  参考网格: Nx={Nx_ref}, Nt={Nt_ref}")
        u_ref = solve_convection_diffusion_upwind(
            x_ref, a, nu, total_time, Nt_ref, u0_ref
        )
    
    return u_ref

def interpolate_to_coarse_grid(x_fine, u_fine, x_coarse):
    """
    将细网格的解插值到粗网格上
    """
    return np.interp(x_coarse, x_fine, u_fine)

def calculate_convergence_metrics(u_num, u_ref, x, label=""):
    """
    计算收敛性指标
    """
    errors = np.abs(u_num - u_ref)
    
    metrics = {
        'max_error': np.max(errors),
        'l1_error': np.trapz(errors, x),
        'l2_error': np.sqrt(np.trapz(errors**2, x)),
        'relative_max_error': np.max(errors) / (np.max(np.abs(u_ref)) + 1e-12),
        'mass_error': abs(np.trapz(u_num, x) - np.trapz(u_ref, x)) / (abs(np.trapz(u_ref, x)) + 1e-12)
    }
    
    if label:
        print(f"\n{label}误差分析:")
        print(f"  最大绝对误差: {metrics['max_error']:.2e}")
        print(f"  L1误差: {metrics['l1_error']:.2e}")
        print(f"  L2误差: {metrics['l2_error']:.2e}")
        print(f"  相对最大误差: {metrics['relative_max_error']:.2e}")
        print(f"  质量相对误差: {metrics['mass_error']:.2e}")
    
    return metrics

def plot_results(x, u_initial, u_num, u_ref, params, output_params, 
                 metrics_cn=None, metrics_upwind=None):
    """
    绘制结果对比图
    """
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    
    a = params['a']
    nu = params['nu']
    total_time = params['total_time']
    
    # 1. 初始条件和最终数值解
    ax = axes[0, 0]
    ax.plot(x, u_initial, 'g-', linewidth=2, alpha=0.7, label='Initial (t=0)')
    ax.plot(x, u_num, 'b-', linewidth=2, label=f'Numerical (t={total_time})')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title(f'Initial vs Numerical\nConvection a={a}, Diffusion ν={nu:.4f}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. 数值解与参考解对比
    ax = axes[0, 1]
    ax.plot(x, u_num, 'b-', linewidth=2, label='Numerical (coarse grid)')
    ax.plot(x, u_ref, 'r--', linewidth=2, alpha=0.7, label='Reference (fine grid)')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title('Numerical vs Reference Solution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 添加误差信息
    if metrics_cn:
        error_text = f'Max error: {metrics_cn["max_error"]:.2e}\nL2 error: {metrics_cn["l2_error"]:.2e}'
        ax.text(0.02, 0.98, error_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 3. 误差分布
    ax = axes[0, 2]
    error = u_num - u_ref
    ax.plot(x, error, 'r-', linewidth=1.5)
    ax.set_xlabel('x')
    ax.set_ylabel('Error (Numerical - Reference)')
    
    if metrics_cn:
        ax.set_title(f'Error Distribution (max={metrics_cn["max_error"]:.2e})')
    else:
        ax.set_title('Error Distribution')
    
    ax.grid(True, alpha=0.3)
    
    # 4. 局部放大（峰值区域）
    ax = axes[1, 0]
    # 找到最大值位置
    idx_max = np.argmax(u_ref)
    x_center = x[idx_max]
    
    # 放大区域
    zoom_width = 0.3
    zoom_mask = (x >= x_center - zoom_width) & (x <= x_center + zoom_width)
    
    if np.any(zoom_mask):
        ax.plot(x[zoom_mask], u_num[zoom_mask], 'b-o', markersize=4, label='Numerical')
        ax.plot(x[zoom_mask], u_ref[zoom_mask], 'r--s', markersize=4, label='Reference')
        ax.set_xlabel('x')
        ax.set_ylabel('u(x,t)')
        ax.set_title(f'Zoom: Peak Region (x≈{x_center:.2f})')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # 5. 对数坐标下的误差
    ax = axes[1, 1]
    abs_error = np.abs(error)
    # 避免log(0)
    abs_error[abs_error < 1e-12] = 1e-12
    
    ax.semilogy(x, abs_error, 'r-', linewidth=1.5)
    ax.set_xlabel('x')
    ax.set_ylabel('|Error| (log scale)')
    ax.set_title('Absolute Error (Log Scale)')
    ax.grid(True, alpha=0.3, which='both')
    
    # 6. 格式对比（如果可用）
    ax = axes[1, 2]
    if metrics_upwind:
        # 绘制两种格式的误差对比
        formats = ['Crank-Nicolson', 'Upwind']
        max_errors = [metrics_cn['max_error'], metrics_upwind['max_error']]
        l2_errors = [metrics_cn['l2_error'], metrics_upwind['l2_error']]
        
        x_pos = np.arange(len(formats))
        width = 0.35
        
        ax.bar(x_pos - width/2, max_errors, width, label='Max Error', alpha=0.8)
        ax.bar(x_pos + width/2, l2_errors, width, label='L2 Error', alpha=0.8)
        
        ax.set_xlabel('Numerical Scheme')
        ax.set_ylabel('Error')
        ax.set_title('Error Comparison: Different Schemes')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(formats)
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_yscale('log')
    else:
        # 显示参考解的质量守恒
        mass_initial = np.trapz(u_initial, x)
        mass_final = np.trapz(u_ref, x)
        mass_change = abs(mass_final - mass_initial) / abs(mass_initial)
        
        ax.text(0.5, 0.7, f'Mass Conservation:\nInitial: {mass_initial:.6f}\nFinal: {mass_final:.6f}\nRelative Change: {mass_change:.2e}',
               ha='center', va='center', fontsize=12,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax.text(0.5, 0.3, 'Note: For pure convection,\nmass should be conserved.\nFor diffusion, mass is conserved.',
               ha='center', va='center', fontsize=10)
        ax.set_title('Mass Conservation Check')
        ax.axis('off')
    
    plt.suptitle(output_params['plot_title'] + '\n(Reference solution from fine grid)', 
                 fontsize=16, y=1.02)
    plt.tight_layout()
    
    if output_params['save_figure']:
        plt.savefig(output_params['figure_name'], dpi=300, bbox_inches='tight')
        print(f"图像已保存为: {output_params['figure_name']}")
    
    plt.show()

def main():
    """主函数"""
    print("="*60)
    print("对流扩散方程数值解验证（使用参考解）")
    print("="*60)
    
    try:
        # 读取配置文件
        print("正在读取配置文件 config_advanced.txt...")
        global_params, init_params, physics_params, output_params, verification_params = read_config('config_advanced.txt')
        
        # 合并参数
        all_params = {**global_params, **init_params, **physics_params}
        
        print(f"计算域: x∈[{global_params['xmin']}, {global_params['xmax']}]")
        print(f"时间: T={global_params['total_time']}")
        print(f"基本网格: Nx={global_params['Nx']}, Nt={global_params['Nt']}")
        print(f"对流速度: a={physics_params['a']}, 扩散系数: ν={physics_params['nu']}")
        print(f"初始条件类型: {global_params['initial_condition']}")
        
        use_ref = verification_params['use_reference_solution']
        ref_Nx_multiplier = verification_params['ref_Nx_multiplier']
        use_high_order = verification_params['use_high_order']
        
        # 创建基本网格
        x = np.linspace(global_params['xmin'], global_params['xmax'], global_params['Nx']+1)
        
        # 计算初始条件
        print("\n正在计算初始条件...")
        u_initial = combined_wave_initial(x, all_params)
        
        # 计算Crank-Nicolson格式数值解
        print("\n计算Crank-Nicolson格式数值解...")
        u_cn = solve_convection_diffusion_cn(
            x, physics_params['a'], physics_params['nu'],
            global_params['total_time'], global_params['Nt'],
            u_initial
        )
        
        # 计算参考解（如果需要）
        u_ref = None
        metrics_cn = None
        
        if use_ref:
            # 创建密网格用于参考解
            Nx_ref = ref_Nx_multiplier * global_params['Nx']
            x_ref = np.linspace(global_params['xmin'], global_params['xmax'], Nx_ref+1)
            
            # 计算参考解
            u_ref = compute_reference_solution(
                x_ref, physics_params['a'], physics_params['nu'],
                global_params['total_time'], combined_wave_initial, all_params,
                use_high_order, verification_params['ref_Nt_multiplier']
            )
            
            # 将参考解插值到基本网格
            u_ref_coarse = interpolate_to_coarse_grid(x_ref, u_ref, x)
            
            # 计算误差指标
            metrics_cn = calculate_convergence_metrics(
                u_cn, u_ref_coarse, x, "Crank-Nicolson格式"
            )
            
            # 可选：计算迎风格式的解进行对比
            print("\n计算迎风格式数值解（对比）...")
            u_upwind = solve_convection_diffusion_upwind(
                x, physics_params['a'], physics_params['nu'],
                global_params['total_time'], global_params['Nt'],
                u_initial
            )
            
            metrics_upwind = calculate_convergence_metrics(
                u_upwind, u_ref_coarse, x, "迎风格式"
            )
        else:
            u_ref_coarse = None
            metrics_upwind = None
            print("\n跳过参考解计算（根据配置设置）")
        
        # 绘制结果
        print("\n正在生成图像...")
        plot_results(x, u_initial, u_cn, u_ref_coarse, all_params, 
                    output_params, metrics_cn, metrics_upwind)
        
        # 网格收敛性测试
        print("\n" + "="*60)
        print("网格收敛性测试")
        print("="*60)
        
        if use_ref and u_ref is not None:
            # 测试不同网格的收敛性
            Nx_list = [global_params['Nx']//4, global_params['Nx']//2, 
                      global_params['Nx'], global_params['Nx']*2]
            
            errors = []
            
            for Nx_test in Nx_list:
                if Nx_test < 20:  # 网格太粗跳过
                    continue
                    
                print(f"\n测试网格: Nx={Nx_test}")
                x_test = np.linspace(global_params['xmin'], global_params['xmax'], Nx_test+1)
                
                # 计算初始条件
                u0_test = combined_wave_initial(x_test, all_params)
                
                # 计算数值解
                Nt_test = max(global_params['Nt'], Nx_test)
                u_test = solve_convection_diffusion_cn(
                    x_test, physics_params['a'], physics_params['nu'],
                    global_params['total_time'], Nt_test, u0_test
                )
                
                # 插值到测试网格并计算误差
                u_ref_test = interpolate_to_coarse_grid(x_ref, u_ref, x_test)
                error_l2 = np.sqrt(np.trapz((u_test - u_ref_test)**2, x_test))
                errors.append((Nx_test, error_l2))
                
                print(f"  L2误差: {error_l2:.2e}")
            
            # 计算收敛阶
            if len(errors) >= 2:
                print("\n收敛阶估计:")
                for i in range(1, len(errors)):
                    Nx1, err1 = errors[i-1]
                    Nx2, err2 = errors[i]
                    order = np.log(err1/err2) / np.log(Nx2/Nx1)
                    print(f"  {Nx1} → {Nx2}: 收敛阶 ≈ {order:.2f}")
        
        print("\n" + "="*60)
        print("验证完成")
        print("="*60)
        
    except FileNotFoundError:
        print("错误: 找不到配置文件 config_advanced.txt")
        print("请创建配置文件或检查文件路径")
    except Exception as e:
        print(f"运行时错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()