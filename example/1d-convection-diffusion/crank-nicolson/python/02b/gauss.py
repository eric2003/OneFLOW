import numpy as np
import matplotlib.pyplot as plt
import configparser

def parse_value(value_str):
    """解析配置值，移除可能的注释"""
    if '#' in value_str:
        value_str = value_str.split('#')[0]
    return value_str.strip()

def read_config(filename='config_convection.txt'):
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
        'case': parse_value(config['Case']['type'])  # 测试案例类型
    }
    
    return params

def initial_condition(x, case_type):
    """
    初始条件函数
    """
    if case_type == 'gaussian':
        # 高斯分布
        t0 = 0.1
        return np.exp(-x**2 / (4*0.1*t0)) / np.sqrt(4*np.pi*0.1*t0)
    
    elif case_type == 'combined_wave':
        # 组合波形（纯对流情况下有解析解）
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
    """
    # 计算在x-at位置上的初始条件
    x_shifted = x - a * t
    return initial_condition(x_shifted, case_type)

def solve_pure_convection_cn(x, a, total_time, Nt, u0):
    """
    使用Crank-Nicolson格式求解纯对流方程
    方程: ∂u/∂t + a ∂u/∂x = 0
    """
    Nx = len(x) - 1
    dx = x[1] - x[0]
    dt = total_time / Nt
    
    u = u0.copy()
    
    # Crank-Nicolson格式时间推进
    for n in range(Nt):
        # 构造三对角矩阵
        N = Nx + 1
        main_diag = np.ones(N)
        lower_diag = np.zeros(N)
        upper_diag = np.zeros(N)
        b = np.zeros(N)
        
        # 对流项系数（中心差分）
        conv_coef = a / (2*dx)
        
        # 内部点
        for i in range(1, Nx):
            # 隐式矩阵系数
            main_diag[i] = 1.0
            lower_diag[i] = -0.5 * dt * conv_coef
            upper_diag[i] = 0.5 * dt * conv_coef
            
            # 右端向量
            b[i] = u[i] + 0.5 * dt * conv_coef * u[i-1] - 0.5 * dt * conv_coef * u[i+1]
        
        # 边界条件：零梯度（周期性边界更合适，但零梯度也可以）
        main_diag[0] = 1.0
        upper_diag[0] = -1.0  # u₀ = u₁
        b[0] = 0.0
        
        main_diag[Nx] = 1.0
        lower_diag[Nx] = -1.0  # u_N = u_{N-1}
        b[Nx] = 0.0
        
        # 构建并求解三对角系统（使用Thomas算法）
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

def solve_pure_convection_upwind(x, a, total_time, Nt, u0):
    """
    使用迎风格式求解纯对流方程
    """
    Nx = len(x) - 1
    dx = x[1] - x[0]
    dt = total_time / Nt
    
    u = u0.copy()
    
    # 时间步进
    for n in range(Nt):
        u_old = u.copy()
        
        if a >= 0:
            # 向后差分
            for i in range(1, Nx+1):
                u[i] = u_old[i] - a * dt / dx * (u_old[i] - u_old[i-1])
        else:
            # 向前差分
            for i in range(0, Nx):
                u[i] = u_old[i] - a * dt / dx * (u_old[i+1] - u_old[i])
        
        # 边界条件
        u[0] = u[1] if a >= 0 else u_old[0]
        u[Nx] = u[Nx-1] if a >= 0 else u_old[Nx]
    
    return u

def plot_comparison(x, u_initial, u_cn, u_upwind, u_exact, params, errors):
    """
    绘制对比图
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    a = params['a']
    total_time = params['total_time']
    case_type = params['case']
    
    # 1. 所有解对比
    ax = axes[0, 0]
    ax.plot(x, u_initial, 'g-', linewidth=2, alpha=0.5, label='Initial (t=0)')
    ax.plot(x, u_cn, 'b-', linewidth=2, label=f'Crank-Nicolson')
    ax.plot(x, u_upwind, 'm-', linewidth=2, alpha=0.8, label='Upwind')
    ax.plot(x, u_exact, 'r--', linewidth=2, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title(f'Pure Convection: a={a}, T={total_time}\nCase: {case_type}')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Crank-Nicolson与精确解对比
    ax = axes[0, 1]
    ax.plot(x, u_cn, 'b-', linewidth=2, label='Crank-Nicolson')
    ax.plot(x, u_exact, 'r--', linewidth=2, alpha=0.7, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title('Crank-Nicolson vs Exact')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 添加误差信息
    error_text = f'Max error: {errors["cn_max"]:.2e}\nL2 error: {errors["cn_l2"]:.2e}'
    ax.text(0.02, 0.98, error_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 3. 迎风格式与精确解对比
    ax = axes[1, 0]
    ax.plot(x, u_upwind, 'm-', linewidth=2, label='Upwind')
    ax.plot(x, u_exact, 'r--', linewidth=2, alpha=0.7, label='Exact')
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')
    ax.set_title('Upwind vs Exact')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    error_text = f'Max error: {errors["upwind_max"]:.2e}\nL2 error: {errors["upwind_l2"]:.2e}'
    ax.text(0.02, 0.98, error_text, transform=ax.transAxes,
           fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 4. 误差分布
    ax = axes[1, 1]
    error_cn = np.abs(u_cn - u_exact)
    error_upwind = np.abs(u_upwind - u_exact)
    
    ax.plot(x, error_cn, 'b-', linewidth=1.5, label='Crank-Nicolson error')
    ax.plot(x, error_upwind, 'm-', linewidth=1.5, label='Upwind error')
    ax.set_xlabel('x')
    ax.set_ylabel('Absolute Error')
    ax.set_title('Error Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # 打印总结
    print("\n" + "="*60)
    print("误差分析总结:")
    print("="*60)
    print(f"Crank-Nicolson格式:")
    print(f"  最大绝对误差: {errors['cn_max']:.2e}")
    print(f"  L2误差: {errors['cn_l2']:.2e}")
    print(f"  质量相对误差: {errors['cn_mass']:.2e}")
    
    print(f"\n迎风格式:")
    print(f"  最大绝对误差: {errors['upwind_max']:.2e}")
    print(f"  L2误差: {errors['upwind_l2']:.2e}")
    print(f"  质量相对误差: {errors['upwind_mass']:.2e}")

def main():
    """
    主函数：纯对流方程验证
    """
    print("="*60)
    print("纯对流方程数值解与精确解对比验证")
    print("方程: ∂u/∂t + a ∂u/∂x = 0")
    print("精确解: u(x,t) = u₀(x-at)")
    print("="*60)
    
    try:
        # 读取配置文件
        print("正在读取配置文件 config_convection.txt...")
        params = read_config('config_convection.txt')
        
        print(f"计算域: x∈[{params['xmin']}, {params['xmax']}]")
        print(f"时间: T={params['total_time']}")
        print(f"网格: Nx={params['Nx']}, Nt={params['Nt']}")
        print(f"对流速度: a={params['a']}")
        print(f"测试案例: {params['case']}")
        
        # 创建空间网格
        x = np.linspace(params['xmin'], params['xmax'], params['Nx']+1)
        
        # 计算初始条件
        print("\n正在计算初始条件...")
        u_initial = initial_condition(x, params['case'])
        
        # 计算精确解
        print("正在计算精确解...")
        u_exact = exact_solution(x, params['total_time'], params['a'], params['case'])
        
        # 计算Crank-Nicolson格式数值解
        print("\n计算Crank-Nicolson格式数值解...")
        u_cn = solve_pure_convection_cn(x, params['a'], params['total_time'], params['Nt'], u_initial)
        
        # 计算迎风格式数值解
        print("计算迎风格式数值解...")
        u_upwind = solve_pure_convection_upwind(x, params['a'], params['total_time'], params['Nt'], u_initial)
        
        # 计算误差
        error_cn = u_cn - u_exact
        error_upwind = u_upwind - u_exact
        
        errors = {
            'cn_max': np.max(np.abs(error_cn)),
            'cn_l2': np.sqrt(np.mean(error_cn**2)),
            'cn_mass': abs(np.trapz(u_cn, x) - np.trapz(u_exact, x)) / abs(np.trapz(u_exact, x)),
            'upwind_max': np.max(np.abs(error_upwind)),
            'upwind_l2': np.sqrt(np.mean(error_upwind**2)),
            'upwind_mass': abs(np.trapz(u_upwind, x) - np.trapz(u_exact, x)) / abs(np.trapz(u_exact, x))
        }
        
        # 绘制对比图
        print("\n正在生成图像...")
        plot_comparison(x, u_initial, u_cn, u_upwind, u_exact, params, errors)
        
        # 网格收敛性测试
        print("\n" + "="*60)
        print("网格收敛性测试")
        print("="*60)
        
        Nx_list = [50, 100, 200, 400]
        errors_cn = []
        errors_upwind = []
        
        for Nx_test in Nx_list:
            print(f"\n测试网格: Nx={Nx_test}")
            x_test = np.linspace(params['xmin'], params['xmax'], Nx_test+1)
            
            # 计算初始条件
            u0_test = initial_condition(x_test, params['case'])
            
            # 计算精确解
            u_exact_test = exact_solution(x_test, params['total_time'], params['a'], params['case'])
            
            # 调整时间步数以保持CFL数
            dx_test = (params['xmax'] - params['xmin']) / Nx_test
            CFL = params['a'] * params['total_time'] / params['Nt'] / dx_test
            Nt_test = int(params['Nt'] * (params['Nx'] / Nx_test))
            
            # 计算Crank-Nicolson解
            u_cn_test = solve_pure_convection_cn(x_test, params['a'], params['total_time'], Nt_test, u0_test)
            error_cn_test = np.sqrt(np.trapz((u_cn_test - u_exact_test)**2, x_test))
            errors_cn.append((Nx_test, error_cn_test))
            
            # 计算迎风格式解
            u_upwind_test = solve_pure_convection_upwind(x_test, params['a'], params['total_time'], Nt_test, u0_test)
            error_upwind_test = np.sqrt(np.trapz((u_upwind_test - u_exact_test)**2, x_test))
            errors_upwind.append((Nx_test, error_upwind_test))
            
            print(f"  Crank-Nicolson L2误差: {error_cn_test:.2e}")
            print(f"  迎风格式 L2误差: {error_upwind_test:.2e}")
        
        # 计算收敛阶
        if len(errors_cn) >= 2:
            print("\nCrank-Nicolson格式收敛阶:")
            for i in range(1, len(errors_cn)):
                Nx1, err1 = errors_cn[i-1]
                Nx2, err2 = errors_cn[i]
                order = np.log(err1/err2) / np.log(Nx2/Nx1)
                print(f"  {Nx1} → {Nx2}: 收敛阶 ≈ {order:.2f}")
        
        if len(errors_upwind) >= 2:
            print("\n迎风格式收敛阶:")
            for i in range(1, len(errors_upwind)):
                Nx1, err1 = errors_upwind[i-1]
                Nx2, err2 = errors_upwind[i]
                order = np.log(err1/err2) / np.log(Nx2/Nx1)
                print(f"  {Nx1} → {Nx2}: 收敛阶 ≈ {order:.2f}")
        
        print("\n" + "="*60)
        print("验证完成")
        print("="*60)
        
    except FileNotFoundError:
        print("错误: 找不到配置文件 config_convection.txt")
        print("正在创建默认配置文件...")
        
        # 创建默认配置文件
        default_config = """# 纯对流方程参数配置

[Domain]
xmin = -1.5
xmax = 1.5

[Time]
total_time = 1.0

[Mesh]
Nx = 200
Nt = 1000

[Physics]
a = 1.0

[Case]
type = combined_wave  # 可选: gaussian, combined_wave, square_wave
"""
        
        with open('config_convection.txt', 'w', encoding='utf-8') as f:
            f.write(default_config)
        
        print("已创建默认配置文件 config_convection.txt")
        print("请重新运行程序")
        
    except Exception as e:
        print(f"运行时错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()