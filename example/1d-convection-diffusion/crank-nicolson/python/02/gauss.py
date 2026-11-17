import numpy as np
import matplotlib.pyplot as plt
import configparser
from scipy import sparse
from scipy.sparse.linalg import spsolve

def read_config(filename='config_advanced.txt'):
    """从配置文件读取参数"""
    config = configparser.ConfigParser()
    config.read(filename, encoding='utf-8')
    
    # 读取全局参数
    global_params = {
        'xmin': float(config['Global']['xmin']),
        'xmax': float(config['Global']['xmax']),
        'total_time': float(config['Global']['total_time']),
        'Nx': int(config['Global']['Nx']),
        'Nt': int(config['Global']['Nt']),
        't0': float(config['Global']['t0']),
        'initial_condition': config['Global']['initial_condition']
    }
    
    # 读取初始条件参数
    init_params = {
        'beta': float(config['InitialCondition']['beta']),
        'z': float(config['InitialCondition']['z']),
        'delta': float(config['InitialCondition']['delta']),
        'alpha': float(config['InitialCondition']['alpha']),
        'a_ellipse': float(config['InitialCondition']['a'])  # 重命名避免冲突
    }
    
    # 读取物理参数
    physics_params = {
        'a': float(config['Physics']['a']),  # 对流速度
        'nu': float(config['Physics']['nu'])  # 扩散系数
    }
    
    # 读取输出设置
    output_params = {
        'plot_title': config['Output']['plot_title'],
        'save_figure': config.getboolean('Output', 'save_figure'),
        'figure_name': config['Output']['figure_name']
    }
    
    return global_params, init_params, physics_params, output_params

def combined_wave_initial(x, params):
    """
    组合波形初始条件
    
    u₀(x) = 
    1) 1/6*(G(x,β,z-δ) + G(x,β,z+δ) + 4G(x,β,z)),  -0.8 ≤ x ≤ -0.6
    2) 1,                                        -0.4 ≤ x ≤ -0.2
    3) 1 - |10(x-0.1)|,                          0 ≤ x ≤ 0.2
    4) 1/6*(F(x,α,a-δ) + F(x,α,a+δ) + 4F(x,α,a)), 0.4 ≤ x ≤ 0.6
    5) 0,                                         otherwise
    """
    beta = params['beta']
    z = params['z']
    delta = params['delta']
    alpha = params['alpha']
    a_ellipse = params['a_ellipse']
    
    u0 = np.zeros_like(x)
    
    # 辅助函数
    def G(x, beta, z):
        """高斯函数"""
        return np.exp(-beta * (x - z)**2)
    
    def F(x, alpha, a):
        """半椭圆函数"""
        return np.sqrt(np.maximum(1 - alpha**2 * (x - a)**2, 0))
    
    # 区域1: 高斯波包 [-0.8, -0.6]
    mask1 = (-0.8 <= x) & (x <= -0.6)
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
    u0[mask4] = (1/6) * (F(x[mask4], alpha, a_ellipse-delta) + 
                        F(x[mask4], alpha, a_ellipse+delta) + 
                        4*F(x[mask4], alpha, a_ellipse))
    
    # 区域5: 其他区域为0（已初始化为0）
    
    return u0

def gaussian_initial(x, params):
    """高斯初始条件（备用）"""
    nu = params.get('nu', 0.1)
    t0 = params.get('t0', 0.1)
    return np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)

def solve_convection_diffusion(x, a, nu, total_time, Nt, u0):
    """
    求解一维对流扩散方程
    Equation: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²
    
    Parameters:
    x: 空间网格
    a: 对流速度
    nu: 扩散系数
    total_time: 总时间
    Nt: 时间步数
    u0: 初始条件
    """
    Nx = len(x) - 1
    dx = x[1] - x[0]
    dt = total_time / Nt
    
    u = u0.copy()
    
    # Crank-Nicolson格式时间推进
    for n in range(Nt):
        N = Nx + 1
        main_diag = np.ones(N)
        lower_diag = np.zeros(N)
        upper_diag = np.zeros(N)
        b = np.zeros(N)
        
        # 内部点
        for i in range(1, Nx):
            # 对流项系数（中心差分）
            conv_coef = a / (2*dx) if a != 0 else 0.0
            # 扩散项系数（中心二阶差分）
            diff_coef = nu / (dx**2) if nu != 0 else 0.0
            
            # 隐式矩阵系数
            main_diag[i] = 1.0 + dt * diff_coef
            lower_diag[i] = -0.5 * dt * (conv_coef + diff_coef)
            upper_diag[i] = 0.5 * dt * (conv_coef - diff_coef)
            
            # 右端向量
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
        
        # 构建并求解稀疏矩阵
        A = sparse.diags([lower_diag[1:], main_diag, upper_diag[:-1]], 
                        [-1, 0, 1], format='csr')
        u = spsolve(A, b)
    
    return u

def pure_convection_solution(x, t, a, u0_func, params):
    """
    纯对流方程的精确解（沿着特征线平移）
    u(x,t) = u₀(x-at)
    
    注意：这只适用于nu=0的纯对流情况
    """
    # 计算在x-at位置上的初始条件
    x_shifted = x - a * t
    return u0_func(x_shifted, params)

def plot_results(x, u_initial, u_numerical, u_exact, params, output_params):
    """
    绘制结果对比图
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 参数信息
    a = params['a']
    nu = params['nu']
    total_time = params['total_time']
    
    # 1. 初始条件和最终数值解
    axes[0, 0].plot(x, u_initial, 'g-', linewidth=2, alpha=0.7, label='Initial (t=0)')
    axes[0, 0].plot(x, u_numerical, 'b-', linewidth=2, label=f'Numerical (t={total_time})')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('u(x,t)')
    axes[0, 0].set_title(f'Combined Wave: Initial vs Numerical\nConvection a={a}, Diffusion ν={nu}')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. 数值解与精确解对比（如果可用）
    if u_exact is not None:
        axes[0, 1].plot(x, u_numerical, 'b-', linewidth=2, label='Numerical')
        axes[0, 1].plot(x, u_exact, 'r--', linewidth=2, alpha=0.7, label='Exact')
        axes[0, 1].set_xlabel('x')
        axes[0, 1].set_ylabel('u(x,t)')
        axes[0, 1].set_title('Numerical vs Exact Solution')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # 计算误差
        error = u_numerical - u_exact
        max_error = np.max(np.abs(error))
        rms_error = np.sqrt(np.mean(error**2))
        
        # 在图上添加误差信息
        error_text = f'Max error: {max_error:.2e}\nRMS error: {rms_error:.2e}'
        axes[0, 1].text(0.02, 0.98, error_text, transform=axes[0, 1].transAxes,
                       fontsize=10, verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    else:
        axes[0, 1].text(0.5, 0.5, 'No exact solution available\nfor diffusion terms',
                       ha='center', va='center', fontsize=12)
        axes[0, 1].set_title('No Exact Solution Available')
    
    # 3. 误差分布（如果精确解可用）
    if u_exact is not None:
        axes[1, 0].plot(x, error, 'r-', linewidth=1.5)
        axes[1, 0].set_xlabel('x')
        axes[1, 0].set_ylabel('Error (Numerical - Exact)')
        axes[1, 0].set_title(f'Error Distribution (max={max_error:.2e})')
        axes[1, 0].grid(True, alpha=0.3)
    
    # 4. 不同时间的演化（数值解）
    # 可以计算中间时间点的数值解来展示演化过程
    axes[1, 1].plot(x, u_initial, 'g-', alpha=0.5, label='t=0')
    axes[1, 1].plot(x, u_numerical, 'b-', linewidth=2, label=f't={total_time}')
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('u(x,t)')
    axes[1, 1].set_title('Time Evolution (Numerical)')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    # 设置总标题
    plt.suptitle(output_params['plot_title'], fontsize=16, y=1.02)
    plt.tight_layout()
    
    # 保存图像
    if output_params['save_figure']:
        plt.savefig(output_params['figure_name'], dpi=300, bbox_inches='tight')
        print(f"图像已保存为: {output_params['figure_name']}")
    
    plt.show()

def main():
    """
    主函数：使用组合波形初始条件进行测试
    """
    print("="*60)
    print("组合波形初始条件的对流扩散方程数值模拟")
    print("="*60)
    
    try:
        # 读取配置文件
        print("正在读取配置文件 config_advanced.txt...")
        global_params, init_params, physics_params, output_params = read_config('config_advanced.txt')
        
        # 合并参数
        all_params = {**global_params, **init_params, **physics_params}
        
        print(f"计算域: x∈[{global_params['xmin']}, {global_params['xmax']}]")
        print(f"时间: T={global_params['total_time']}, 网格: Nx={global_params['Nx']}, Nt={global_params['Nt']}")
        print(f"对流速度: a={physics_params['a']}, 扩散系数: ν={physics_params['nu']}")
        print(f"初始条件类型: {global_params['initial_condition']}")
        
        # 创建空间网格
        x = np.linspace(global_params['xmin'], global_params['xmax'], global_params['Nx']+1)
        
        # 设置初始条件
        if global_params['initial_condition'] == 'combined_wave':
            u0_func = combined_wave_initial
        elif global_params['initial_condition'] == 'gaussian':
            u0_func = gaussian_initial
        else:
            print(f"警告: 未知的初始条件类型 '{global_params['initial_condition']}'，使用组合波形")
            u0_func = combined_wave_initial
        
        u_initial = u0_func(x, all_params)
        
        # 计算数值解
        print("\n正在计算数值解...")
        u_numerical = solve_convection_diffusion(
            x, physics_params['a'], physics_params['nu'],
            global_params['total_time'], global_params['Nt'],
            u_initial
        )
        
        # 计算精确解（如果可能）
        u_exact = None
        if physics_params['nu'] == 0:
            print("注意: ν=0，这是纯对流方程，可以使用特征线法得到精确解")
            u_exact = pure_convection_solution(
                x, global_params['total_time'], physics_params['a'],
                u0_func, all_params
            )
        else:
            print("注意: ν≠0，这是对流扩散方程，没有简单的解析解")
        
        # 绘制结果
        print("\n正在生成图像...")
        plot_results(x, u_initial, u_numerical, u_exact, all_params, output_params)
        
        # 打印总结信息
        print("\n" + "="*60)
        print("模拟完成总结:")
        print("="*60)
        print(f"初始条件: {global_params['initial_condition']}")
        print(f"对流速度: a={physics_params['a']}")
        print(f"扩散系数: ν={physics_params['nu']}")
        print(f"模拟时间: T={global_params['total_time']}")
        
        if physics_params['nu'] == 0 and u_exact is not None:
            error = u_numerical - u_exact
            max_error = np.max(np.abs(error))
            rms_error = np.sqrt(np.mean(error**2))
            print(f"最大绝对误差: {max_error:.2e}")
            print(f"RMS误差: {rms_error:.2e}")
        else:
            print("误差分析: 无精确解可用于误差计算")
        
    except FileNotFoundError:
        print("错误: 找不到配置文件 config_advanced.txt")
        print("请确保当前目录下存在 config_advanced.txt 文件")
    except Exception as e:
        print(f"运行时错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()