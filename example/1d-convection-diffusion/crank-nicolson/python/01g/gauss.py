import numpy as np
import matplotlib.pyplot as plt
import configparser
from scipy import sparse
from scipy.sparse.linalg import spsolve

def read_config(filename='config.txt'):
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
        't0': float(config['Global']['t0'])
    }
    
    # 读取各案例参数
    cases = []
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown']
    
    for i, section in enumerate(config.sections()):
        if section.startswith('Case'):
            case_params = {
                'a': float(config[section]['a']),
                'nu': float(config[section]['nu']),
                'label': config[section]['label'],
                'color': colors[i-1] if i-1 < len(colors) else 'black'
            }
            cases.append(case_params)
    
    return global_params, cases

def solve_convection_diffusion(a, nu, global_params):
    """
    求解一维对流扩散方程
    Equation: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²
    """
    # 解包全局参数
    xmin = global_params['xmin']
    xmax = global_params['xmax']
    total_time = global_params['total_time']
    Nx = global_params['Nx']
    Nt = global_params['Nt']
    t0 = global_params['t0']
    
    # 网格设置
    dx = (xmax - xmin) / Nx
    dt = total_time / Nt
    x = np.linspace(xmin, xmax, Nx+1)
    
    # 初始条件：高斯分布
    if nu > 0:
        u = np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)
    else:
        # 纯对流情况，使用小的扩散系数以保持稳定性
        nu_small = 0.01
        u = np.exp(-x**2 / (4*nu_small*t0)) / np.sqrt(4*np.pi*nu_small*t0)
    
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
    
    return x, u

def exact_solution(x, t, a, nu, t0):
    """
    精确解
    u(x,t) = 1/√[4πν(t+t0)] * exp[-(x-at)²/(4ν(t+t0))]
    """
    if nu == 0:
        # 纯对流情况，使用小的扩散系数
        nu_small = 0.01
        return 1.0 / np.sqrt(4*np.pi*nu_small*(t+t0)) * np.exp(-(x-a*t)**2/(4*nu_small*(t+t0)))
    else:
        return 1.0 / np.sqrt(4*np.pi*nu*(t+t0)) * np.exp(-(x-a*t)**2/(4*nu*(t+t0)))

def plot_comparison(global_params, cases):
    """
    绘制所有案例的精确解与计算解对比图
    """
    # 根据案例数量确定子图布局
    n_cases = len(cases)
    if n_cases <= 3:
        nrows, ncols = 1, n_cases
    elif n_cases <= 6:
        nrows, ncols = 2, 3
    else:
        nrows, ncols = 2, 4
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    if n_cases == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    # 设置全局参数
    total_time = global_params['total_time']
    t0 = global_params['t0']
    
    # 存储误差信息用于汇总
    error_summary = []
    
    for idx, case in enumerate(cases):
        if idx >= len(axes):
            break
            
        ax = axes[idx]
        a = case['a']
        nu = case['nu']
        label = case['label']
        color = case['color']
        
        print(f"计算案例 {idx+1}: {label} (a={a}, ν={nu})")
        
        # 计算数值解
        x, u_num = solve_convection_diffusion(a, nu, global_params)
        
        # 计算精确解
        u_exact = exact_solution(x, total_time, a, nu, t0)
        
        # 计算误差
        error = np.abs(u_num - u_exact)
        max_error = np.max(error)
        rms_error = np.sqrt(np.mean(error**2))
        
        # 记录误差
        error_summary.append({
            'label': label,
            'max_error': max_error,
            'rms_error': rms_error
        })
        
        # 绘制对比图
        ax.plot(x, u_num, color=color, linestyle='-', linewidth=2, label='Numerical')
        ax.plot(x, u_exact, color=color, linestyle='--', linewidth=2, alpha=0.7, label='Exact')
        
        # 设置子图属性
        ax.set_xlabel('x')
        ax.set_ylabel('u(x,T)')
        ax.set_title(f'{label}\na={a}, ν={nu}')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # 在图上添加误差信息
        error_text = f'Max error: {max_error:.2e}\nRMS error: {rms_error:.2e}'
        ax.text(0.02, 0.98, error_text, transform=ax.transAxes, 
                fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 隐藏多余的子图
    for idx in range(len(cases), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    
    # 保存图像
    plt.savefig('comparison_results.png', dpi=300, bbox_inches='tight')
    print("\n对比图已保存为: comparison_results.png")
    
    # 显示图像
    plt.show()
    
    # 打印误差汇总
    print("\n" + "="*60)
    print("误差汇总:")
    print("="*60)
    for summary in error_summary:
        print(f"{summary['label']:25s} Max error: {summary['max_error']:.2e}  RMS error: {summary['rms_error']:.2e}")
    
    return error_summary

def main():
    """
    主函数：从配置文件读取参数并绘制对比图
    """
    print("="*60)
    print("对流扩散方程数值解与精确解对比")
    print("="*60)
    
    try:
        # 读取配置文件
        print("正在读取配置文件 config.txt...")
        global_params, cases = read_config('config.txt')
        
        print(f"找到 {len(cases)} 个测试案例")
        print(f"全局参数: x∈[{global_params['xmin']}, {global_params['xmax']}], "
              f"T={global_params['total_time']}, Nx={global_params['Nx']}, Nt={global_params['Nt']}")
        
        # 绘制对比图
        error_summary = plot_comparison(global_params, cases)
        
    except FileNotFoundError:
        print("错误: 找不到配置文件 config.txt")
        print("请确保当前目录下存在 config.txt 文件")
        print("正在使用默认参数运行...")
        
        # 使用默认参数
        global_params = {
            'xmin': -5.0,
            'xmax': 10.0,
            'total_time': 1.0,
            'Nx': 200,
            'Nt': 1000,
            't0': 0.1
        }
        
        cases = [
            {'a': 0.0, 'nu': 0.1, 'label': 'Pure Diffusion', 'color': 'blue'},
            {'a': 1.0, 'nu': 0.01, 'label': 'Pure Convection', 'color': 'red'},
            {'a': 1.0, 'nu': 0.1, 'label': 'Convection+Diffusion', 'color': 'green'}
        ]
        
        plot_comparison(global_params, cases)

if __name__ == "__main__":
    main()