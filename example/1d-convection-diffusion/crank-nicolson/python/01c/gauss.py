import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve

def solve_convection_diffusion_cn(a, nu, t0, xmin, xmax, T, Nx, Nt):
    """
    使用Crank-Nicolson格式求解一维对流扩散方程
    方程: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²
    """
    dx = (xmax - xmin) / Nx
    dt = T / Nt
    x = np.linspace(xmin, xmax, Nx+1)
    
    # 初始条件
    u = np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)
    
    # 时间步进
    for n in range(Nt):
        # 构造线性系统: (I + 0.5*dt*L)u^{n+1} = (I - 0.5*dt*L)u^n
        
        N = Nx + 1  # 总点数
        main_diag = np.ones(N)
        lower_diag = np.zeros(N)
        upper_diag = np.zeros(N)
        b = np.zeros(N)
        
        # 内部点 (1 到 Nx-1)
        for i in range(1, Nx):
            # 对流项系数 (中心差分)
            conv_coef = a / (2*dx)
            # 扩散项系数 (中心二阶差分)
            diff_coef = nu / (dx**2)
            
            # 隐式矩阵系数 (I + 0.5*dt*L)
            main_diag[i] = 1.0 + dt * diff_coef
            lower_diag[i] = -0.5 * dt * (conv_coef + diff_coef)
            upper_diag[i] = 0.5 * dt * (conv_coef - diff_coef)
            
            # 右端向量 ((I - 0.5*dt*L)u^n)
            b[i] = (1.0 - dt * diff_coef) * u[i] + \
                   0.5 * dt * (conv_coef + diff_coef) * u[i-1] + \
                   0.5 * dt * (-conv_coef + diff_coef) * u[i+1]
        
        # 边界条件: 零梯度 (Neumann边界)
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

def exact_solution_convection_diffusion(x, t, a, nu, t0):
    """
    对流扩散方程的解析解
    u(x,t) = 1/√[4πν(t+t0)] * exp[-(x-at)²/(4ν(t+t0))]
    """
    return 1.0 / np.sqrt(4*np.pi*nu*(t+t0)) * np.exp(-(x-a*t)**2/(4*nu*(t+t0)))

def compute_metrics(x, u, a, nu, t, t0):
    """
    计算数值解的统计特征
    """
    # 使用trapezoid进行积分（避免弃用警告）
    if hasattr(np, 'trapezoid'):
        mass = np.trapezoid(u, x)
        mean = np.trapezoid(x * u, x) / mass
        variance = np.trapezoid((x - mean)**2 * u, x) / mass
    else:
        mass = np.trapz(u, x)
        mean = np.trapz(x * u, x) / mass
        variance = np.trapz((x - mean)**2 * u, x) / mass
    
    # 峰值位置和高度
    idx_max = np.argmax(u)
    x_peak = x[idx_max]
    u_max = u[idx_max]
    
    # 解析解的理论值
    x_peak_exact = a * t
    u_max_exact = 1.0 / np.sqrt(4*np.pi*nu*(t+t0))
    std_exact = np.sqrt(2*nu*(t+t0))
    
    return {
        'mass': mass,
        'x_peak': x_peak,
        'u_max': u_max,
        'std': np.sqrt(variance),
        'x_peak_exact': x_peak_exact,
        'u_max_exact': u_max_exact,
        'std_exact': std_exact
    }

def run_convergence_test(a=1.0, nu=0.1, t0=0.1, T=2.0, xmin=-5.0, xmax=10.0):
    """运行网格收敛性测试"""
    print("\n" + "="*60)
    print("网格收敛性测试")
    print("="*60)
    
    Nx_list = [50, 100, 200, 400]
    errors = []
    
    # 用于插值的参考解（最细网格）
    Nx_ref = 800
    Nt_ref = Nx_ref * 10
    x_ref, u_ref = solve_convection_diffusion_cn(a, nu, t0, xmin, xmax, T, Nx_ref, Nt_ref)
    u_exact_ref = exact_solution_convection_diffusion(x_ref, T, a, nu, t0)
    
    for Nx in Nx_list:
        Nt = Nx * 5  # 保持CFL数大致不变
        dx = (xmax - xmin) / Nx
        
        # 计算数值解
        x, u_num = solve_convection_diffusion_cn(a, nu, t0, xmin, xmax, T, Nx, Nt)
        
        # 计算解析解在相同网格上
        u_exact = exact_solution_convection_diffusion(x, T, a, nu, t0)
        
        # 计算L2误差
        error = np.sqrt(np.trapezoid((u_num - u_exact)**2, x))
        errors.append(error)
        
        print(f"  Nx={Nx:3d}, Nt={Nt:4d}, dx={dx:.4f}, Error={error:.2e}")
    
    # 计算收敛阶
    if len(errors) >= 2:
        print("\n收敛阶估计:")
        for i in range(1, len(errors)):
            ratio = errors[i-1] / errors[i]
            order = np.log2(ratio) / np.log2(Nx_list[i]/Nx_list[i-1])
            print(f"  Nx {Nx_list[i-1]} → {Nx_list[i]}: 阶数 ≈ {order:.2f}")
    
    # 绘制收敛图
    plt.figure(figsize=(8, 6))
    dx_list = [(xmax - xmin)/Nx for Nx in Nx_list]
    
    plt.loglog(dx_list, errors, 'bo-', linewidth=2, markersize=8, label='Crank-Nicolson')
    
    # 添加参考线（一阶和二阶收敛）
    dx_min, dx_max = min(dx_list), max(dx_list)
    error_ref = errors[0]
    
    # 一阶收敛参考线
    x_ref_1 = np.linspace(dx_min, dx_max, 10)
    y_ref_1 = error_ref * (x_ref_1 / dx_list[0])
    plt.loglog(x_ref_1, y_ref_1, 'r--', label='一阶收敛 (slope=1)')
    
    # 二阶收敛参考线
    y_ref_2 = error_ref * (x_ref_1 / dx_list[0])**2
    plt.loglog(x_ref_1, y_ref_2, 'g--', label='二阶收敛 (slope=2)')
    
    plt.xlabel('网格间距 dx')
    plt.ylabel('L2误差')
    plt.title('Crank-Nicolson格式收敛性测试')
    plt.legend()
    plt.grid(True, alpha=0.3, which='both')
    plt.show()
    
    return errors, dx_list

def main():
    # 参数设置
    a = 1.0      # 对流速度
    nu = 0.1     # 扩散系数
    t0 = 0.1     # 初始分布宽度参数
    T = 2.0      # 总时间
    
    # 计算域
    xmin, xmax = -5.0, 10.0
    Nx, Nt = 300, 2000
    
    print("="*60)
    print("高斯初始分布的对流扩散算例验证")
    print(f"参数: a={a}, ν={nu}, t0={t0}, T={T}")
    print("="*60)
    
    # 数值解
    x, u_num = solve_convection_diffusion_cn(a, nu, t0, xmin, xmax, T, Nx, Nt)
    
    # 解析解
    u_exact = exact_solution_convection_diffusion(x, T, a, nu, t0)
    u_initial = exact_solution_convection_diffusion(x, 0, a, nu, t0)
    
    # 计算误差和统计特征
    error = u_num - u_exact
    metrics = compute_metrics(x, u_num, a, nu, T, t0)
    
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号    
    
    # ==================== 绘图 ====================
    plt.figure(figsize=(12, 8))
    
    # 主图：对比
    plt.subplot(2, 2, 1)
    plt.plot(x, u_initial, 'g-', alpha=0.7, label='Initial (t=0)')
    plt.plot(x, u_num, 'b-', linewidth=2, label='Numerical (t=T)')
    plt.plot(x, u_exact, 'r--', linewidth=2, label='Exact (t=T)')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.title(f'Convection-Diffusion: a={a}, ν={nu}, T={T}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 误差分布
    plt.subplot(2, 2, 2)
    plt.plot(x, error, 'r-', linewidth=1.5)
    plt.xlabel('x')
    plt.ylabel('Error (Numerical - Exact)')
    plt.title(f'Absolute Error (max={np.max(np.abs(error)):.2e})')
    plt.grid(True, alpha=0.3)
    
    # 峰值区域放大
    plt.subplot(2, 2, 3)
    peak_region = np.abs(x - metrics['x_peak_exact']) < 1.0
    plt.plot(x[peak_region], u_num[peak_region], 'b-o', markersize=4, label='Numerical')
    plt.plot(x[peak_region], u_exact[peak_region], 'r--s', markersize=4, label='Exact')
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.title(f'Peak Region (x ≈ {metrics["x_peak_exact"]:.2f})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 对数坐标
    plt.subplot(2, 2, 4)
    plt.semilogy(x, u_initial, 'g-', alpha=0.7, label='Initial')
    plt.semilogy(x, u_num, 'b-', label='Numerical')
    plt.semilogy(x, u_exact, 'r--', label='Exact')
    plt.xlabel('x')
    plt.ylabel('log u(x,t)')
    plt.title('Solution in Log Scale')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # ==================== 结果分析 ====================
    print("\n" + "="*60)
    print("验证结果")
    print("="*60)
    
    print(f"\n1. 质量守恒验证:")
    print(f"   数值解质量: {metrics['mass']:.6f}")
    print(f"   理论值: 1.000000")
    print(f"   相对误差: {abs(metrics['mass']-1):.2e}")
    
    print(f"\n2. 对流效应验证 (峰值位置):")
    print(f"   理论位移: x = aT = {a} * {T} = {a*T:.4f}")
    print(f"   数值解峰值位置: {metrics['x_peak']:.4f}")
    print(f"   绝对误差: {abs(metrics['x_peak'] - metrics['x_peak_exact']):.4f}")
    
    print(f"\n3. 扩散效应验证:")
    print(f"   峰值高度 (耗散):")
    print(f"     理论值: {metrics['u_max_exact']:.6f}")
    print(f"     数值解: {metrics['u_max']:.6f}")
    print(f"     相对误差: {abs(metrics['u_max'] - metrics['u_max_exact'])/metrics['u_max_exact']:.2e}")
    
    print(f"\n   分布宽度 (扩散):")
    print(f"     理论宽度: σ = √[2ν(T+t0)] = {metrics['std_exact']:.6f}")
    print(f"     数值宽度: {metrics['std']:.6f}")
    print(f"     相对误差: {abs(metrics['std'] - metrics['std_exact'])/metrics['std_exact']:.2e}")
    
    print(f"\n4. 误差统计:")
    print(f"   最大绝对误差: {np.max(np.abs(error)):.2e}")
    print(f"   RMS误差: {np.sqrt(np.mean(error**2)):.2e}")
    print(f"   平均绝对误差: {np.mean(np.abs(error)):.2e}")
    
    # 运行收敛性测试
    run_convergence_test(a, nu, t0, T, xmin, xmax)
    
    # 不同参数的测试
    print("\n" + "="*60)
    print("不同参数下的性能测试")
    print("="*60)
    
    test_params = [
        ("强对流", 2.0, 0.1),
        ("中等对流", 0.5, 0.1),
        ("小扩散", 1.0, 0.01),
        ("大扩散", 1.0, 0.5),
    ]
    
    for name, a_test, nu_test in test_params:
        print(f"\n--- {name} (a={a_test}, ν={nu_test}) ---")
        
        # 调整网格以适应不同参数
        if nu_test < 0.05:  # 小扩散需要更细网格
            Nx_test, Nt_test = 500, 3000
        else:
            Nx_test, Nt_test = 200, 1000
        
        x_test, u_test = solve_convection_diffusion_cn(a_test, nu_test, t0, xmin, xmax, T, Nx_test, Nt_test)
        u_exact_test = exact_solution_convection_diffusion(x_test, T, a_test, nu_test, t0)
        
        metrics_test = compute_metrics(x_test, u_test, a_test, nu_test, T, t0)
        error_test = u_test - u_exact_test
        
        print(f"   峰值位置误差: {abs(metrics_test['x_peak'] - a_test*T):.4f}")
        print(f"   峰值高度相对误差: {abs(metrics_test['u_max'] - metrics_test['u_max_exact'])/metrics_test['u_max_exact']:.2e}")
        print(f"   RMS误差: {np.sqrt(np.mean(error_test**2)):.2e}")

if __name__ == "__main__":
    main()