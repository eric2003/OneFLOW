import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve

def solve_convection_diffusion_cn(a, nu, t0, xmin, xmax, T, Nx, Nt):
    """
    使用Crank-Nicolson格式求解一维对流扩散方程
    
    方程: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²
    
    参数:
    a: 对流速度
    nu: 扩散系数
    t0: 高斯初始分布的宽度参数
    xmin, xmax: 空间域边界
    T: 总时间
    Nx: 空间网格点数
    Nt: 时间步数
    """
    # 网格参数
    dx = (xmax - xmin) / Nx
    dt = T / Nt
    x = np.linspace(xmin, xmax, Nx+1)
    
    # Courant数和扩散数
    C = a * dt / dx
    d = nu * dt / (dx**2)
    
    print(f"参数: a={a}, nu={nu}, t0={t0}")
    print(f"网格: dx={dx:.4f}, dt={dt:.4f}")
    print(f"CFL: C={C:.4f}, d={d:.4f}")
    
    # 初始条件
    u = np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)
    
    # 时间步进
    for n in range(Nt):
        t_n = n * dt
        
        # 构造线性系统: (I + 0.5*dt*L)u^{n+1} = (I - 0.5*dt*L)u^n
        # L = -a*d/dx + nu*d²/dx²
        
        # 内部点索引 (1 到 Nx-1)
        N = Nx + 1  # 总点数
        interior = np.arange(1, Nx)
        
        # 系数矩阵 A (隐式部分: I + 0.5*dt*L)
        main_diag = np.ones(N)
        lower_diag = np.zeros(N)
        upper_diag = np.zeros(N)
        
        # 右端向量 b (显式部分: (I - 0.5*dt*L)u^n)
        b = np.zeros(N)
        
        # 内部点
        for i in interior:
            # 对流项系数 (中心差分)
            conv_coef = a / (2*dx)
            # 扩散项系数 (中心二阶差分)
            diff_coef = nu / (dx**2)
            
            # 隐式矩阵系数 (I + 0.5*dt*L)
            main_diag[i] = 1.0 + dt * diff_coef  # 主对角
            lower_diag[i] = -0.5 * dt * (conv_coef + diff_coef)  # 下对角
            upper_diag[i] = 0.5 * dt * (conv_coef - diff_coef)   # 上对角
            
            # 右端向量 ((I - 0.5*dt*L)u^n)
            b[i] = (1.0 - dt * diff_coef) * u[i] + \
                   0.5 * dt * (conv_coef + diff_coef) * u[i-1] + \
                   0.5 * dt * (-conv_coef + diff_coef) * u[i+1]
        
        # 边界条件: 使用周期边界（或零梯度）
        # 左边界 (i=0): u_0 = u_1
        main_diag[0] = 1.0
        upper_diag[0] = -1.0
        b[0] = 0.0
        
        # 右边界 (i=Nx): u_N = u_{N-1}
        main_diag[Nx] = 1.0
        lower_diag[Nx] = -1.0
        b[Nx] = 0.0
        
        # 构建稀疏矩阵
        A = sparse.diags([lower_diag[1:], main_diag, upper_diag[:-1]], 
                        [-1, 0, 1], format='csr')
        
        # 求解
        u_new = spsolve(A, b)
        u = u_new.copy()
    
    return x, u

def solve_convection_diffusion_upwind(a, nu, t0, xmin, xmax, T, Nx, Nt):
    """
    使用迎风格式（一阶精度）求解，更稳定但耗散较大
    """
    dx = (xmax - xmin) / Nx
    dt = T / Nt
    x = np.linspace(xmin, xmax, Nx+1)
    
    C = a * dt / dx
    d = nu * dt / (dx**2)
    
    print(f"使用迎风格式: C={C:.4f}, d={d:.4f}")
    
    # 初始条件
    u = np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)
    
    # 时间步进
    for n in range(Nt):
        u_old = u.copy()
        
        # 内部点
        for i in range(1, Nx):
            if a > 0:
                # 迎风格式：向后差分
                conv = a * (u_old[i] - u_old[i-1]) / dx
            else:
                # 如果a<0则用向前差分
                conv = a * (u_old[i+1] - u_old[i]) / dx
            
            # 扩散项：中心差分
            diff = nu * (u_old[i+1] - 2*u_old[i] + u_old[i-1]) / (dx**2)
            
            u[i] = u_old[i] - dt * conv + dt * diff
        
        # 边界条件（零梯度）
        u[0] = u[1]
        u[Nx] = u[Nx-1]
    
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
    dx = x[1] - x[0]
    
    # 质量（积分）- 使用trapezoid替代trapz
    try:
        mass = np.trapz(u, x)  # 旧版兼容
    except:
        mass = np.trapezoid(u, x)  # 新版numpy
    
    # 峰值位置和高度
    idx_max = np.argmax(u)
    x_peak = x[idx_max]
    u_max = u[idx_max]
    
    # 解析解的峰值
    x_peak_exact = a * t
    u_max_exact = 1.0 / np.sqrt(4*np.pi*nu*(t+t0))
    
    # 计算标准差
    if mass > 0:
        # 使用梯形积分计算一阶和二阶矩
        try:
            mean = np.trapz(x * u, x) / mass
            variance = np.trapz((x - mean)**2 * u, x) / mass
        except:
            mean = np.trapezoid(x * u, x) / mass
            variance = np.trapezoid((x - mean)**2 * u, x) / mass
        std = np.sqrt(variance)
    else:
        mean = 0
        std = 0
    
    # 解析解的标准差
    std_exact = np.sqrt(2*nu*(t+t0))
    
    return {
        'mass': mass,
        'x_peak': x_peak,
        'u_max': u_max,
        'mean': mean,
        'std': std,
        'x_peak_exact': x_peak_exact,
        'u_max_exact': u_max_exact,
        'std_exact': std_exact
    }

def main():
    # 参数设置
    a = 1.0      # 对流速度
    nu = 0.1     # 扩散系数
    t0 = 0.1     # 初始分布宽度参数
    
    # 空间域（需要足够大以包含移动后的分布）
    xmin = -5.0
    xmax = 10.0  # 注意：峰值会移动到 x = aT = 2.0 处
    Nx = 300
    
    # 时间参数
    T = 2.0
    Nt = 2000    # 增加时间步数以提高精度
    
    print("="*60)
    print("高斯初始分布的对流扩散算例")
    print("="*60)
    
    # 方法1：Crank-Nicolson格式（中心差分）
    print("\n--- 方法1: Crank-Nicolson格式 ---")
    x_cn, u_cn = solve_convection_diffusion_cn(a, nu, t0, xmin, xmax, T, Nx, Nt)
    
    # 方法2：迎风格式
    print("\n--- 方法2: 迎风格式 ---")
    x_up, u_up = solve_convection_diffusion_upwind(a, nu, t0, xmin, xmax, T, Nx, Nt)
    
    # 解析解
    u_exact = exact_solution_convection_diffusion(x_cn, T, a, nu, t0)
    u_initial = exact_solution_convection_diffusion(x_cn, 0, a, nu, t0)
    
    # 计算统计特征
    metrics_cn = compute_metrics(x_cn, u_cn, a, nu, T, t0)
    metrics_up = compute_metrics(x_up, u_up, a, nu, T, t0)
    
    # ==================== 绘图 ====================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. 不同方法的对比
    axes[0, 0].plot(x_cn, u_initial, 'g-', linewidth=1.5, alpha=0.7, label='Initial (t=0)')
    axes[0, 0].plot(x_cn, u_cn, 'b-', linewidth=2, label=f'Crank-Nicolson')
    axes[0, 0].plot(x_up, u_up, 'm-', linewidth=2, alpha=0.8, label='Upwind')
    axes[0, 0].plot(x_cn, u_exact, 'r--', linewidth=2, label='Exact (t=T)')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('u(x,t)')
    axes[0, 0].set_title(f'Comparison of Numerical Methods (a={a}, ν={nu}, T={T})')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. 误差分布
    error_cn = u_cn - u_exact
    error_up = np.interp(x_cn, x_up, u_up) - u_exact
    
    axes[0, 1].plot(x_cn, error_cn, 'b-', linewidth=1.5, label='Crank-Nicolson Error')
    axes[0, 1].plot(x_cn, error_up, 'm-', linewidth=1.5, label='Upwind Error')
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('Error (Numerical - Exact)')
    axes[0, 1].set_title('Error Distribution')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. 峰值区域放大
    peak_window = 1.0
    mask = np.abs(x_cn - metrics_cn['x_peak_exact']) < peak_window
    axes[1, 0].plot(x_cn[mask], u_cn[mask], 'b-o', markersize=4, label='Crank-Nicolson')
    axes[1, 0].plot(x_cn[mask], u_exact[mask], 'r--s', markersize=4, label='Exact')
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('u(x,t)')
    axes[1, 0].set_title(f'Peak Region Zoom (x ≈ {metrics_cn["x_peak_exact"]:.2f})')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. 演化过程动画（关键时间点）
    axes[1, 1].plot(x_cn, u_initial, 'g-', alpha=0.5, label='t=0')
    
    # 计算中间时间点的解析解
    times = [T/4, T/2, 3*T/4, T]
    colors = ['orange', 'purple', 'brown', 'red']
    linestyles = ['--', '--', '--', '--']
    
    for t_val, color, ls in zip(times, colors, linestyles):
        u_mid = exact_solution_convection_diffusion(x_cn, t_val, a, nu, t0)
        axes[1, 1].plot(x_cn, u_mid, color=color, linestyle=ls, alpha=0.5, label=f't={t_val:.1f}')
    
    axes[1, 1].plot(x_cn, u_exact, 'r-', linewidth=2, label=f't={T} (exact)')
    axes[1, 1].plot(x_cn, u_cn, 'b-', linewidth=2, alpha=0.8, label=f't={T} (Crank-Nicolson)')
    
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('u(x,t)')
    axes[1, 1].set_title('Time Evolution (Analytical)')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # ==================== 结果分析 ====================
    print("\n" + "="*60)
    print("结果验证")
    print("="*60)
    
    print(f"\n1. Crank-Nicolson格式结果:")
    print(f"   质量: {metrics_cn['mass']:.6f} (理论值: 1.000000)")
    print(f"   峰值位置: {metrics_cn['x_peak']:.6f} (理论: {metrics_cn['x_peak_exact']:.6f})")
    print(f"   峰值高度: {metrics_cn['u_max']:.6f} (理论: {metrics_cn['u_max_exact']:.6f})")
    print(f"   分布宽度(std): {metrics_cn['std']:.6f} (理论: {metrics_cn['std_exact']:.6f})")
    
    print(f"\n2. 迎风格式结果:")
    print(f"   质量: {metrics_up['mass']:.6f} (理论值: 1.000000)")
    print(f"   峰值位置: {metrics_up['x_peak']:.6f} (理论: {metrics_up['x_peak_exact']:.6f})")
    print(f"   峰值高度: {metrics_up['u_max']:.6f} (理论: {metrics_up['u_max_exact']:.6f})")
    print(f"   分布宽度(std): {metrics_up['std']:.6f} (理论: {metrics_up['std_exact']:.6f})")
    
    print(f"\n3. 误差分析:")
    print(f"   Crank-Nicolson L2误差: {np.sqrt(np.mean(error_cn**2)):.2e}")
    print(f"   迎风格式 L2误差: {np.sqrt(np.mean(error_up**2)):.2e}")
    
    print(f"\n4. 对流效应验证:")
    print(f"   理论位移: x = aT = {a} * {T} = {a*T:.2f}")
    print(f"   CN格式位移: {metrics_cn['x_peak']:.2f}")
    print(f"   迎风格式位移: {metrics_up['x_peak']:.2f}")
    
    print(f"\n5. 扩散效应验证:")
    print(f"   理论宽度增长: σ = √[2ν(t+t0)]")
    print(f"   初始宽度: √[2ν*t0] = {np.sqrt(2*nu*t0):.4f}")
    print(f"   最终宽度: √[2ν*(T+t0)] = {np.sqrt(2*nu*(T+t0)):.4f}")
    print(f"   CN格式宽度: {metrics_cn['std']:.4f}")
    print(f"   迎风格式宽度: {metrics_up['std']:.4f}")
    
    # 测试不同参数
    print("\n" + "="*60)
    print("参数敏感性测试")
    print("="*60)
    
    test_cases = [
        (2.0, 0.1, "强对流"),
        (0.5, 0.1, "中等对流"),
        (1.0, 0.01, "小扩散"),
        (1.0, 0.5, "大扩散"),
    ]
    
    for a_test, nu_test, desc in test_cases:
        print(f"\n--- {desc} (a={a_test}, ν={nu_test}) ---")
        x_test, u_test = solve_convection_diffusion_cn(a_test, nu_test, t0, xmin, xmax, T, Nx, Nt)
        u_exact_test = exact_solution_convection_diffusion(x_test, T, a_test, nu_test, t0)
        metrics_test = compute_metrics(x_test, u_test, a_test, nu_test, T, t0)
        
        print(f"   峰值位置误差: {abs(metrics_test['x_peak'] - a_test*T):.4f}")
        print(f"   峰值高度误差: {abs(metrics_test['u_max'] - 1/np.sqrt(4*np.pi*nu_test*(T+t0))):.4f}")

if __name__ == "__main__":
    main()