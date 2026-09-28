import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve

def solve_convection_diffusion(a, nu, t0, xmin, xmax, T, Nx, Nt, 
                               scheme='CN'):
    """
    求解一维对流扩散方程
    
    方程: ∂u/∂t + a ∂u/∂x = ν ∂²u/∂x²
    
    参数:
    a: 对流速度
    nu: 扩散系数
    t0: 高斯初始分布的宽度参数
    xmin, xmax: 空间域边界
    T: 总时间
    Nx: 空间网格点数
    Nt: 时间步数
    scheme: 数值格式 ('CN' for Crank-Nicolson)
    """
    # 网格参数
    dx = (xmax - xmin) / Nx
    dt = T / Nt
    x = np.linspace(xmin, xmax, Nx+1)  # 包括边界点
    
    # Courant数和扩散数
    C = a * dt / dx  # Courant数
    d = nu * dt / (dx**2)  # 扩散数
    
    print(f"参数: a={a}, nu={nu}, t0={t0}")
    print(f"网格: dx={dx:.4f}, dt={dt:.4f}")
    print(f"CFL条件: C={C:.4f}, d={d:.4f}")
    
    # 初始条件
    u = np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)
    
    if scheme == 'CN':
        # Crank-Nicolson格式
        # 构造三对角矩阵（中心差分）
        
        # 内部点索引
        interior = np.arange(1, Nx)
        
        for n in range(Nt):
            # 显式部分（时间层 n）
            u_explicit = u.copy()
            
            # 使用周期边界条件（或零梯度边界）
            # 左边界
            u_explicit[0] = u[1]  # Neumann边界
            # 右边界
            u_explicit[Nx] = u[Nx-1]  # Neumann边界
            
            # 构造系数矩阵 A (隐式部分)
            diag_main = np.ones(Nx+1)  # 主对角线
            diag_lower = np.zeros(Nx+1)  # 下对角线
            diag_upper = np.zeros(Nx+1)  # 上对角线
            rhs = np.zeros(Nx+1)  # 右端项
            
            # 内部点
            for i in interior:
                # 对流项：中心差分
                conv_coef = a * dt / (4 * dx)
                # 扩散项：中心差分
                diff_coef = nu * dt / (2 * dx**2)
                
                # 隐式矩阵系数
                diag_main[i] = 1.0 + diff_coef * 2
                diag_lower[i] = -diff_coef + conv_coef
                diag_upper[i] = -diff_coef - conv_coef
                
                # 右端项（显式部分）
                rhs[i] = (1.0 - diff_coef * 2) * u[i] + \
                        (diff_coef + conv_coef) * u[i-1] + \
                        (diff_coef - conv_coef) * u[i+1]
            
            # 边界条件
            # Neumann边界：u_0 = u_1, u_N = u_{N-1}
            diag_main[0] = 1.0
            diag_upper[0] = -1.0  # u0 - u1 = 0
            rhs[0] = 0.0
            
            diag_main[Nx] = 1.0
            diag_lower[Nx] = -1.0  # uN - u_{N-1} = 0
            rhs[Nx] = 0.0
            
            # 构造三对角矩阵
            A = sparse.diags([diag_lower[1:], diag_main, diag_upper[:-1]], 
                            [-1, 0, 1], format='csr')
            
            # 求解线性系统
            u = spsolve(A, rhs)
    
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
    
    # 质量（积分）
    mass = np.sum(u) * dx
    
    # 峰值位置和高度
    idx_max = np.argmax(u)
    x_peak = x[idx_max]
    u_max = u[idx_max]
    
    # 解析解的峰值
    x_peak_exact = a * t
    u_max_exact = 1.0 / np.sqrt(4*np.pi*nu*(t+t0))
    
    # 标准差（分布宽度）
    if mass > 0:
        mean = np.sum(x * u) * dx / mass
        variance = np.sum((x - mean)**2 * u) * dx / mass
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
        'x_peak_exact': x_peak_exact,
        'u_max_exact': u_max_exact,
        'std': std,
        'std_exact': std_exact
    }

def main():
    # 参数设置
    a = 1.0      # 对流速度
    nu = 0.1     # 扩散系数
    t0 = 0.1     # 初始分布宽度参数
    
    # 空间域（需要足够大以包含高斯分布）
    xmin = -5.0
    xmax = 10.0
    Nx = 300     # 空间网格点数
    
    # 时间参数
    T = 2.0      # 总时间
    Nt = 1000    # 时间步数
    
    print("="*60)
    print("高斯初始分布的对流扩散算例")
    print("="*60)
    
    # 数值解
    x, u_num = solve_convection_diffusion(a, nu, t0, xmin, xmax, T, Nx, Nt)
    
    # 解析解
    u_exact = exact_solution_convection_diffusion(x, T, a, nu, t0)
    
    # 初始条件（用于比较）
    u_initial = exact_solution_convection_diffusion(x, 0, a, nu, t0)
    
    # 计算统计特征
    metrics = compute_metrics(x, u_num, a, nu, T, t0)
    
    # ==================== 绘图 ====================
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. 数值解与解析解对比
    axes[0, 0].plot(x, u_initial, 'g-', linewidth=1.5, alpha=0.7, label='Initial (t=0)')
    axes[0, 0].plot(x, u_num, 'b-', linewidth=2, label='Numerical (t=T)')
    axes[0, 0].plot(x, u_exact, 'r--', linewidth=2, label='Exact (t=T)')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('u(x,t)')
    axes[0, 0].set_title('Numerical vs Exact Solution')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    
    # 2. 误差分布
    error = u_num - u_exact
    axes[0, 1].plot(x, error, 'r-', linewidth=1.5)
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('Error (Numerical - Exact)')
    axes[0, 1].set_title(f'Absolute Error (max={np.max(np.abs(error)):.2e})')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. 峰值区域放大
    peak_region = np.abs(x - metrics['x_peak_exact']) < 2.0
    axes[1, 0].plot(x[peak_region], u_num[peak_region], 'b-o', markersize=4, label='Numerical')
    axes[1, 0].plot(x[peak_region], u_exact[peak_region], 'r--s', markersize=4, label='Exact')
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('u(x,t)')
    axes[1, 0].set_title('Peak Region (Zoomed)')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # 4. 对数坐标下的分布
    axes[1, 1].semilogy(x, u_initial, 'g-', alpha=0.7, label='Initial')
    axes[1, 1].semilogy(x, u_num, 'b-', label='Numerical')
    axes[1, 1].semilogy(x, u_exact, 'r--', label='Exact')
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('log u(x,t)')
    axes[1, 1].set_title('Solution in Log Scale')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    # ==================== 结果分析 ====================
    print("\n" + "="*60)
    print("结果验证")
    print("="*60)
    
    print(f"\n1. 质量守恒 (应接近1):")
    print(f"   数值解质量: {metrics['mass']:.6f}")
    print(f"   解析解质量: 1.000000")
    print(f"   相对误差: {abs(metrics['mass']-1):.2e}")
    
    print(f"\n2. 峰值位置 (对流效应):")
    print(f"   数值解峰值位置: {metrics['x_peak']:.6f}")
    print(f"   解析解峰值位置: {metrics['x_peak_exact']:.6f} (at)")
    print(f"   绝对误差: {abs(metrics['x_peak'] - metrics['x_peak_exact']):.2e}")
    print(f"   相对误差: {abs(metrics['x_peak'] - metrics['x_peak_exact'])/abs(metrics['x_peak_exact']):.2e}")
    
    print(f"\n3. 峰值高度 (扩散效应):")
    print(f"   数值解峰值高度: {metrics['u_max']:.6f}")
    print(f"   解析解峰值高度: {metrics['u_max_exact']:.6f}")
    print(f"   绝对误差: {abs(metrics['u_max'] - metrics['u_max_exact']):.2e}")
    print(f"   相对误差: {abs(metrics['u_max'] - metrics['u_max_exact'])/metrics['u_max_exact']:.2e}")
    
    print(f"\n4. 分布宽度 (扩散效应):")
    print(f"   数值解标准差: {metrics['std']:.6f}")
    print(f"   解析解标准差: {metrics['std_exact']:.6f}")
    print(f"   绝对误差: {abs(metrics['std'] - metrics['std_exact']):.2e}")
    print(f"   相对误差: {abs(metrics['std'] - metrics['std_exact'])/metrics['std_exact']:.2e}")
    
    print(f"\n5. 数值误差统计:")
    print(f"   最大绝对误差: {np.max(np.abs(error)):.2e}")
    print(f"   L2误差范数: {np.sqrt(np.sum(error**2)/len(error)):.2e}")
    print(f"   平均绝对误差: {np.mean(np.abs(error)):.2e}")
    
    # 检查CFL稳定性条件
    print(f"\n6. CFL稳定性检查:")
    print(f"   对流CFL数 (C): {a * (T/Nt) / ((xmax-xmin)/Nx):.4f}")
    print(f"   扩散CFL数 (d): {nu * (T/Nt) / (((xmax-xmin)/Nx)**2):.4f}")
    print(f"   建议: C < 1 且 d < 0.5 以确保稳定性")

if __name__ == "__main__":
    main()