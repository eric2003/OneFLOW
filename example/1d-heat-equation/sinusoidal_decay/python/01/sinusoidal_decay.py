import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation(alpha, L, T, Nx, Nt, initial_condition='gaussian', 
                       boundary_left='neumann', boundary_right='neumann'):
    """
    求解一维热传导方程
    
    参数:
    alpha: 扩散系数
    L: 空间域长度 [0, L]
    T: 总时间
    Nx: 空间网格点数
    Nt: 时间步数
    initial_condition: 'gaussian' 或 'sine'
    boundary_left: 左边界条件 'neumann' 或 'dirichlet'
    boundary_right: 右边界条件 'neumann' 或 'dirichlet'
    """
    dx = L / Nx
    dt = T / Nt
    x = np.linspace(0, L, Nx+1)
    r = alpha * dt / (dx**2)
    
    # 设置初始条件
    if initial_condition == 'gaussian':
        u = np.exp(-x**2)  # 高斯脉冲
    elif initial_condition == 'sine':
        u = -np.sin(np.pi * x / L)  # 正弦函数，在x=L处为0
    else:
        raise ValueError("初始条件必须是 'gaussian' 或 'sine'")
    
    # 三对角矩阵系数
    a = np.full(Nx-1, -r/2.0)  # 下对角线
    b = np.full(Nx-1, 1.0 + r) # 主对角线  
    c = np.full(Nx-1, -r/2.0)  # 上对角线
    
    # 时间推进
    for _ in range(Nt):
        un = u.copy()
        d = np.zeros(Nx-1)  # 右端项
        
        # 构建右端项
        for i in range(1, Nx):
            j = i - 1
            if i == 1:  # 第一个内点
                term = (1.0 - r) * un[i] + (r/2.0) * un[i+1]
                if boundary_left == 'neumann':
                    term += (r/2.0) * un[i-1]  # Neumann边界：u_0 ≈ u_1
                d[j] = term
            elif i == Nx-1:  # 最后一个内点
                term = (r/2.0) * un[i-1] + (1.0 - r) * un[i]
                if boundary_right == 'neumann':
                    term += (r/2.0) * un[i+1]  # Neumann边界：u_N ≈ u_{N-1}
                d[j] = term
            else:  # 内部点
                d[j] = (r/2.0) * un[i-1] + (1.0 - r) * un[i] + (r/2.0) * un[i+1]
        
        # 应用边界条件修正
        if boundary_left == 'neumann':
            # u[0] = u[1] 的隐式处理：修正第一个方程
            d[0] += (r/2.0) * un[0]
        elif boundary_left == 'dirichlet':
            # Dirichlet边界：u[0]固定，修正第一个方程
            d[0] += (r/2.0) * un[0]
        
        # 追赶法求解三对角系统
        def thomas_solve(a, b, c, d):
            n = len(d)
            c_prime = np.zeros(n-1)
            d_prime = np.zeros(n)
            
            # 前向消元
            c_prime[0] = c[0] / b[0]
            d_prime[0] = d[0] / b[0]
            
            for i in range(1, n):
                if i < n-1:
                    denom = b[i] - a[i] * c_prime[i-1]
                    c_prime[i] = c[i] / denom
                denom = b[i] - a[i] * c_prime[i-1] if i > 0 else b[0]
                d_prime[i] = (d[i] - a[i] * d_prime[i-1]) / denom
            
            # 回代
            x_sol = np.zeros(n)
            x_sol[-1] = d_prime[-1]
            for i in range(n-2, -1, -1):
                x_sol[i] = d_prime[i] - c_prime[i] * x_sol[i+1]
            
            return x_sol
        
        u_interior = thomas_solve(a, b, c, d)
        u[1:Nx] = u_interior
        
        # 更新边界值
        if boundary_left == 'neumann':
            u[0] = u[1]
        elif boundary_left == 'dirichlet':
            u[0] = 0.0  # 齐次Dirichlet
        
        if boundary_right == 'neumann':
            u[Nx] = u[Nx-1]
        elif boundary_right == 'dirichlet':
            u[Nx] = 0.0  # 齐次Dirichlet
    
    return x, u

def exact_solution(x, t, alpha, L, initial_condition='gaussian'):
    """计算解析解"""
    if initial_condition == 'gaussian':
        # 高斯脉冲的解析解（全空间近似）
        return np.exp(-x**2 / (4 * alpha * t + 1)) / np.sqrt(4 * alpha * t + 1)
    elif initial_condition == 'sine':
        # 正弦函数的解析解（有限区间[0,L]，齐次Dirichlet边界）
        return -np.exp(-alpha * (np.pi/L)**2 * t) * np.sin(np.pi * x / L)
    else:
        raise ValueError("初始条件必须是 'gaussian' 或 'sine'")

# ==================== 高斯脉冲测试 ====================
print("=== 高斯脉冲测试 ===")
alpha = 0.1
L = 10.0
T = 1.0
Nt = 1000
Nx = 200

# 数值解
x_gauss, u_gauss_num = solve_heat_equation(
    alpha=alpha, L=L, T=T, Nx=Nx, Nt=Nt,
    initial_condition='gaussian',
    boundary_left='neumann', boundary_right='neumann'
)

# 解析解
u_gauss_exact = exact_solution(x_gauss, T, alpha, L, 'gaussian')

# ==================== 正弦函数测试 ====================
print("\n=== 正弦函数测试 ===")
# 注意：对于正弦函数，通常用Dirichlet边界条件 u(0)=u(L)=0
L_sine = 2.0  # 调整区间使得有完整周期
T_sine = 0.5
Nx_sine = 100
Nt_sine = 500

x_sine, u_sine_num = solve_heat_equation(
    alpha=alpha, L=L_sine, T=T_sine, Nx=Nx_sine, Nt=Nt_sine,
    initial_condition='sine',
    boundary_left='dirichlet', boundary_right='dirichlet'
)

# 解析解
u_sine_exact = exact_solution(x_sine, T_sine, alpha, L_sine, 'sine')

# ==================== 绘图 ====================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 高斯脉冲
axes[0].plot(x_gauss, u_gauss_num, 'b-', linewidth=2, label='Numerical')
axes[0].plot(x_gauss, u_gauss_exact, 'r--', linewidth=2, label='Exact')
axes[0].set_xlabel('x')
axes[0].set_ylabel('u(x,T)')
axes[0].set_title('Gaussian Pulse Diffusion')
axes[0].legend()
axes[0].grid(True)

# 正弦函数
axes[1].plot(x_sine, u_sine_num, 'b-', linewidth=2, label='Numerical')
axes[1].plot(x_sine, u_sine_exact, 'r--', linewidth=2, label='Exact')
axes[1].set_xlabel('x')
axes[1].set_ylabel('u(x,T)')
axes[1].set_title('Sine Function Decay')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()

# ==================== 误差分析 ====================
print(f"\n高斯脉冲最大绝对误差: {np.max(np.abs(u_gauss_num - u_gauss_exact)):.6e}")
print(f"正弦函数最大绝对误差: {np.max(np.abs(u_sine_num - u_sine_exact)):.6e}")