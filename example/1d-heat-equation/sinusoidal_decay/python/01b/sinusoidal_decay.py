import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation_general(alpha, xmin, xmax, T, Nx, Nt, 
                               initial_condition='gaussian',
                               boundary_left='neumann', boundary_right='neumann'):
    """
    求解一维热传导方程（通用区间版本）
    
    参数:
    alpha: 扩散系数
    xmin: 空间域左边界
    xmax: 空间域右边界
    T: 总时间
    Nx: 空间网格点数
    Nt: 时间步数
    initial_condition: 'gaussian' 或 'sine' 或自定义函数
    boundary_left: 左边界条件 'neumann' 或 'dirichlet'
    boundary_right: 右边界条件 'neumann' 或 'dirichlet'
    """
    L = xmax - xmin  # 区间长度
    dx = L / Nx
    dt = T / Nt
    x = np.linspace(xmin, xmax, Nx+1)
    r = alpha * dt / (dx**2)
    
    # 设置初始条件
    if initial_condition == 'gaussian':
        # 高斯脉冲，中心在x=0附近
        u = np.exp(-x**2)
    elif initial_condition == 'sine':
        # 正弦函数，周期与区间长度匹配
        period = L  # 完整周期长度为L
        u = -np.sin(2*np.pi * (x) / period)
    elif callable(initial_condition):
        # 自定义初始条件函数
        u = initial_condition(x)
    else:
        raise ValueError("初始条件必须是 'gaussian' 或 'sine' 或函数")
    
    # 三对角矩阵系数（内部点）
    n_interior = Nx - 1  # 内部点数
    a = np.full(n_interior, -r/2.0)  # 下对角线
    b = np.full(n_interior, 1.0 + r) # 主对角线  
    c = np.full(n_interior, -r/2.0)  # 上对角线
    
    # 时间推进
    for _ in range(Nt):
        un = u.copy()
        d = np.zeros(n_interior)  # 右端项
        
        # 构建右端项（CN格式的显式部分）
        for i in range(1, Nx):
            j = i - 1  # 内部点索引
            
            # 基本项（标准三点格式）
            if 1 < i < Nx-1:
                d[j] = (r/2.0) * un[i-1] + (1.0 - r) * un[i] + (r/2.0) * un[i+1]
            
            # 左边界邻近点 (i=1)
            elif i == 1:
                term = (1.0 - r) * un[i] + (r/2.0) * un[i+1]
                if boundary_left == 'neumann':
                    # Neumann边界：u[xmin] ≈ u[xmin+dx]
                    term += (r/2.0) * un[i-1]
                elif boundary_left == 'dirichlet':
                    # Dirichlet边界：u[xmin]固定，已包含在un[i-1]中
                    term += (r/2.0) * un[i-1]
                d[j] = term
            
            # 右边界邻近点 (i=Nx-1)
            elif i == Nx-1:
                term = (r/2.0) * un[i-1] + (1.0 - r) * un[i]
                if boundary_right == 'neumann':
                    # Neumann边界：u[xmax] ≈ u[xmax-dx]
                    term += (r/2.0) * un[i+1]
                elif boundary_right == 'dirichlet':
                    # Dirichlet边界：u[xmax]固定，已包含在un[i+1]中
                    term += (r/2.0) * un[i+1]
                d[j] = term
        
        # 应用边界条件到右端项
        if boundary_left == 'neumann':
            # 隐式Neumann边界：u_0 = u_1
            # 对应第一个方程的修正
            d[0] += (r/2.0) * un[0]
        # Dirichlet边界不需要额外修正，因为边界值在un中已经固定
        
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
        
        # 求解内部点
        u_interior = thomas_solve(a, b, c, d)
        u[1:Nx] = u_interior
        
        # 更新边界值
        if boundary_left == 'neumann':
            u[0] = u[1]  # Neumann边界条件
        elif boundary_left == 'dirichlet':
            u[0] = 0.0  # 齐次Dirichlet边界，可根据需要修改
        
        if boundary_right == 'neumann':
            u[Nx] = u[Nx-1]  # Neumann边界条件
        elif boundary_right == 'dirichlet':
            u[Nx] = 0.0  # 齐次Dirichlet边界，可根据需要修改
    
    return x, u

def exact_solution_general(x, t, alpha, xmin, xmax, initial_condition='gaussian'):
    """计算解析解（通用区间版本）"""
    L = xmax - xmin
    
    if initial_condition == 'gaussian':
        # 高斯脉冲的解析解（全空间近似）
        return np.exp(-x**2 / (4 * alpha * t + 1)) / np.sqrt(4 * alpha * t + 1)
    elif initial_condition == 'sine':
        # 正弦函数的解析解（有限区间[xmin, xmax]，齐次Dirichlet边界）
        # u(x,0) = -sin(2π(x)/L), 在x=xmin和x=xmax处为0
        return -np.exp(-alpha * (2*np.pi/L)**2 * t) * np.sin(2*np.pi * (x) / L)
    else:
        raise ValueError("初始条件必须是 'gaussian' 或 'sine'")

# ==================== 示例1：高斯脉冲（对称区间） ====================
print("=== 示例1：高斯脉冲（对称区间） ===")
alpha = 0.1
#xmin = -5.0  # 从负值开始
xmin = 0.0  # 从负值开始
xmax = 5.0   # 到正值结束
T = 1.0
Nx = 200
Nt = 1000

# 数值解
x_gauss, u_gauss_num = solve_heat_equation_general(
    alpha=alpha, xmin=xmin, xmax=xmax, T=T, Nx=Nx, Nt=Nt,
    initial_condition='gaussian',
    boundary_left='neumann', boundary_right='neumann'
)

# 解析解
u_gauss_exact = exact_solution_general(x_gauss, T, alpha, xmin, xmax, 'gaussian')

# ==================== 示例2：正弦函数（任意区间） ====================
print("\n=== 示例2：正弦函数（任意区间） ===")
xmin_sine = -1.0
xmax_sine = 1.0  # 区间长度为2.0
T_sine = 0.5
Nx_sine = 50
Nt_sine = 500

# 数值解
x_sine, u_sine_num = solve_heat_equation_general(
    alpha=alpha, xmin=xmin_sine, xmax=xmax_sine, T=T_sine, 
    Nx=Nx_sine, Nt=Nt_sine,
    initial_condition='sine',
    boundary_left='dirichlet', boundary_right='dirichlet'
)

# 解析解
u_sine_exact = exact_solution_general(x_sine, T_sine, alpha, 
                                     xmin_sine, xmax_sine, 'sine')

# ==================== 绘图 ====================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 高斯脉冲（对称区间）
axes[0].plot(x_gauss, u_gauss_num, 'b-', linewidth=2, label='Numerical')
axes[0].plot(x_gauss, u_gauss_exact, 'r--', linewidth=2, label='Exact')
axes[0].set_xlabel('x')
axes[0].set_ylabel('u(x,T)')
axes[0].set_title(f'Gaussian Pulse (x∈[{xmin},{xmax}])')
axes[0].legend()
axes[0].grid(True)

# 正弦函数（任意区间）
axes[1].plot(x_sine, u_sine_num, 'b-', linewidth=2, label='Numerical')
axes[1].plot(x_sine, u_sine_exact, 'r--', linewidth=2, label='Exact')
axes[1].set_xlabel('x')
axes[1].set_ylabel('u(x,T)')
axes[1].set_title(f'Sine Function (x∈[{xmin_sine},{xmax_sine}])')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.show()

# ==================== 误差分析 ====================
print(f"\n高斯脉冲最大绝对误差: {np.max(np.abs(u_gauss_num - u_gauss_exact)):.6e}")
print(f"正弦函数最大绝对误差: {np.max(np.abs(u_sine_num - u_sine_exact)):.6e}")

# ==================== 验证区间参数 ====================
print(f"\n验证区间设置:")
print(f"高斯脉冲: x范围 = [{x_gauss[0]:.2f}, {x_gauss[-1]:.2f}], 点数 = {len(x_gauss)}")
print(f"正弦函数: x范围 = [{x_sine[0]:.2f}, {x_sine[-1]:.2f}], 点数 = {len(x_sine)}")