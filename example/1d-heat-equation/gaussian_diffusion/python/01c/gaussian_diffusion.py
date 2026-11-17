import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation_0_to_10(alpha, L, T, Nx, Nt):
    dx = L / Nx  # 空间范围 [0, L]
    dt = T / Nt
    x = np.linspace(0, L, Nx+1)
    u = np.exp(-x**2)  # 初始高斯脉冲
    
    for _ in range(Nt):
        u_prev = u.copy()
        # Crank-Nicolson 格式的三对角系数（Neumann 边界条件在 x=0）
        a = np.ones(Nx-1) * (-alpha * dt / (2 * dx**2))
        b = np.ones(Nx-1) * (1 + alpha * dt / dx**2)
        c = np.ones(Nx-1) * (-alpha * dt / (2 * dx**2))
        d = u_prev[1:-1] + alpha * dt / (2 * dx**2) * (u_prev[2:] - 2*u_prev[1:-1] + u_prev[:-2])
        
        # 处理 x=0 的 Neumann 边界条件 (du/dx=0)
        # 离散化: (u[1] - u[-1])/(2*dx) = 0 => u[-1] = u[1]（虚拟点法）
        # 修改第一行方程： (2*b[0]*u[1] - a[0]*u[2] - c[0]*u[0]) = 2*d[0]
        # 但更简单的方法是直接修改矩阵（此处简化处理）
        d[0] += alpha * dt / (2 * dx**2) * u_prev[0]  # 近似处理
        
        # 追赶法求解
        def thomas(a, b, c, d):
            n = len(d)
            x_sol = np.zeros(n)
            for i in range(1, n):
                w = a[i] / b[i-1]
                b[i] -= w * c[i-1]
                d[i] -= w * d[i-1]
            x_sol[-1] = d[-1] / b[-1]
            for i in range(n-2, -1, -1):
                x_sol[i] = (d[i] - c[i] * x_sol[i+1]) / b[i]
            return x_sol
        
        u[1:-1] = thomas(a, b, c, d)
        # 强制满足 Neumann 条件（简化版）
        u[0] = u[1]  # 粗略近似
    
    return x, u

# 参数设置
alpha = 0.1
L = 10.0
T = 1.0
Nt = 1000

# 不同网格密度
Nx_list = [50, 100, 200]
solutions_0_to_10 = []
for Nx in Nx_list:
    x, u = solve_heat_equation_0_to_10(alpha, L, T, Nx, Nt)
    solutions_0_to_10.append((x, u))

# 解析解（对称区间 [-L, L] 的解在 [0, L] 的部分）
t_exact = T
x_exact = np.linspace(0, L, 1000)
u_exact = np.exp(-x_exact**2 / (4 * alpha * t_exact + 1)) / np.sqrt(4 * alpha * t_exact + 1)

# 绘图对比
plt.figure(figsize=(10, 6))
for i, (Nx, (x, u)) in enumerate(zip(Nx_list, solutions_0_to_10)):
    plt.plot(x, u, 'o-', label=f'Numerical (Nx={Nx})', markersize=4)
plt.plot(x_exact, u_exact, 'k--', linewidth=2, label='Exact Solution (symmetric)')
plt.xlabel('x')
plt.ylabel('u(x,T)')
plt.title('Gaussian Pulse Diffusion: 0-10 Interval Grid Density Comparison')
plt.legend()
plt.grid(True)
plt.show()