import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation(alpha, L, T, Nx, Nt):
    dx = L / Nx
    dt = T / Nt
    x = np.linspace(0, L, Nx+1)
    u = np.exp(-x**2)  # 初始高斯脉冲
    
    for _ in range(Nt):
        u_prev = u.copy()
        # 构造三对角矩阵的系数（Crank-Nicolson）
        a = np.ones(Nx-1) * (-alpha * dt / (2 * dx**2))
        b = np.ones(Nx-1) * (1 + alpha * dt / dx**2)
        c = np.ones(Nx-1) * (-alpha * dt / (2 * dx**2))
        d = u_prev[1:-1] + alpha * dt / (2 * dx**2) * (u_prev[2:] - 2*u_prev[1:-1] + u_prev[:-2])
        
        # 追赶法求解
        def thomas(a, b, c, d):
            n = len(d)
            x = np.zeros(n)
            for i in range(1, n):
                w = a[i] / b[i-1]
                b[i] -= w * c[i-1]
                d[i] -= w * d[i-1]
            x[-1] = d[-1] / b[-1]
            for i in range(n-2, -1, -1):
                x[i] = (d[i] - c[i] * x[i+1]) / b[i]
            return x
        
        u[1:-1] = thomas(a, b, c, d)
    
    return x, u

# 参数设置
alpha = 0.1  # 热扩散系数
L = 10.0     # 空间范围
T = 1.0      # 总时间
Nx = 100     # 空间网格数
Nt = 1000    # 时间步数

x, u = solve_heat_equation(alpha, L, T, Nx, Nt)

# 解析解
t_exact = T
u_exact = np.exp(-x**2 / (4 * alpha * t_exact + 1)) / np.sqrt(4 * alpha * t_exact + 1)

# 绘图对比
plt.plot(x, u, 'b-', label='Numerical')
plt.plot(x, u_exact, 'r--', label='Exact')
plt.xlabel('x')
plt.ylabel('u(x,T)')
plt.title('Gaussian Pulse Diffusion (1D Heat Equation)')
plt.legend()
plt.show()