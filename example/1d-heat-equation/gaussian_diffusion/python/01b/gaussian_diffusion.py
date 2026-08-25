import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation(alpha, L, T, Nx, Nt):
    dx = 2 * L / Nx  # 空间范围对称 [-L, L]
    dt = T / Nt
    x = np.linspace(-L, L, Nx+1)
    u = np.exp(-x**2)  # 初始高斯脉冲
    
    for _ in range(Nt):
        u_prev = u.copy()
        # Crank-Nicolson 格式的三对角系数
        a = np.ones(Nx-1) * (-alpha * dt / (2 * dx**2))
        b = np.ones(Nx-1) * (1 + alpha * dt / dx**2)
        c = np.ones(Nx-1) * (-alpha * dt / (2 * dx**2))
        d = u_prev[1:-1] + alpha * dt / (2 * dx**2) * (u_prev[2:] - 2*u_prev[1:-1] + u_prev[:-2])
        
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
    
    return x, u

# 参数设置
alpha = 0.1
L = 10.0
T = 1.0
Nt = 1000

# 不同网格密度
Nx_list = [50, 100, 200]
solutions = []
for Nx in Nx_list:
    x, u = solve_heat_equation(alpha, L, T, Nx, Nt)
    solutions.append((x, u))

# 解析解
t_exact = T
x_exact = np.linspace(-L, L, 1000)
u_exact = np.exp(-x_exact**2 / (4 * alpha * t_exact + 1)) / np.sqrt(4 * alpha * t_exact + 1)

# 绘图对比
plt.figure(figsize=(10, 6))
for i, (Nx, (x, u)) in enumerate(zip(Nx_list, solutions)):
    plt.plot(x, u, 'o-', label=f'Numerical (Nx={Nx})', markersize=4)
plt.plot(x_exact, u_exact, 'k--', linewidth=2, label='Exact Solution')
plt.xlabel('x')
plt.ylabel('u(x,T)')
plt.title('Gaussian Pulse Diffusion: Grid Density Comparison')
plt.legend()
plt.grid(True)
plt.show()