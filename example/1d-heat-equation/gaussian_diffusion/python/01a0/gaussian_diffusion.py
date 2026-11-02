import numpy as np
import matplotlib.pyplot as plt

def thomas_algorithm(a, b, c, d):
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

def test_gaussian_diffusion():
    L = 10.0
    T = 1.0
    Nx = 100
    Nt = 1000
    alpha = 0.1
    dx = L / Nx
    dt = T / Nt
    
    x = np.linspace(0, L, Nx+1)
    u = np.exp(-x**2)  # 初始条件
    
    for _ in range(Nt):
        u_prev = u.copy()
        
        # 处理 Neumann 边界条件 (du/dx=0 at x=0)
        # 使用虚拟点法: u[-1] = u[1] (对于左边界)
        # 这里我们通过修改方程来体现
        
        a = np.ones(Nx-1) * (-alpha * dt / dx**2)
        b = np.ones(Nx-1) * (1 + 2 * alpha * dt / dx**2)
        c = np.ones(Nx-1) * (-alpha * dt / dx**2)
        
        d = u_prev[1:-1].copy()
        d += alpha * dt / (2 * dx**2) * (u_prev[2:] - 2*u_prev[1:-1] + u_prev[:-2])
        
        # 修改第一行和最后一行以体现 Neumann 条件
        # 对于左边界 (i=0): 2*b[0]*u[1] - 2*a[0]*u[2]? (更精确的实现需要调整)
        # 简化处理: 使用前向差分近似 du/dx=0 => u[0] ≈ u[1]
        # 所以我们可以将 u[0] 从方程中消除
        # 这里我们采用更简单的方法: 在计算 d 时考虑边界影响
        
        # 更精确的边界处理需要修改矩阵的第一个和最后一个方程
        # 但为了简单,我们可以这样做:
        d[0] += alpha * dt / (2 * dx**2) * u_prev[0]  # 近似处理左边界
        # 类似地处理右边界
        d[-1] += alpha * dt / (2 * dx**2) * u_prev[-1]
        
        u_interior = thomas_algorithm(a, b, c, d)
        u[1:-1] = u_interior
        
        # 强制满足 Neumann 条件 (简化版)
        u[0] = u[1]  # 粗略近似
        u[-1] = u[-2]  # 粗略近似
    
    # 理论解
    u_exact = np.exp(-x**2 / (4 * alpha * T + 1)) / np.sqrt(4 * alpha * T + 1)
    
    plt.plot(x, u, 'b-', label='Numerical')
    plt.plot(x, u_exact, 'r--', label='Exact')
    plt.legend()
    plt.show()

test_gaussian_diffusion()