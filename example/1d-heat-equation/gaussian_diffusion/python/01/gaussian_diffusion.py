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
        # 使用Thomas算法求解
        a = np.ones(Nx-1) * (-alpha * dt / dx**2)
        b = np.ones(Nx-1) * (1 + 2 * alpha * dt / dx**2)
        c = np.ones(Nx-1) * (-alpha * dt / dx**2)
        d = u_prev[1:-1] + alpha * dt / dx**2 * (u_prev[2:] - 2*u_prev[1:-1] + u_prev[:-2])
        u[1:-1] = thomas_algorithm(a, b, c, d)
    
    # 理论解
    u_exact = np.exp(-x**2 / (4 * alpha * T + 1)) / np.sqrt(4 * alpha * T + 1)
    
    plt.plot(x, u, 'b-', label='Numerical')
    plt.plot(x, u_exact, 'r--', label='Exact')
    plt.legend()
    plt.show()

test_gaussian_diffusion()