import numpy as np
import matplotlib.pyplot as plt

def simple_convection_diffusion(a=1.0, nu=0.1, t0=0.1):
    """简单的对流扩散验证"""
    xmin, xmax = -5, 10
    Nx = 200
    T = 2.0
    Nt = 1000
    
    dx = (xmax - xmin) / Nx
    dt = T / Nt
    x = np.linspace(xmin, xmax, Nx+1)
    
    # 初始条件
    u = np.exp(-x**2 / (4*nu*t0)) / np.sqrt(4*np.pi*nu*t0)
    
    # 简单的时间步进（迎风格式）
    for n in range(Nt):
        u_old = u.copy()
        for i in range(1, Nx):
            # 对流项：迎风
            if a > 0:
                conv = a * (u_old[i] - u_old[i-1]) / dx
            else:
                conv = a * (u_old[i+1] - u_old[i]) / dx
            
            # 扩散项：中心差分
            diff = nu * (u_old[i+1] - 2*u_old[i] + u_old[i-1]) / (dx**2)
            
            u[i] = u_old[i] - dt * conv + dt * diff
        
        # 边界
        u[0] = u[1]
        u[Nx] = u[Nx-1]
    
    # 解析解
    u_exact = np.exp(-(x - a*T)**2 / (4*nu*(T+t0))) / np.sqrt(4*np.pi*nu*(T+t0))
    
    # 绘图
    plt.figure(figsize=(10, 6))
    plt.plot(x, u, 'b-', label='Numerical')
    plt.plot(x, u_exact, 'r--', label='Exact')
    plt.xlabel('x')
    plt.ylabel('u(x,T)')
    plt.title(f'Convection-Diffusion: a={a}, ν={nu}, T={T}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
    
    # 分析
    idx_num = np.argmax(u)
    idx_exact = np.argmax(u_exact)
    print(f"Numerical peak at x = {x[idx_num]:.4f}")
    print(f"Exact peak at x = {x[idx_exact]:.4f}")
    print(f"Theoretical: x = aT = {a*T:.4f}")

if __name__ == "__main__":
    simple_convection_diffusion()