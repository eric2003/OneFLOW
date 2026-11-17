import numpy as np
import matplotlib.pyplot as plt

# 参数设置
a = 1.0          # 对流速度
L = 1.0          # 计算域长度
T = 0.4          # 总时间
nx = 101         # 空间网格数
dx = L / (nx - 1) # 空间步长
cfl = 0.8        # CFL 数
dt = cfl * dx / a # 时间步长
nt = int(T / dt) + 1 # 时间步数

# 空间和时间网格
x = np.linspace(0, L, nx)
t = np.linspace(0, T, nt)

# 初始条件
def initial_condition(x):
    u0 = np.zeros_like(x)
    for i in range(len(x)):
        if 0.1 <= x[i] <= 0.3:
            u0[i] = 1.0
    return u0

# 理论解
def analytical_solution(x, t, a):
    # 初始条件沿 x - at 平移
    x_shifted = x - a * t
    return initial_condition((x_shifted + L) % L)  # 周期边界条件

# Lax-Friedrichs 方法
def lax_friedrichs(u0, nx, nt, dx, dt, a):
    u = u0.copy()
    u_new = np.zeros_like(u)
    
    for n in range(nt):
        for i in range(1, nx-1):
            # Lax-Friedrichs 通量
            f_left = 0.5 * (a * u[i] + a * u[i-1]) - 0.5 * (dx / dt) * (u[i] - u[i-1])
            f_right = 0.5 * (a * u[i+1] + a * u[i]) - 0.5 * (dx / dt) * (u[i+1] - u[i])
            # 更新 u
            u_new[i] = u[i] - (dt / dx) * (f_right - f_left)
        
        # 边界条件（周期性）
        u_new[0] = u_new[-2]
        u_new[-1] = u_new[1]
        
        # 更新时间步
        u = u_new.copy()
    
    return u

# 计算初始条件
u0 = initial_condition(x)

# 计算数值解
u_numerical = lax_friedrichs(u0, nx, nt, dx, dt, a)

# 计算理论解
u_analytical = analytical_solution(x, T, a)

# 可视化
plt.figure(figsize=(10, 6))
plt.plot(x, u_numerical, 'b-', label='Numerical (Lax-Friedrichs)')
plt.plot(x, u_analytical, 'r--', label='Analytical')
plt.xlabel('x')
plt.ylabel('u')
plt.title(f'1D Convection Equation at t = {T:.2f}')
plt.legend()
plt.grid(True)
plt.show()

# 计算误差
error = np.abs(u_numerical - u_analytical)
print(f"最大误差: {np.max(error):.6f}")
print(f"L2 范数误差: {np.sqrt(np.sum(error**2) * dx):.6f}")