import numpy as np
import matplotlib.pyplot as plt

# Rusanov flux 函数
def rusanov_flux(u_L, u_R, a):
    F_L = a * u_L  # 左状态通量
    F_R = a * u_R  # 右状态通量
    alpha = abs(a) # 最大波速
    flux = 0.5 * (F_L + F_R) - 0.5 * alpha * (u_R - u_L)
    return flux
    
# 初始条件
def initial_condition(x):
    u0 = np.zeros_like(x)
    for i in range(len(x)):
        if 0.5 <= x[i] <= 1.0:
            u0[i] = 2.0
        else:
            u0[i] = 1.0
    return u0

# 理论解
def analytical_solution(x, t, a):
    # 初始条件沿 x - at 平移
    x_shifted = x - a * t
    return initial_condition((x_shifted + L) % L)  # 周期边界条件    

# rusanov 方法
def rusanov(u0, nx, nt, dx, dt, c):
    un = u.copy()
    u_new = np.zeros_like(u)
    
    uL = np.zeros_like(u)
    uR = np.zeros_like(u)
    
    for n in range(nt):
        for i in range(0, nx-1):
            uL[i] = un[i]
            uR[i] = un[i+1]
        for i in range(1, nx-1):
            f_left = rusanov_flux(uL[i-1], uR[i-1], c)
            f_right = rusanov_flux(uL[i], uR[i], c)
            # 更新 u
            u[i] = un[i] - (dt / dx) * (f_right - f_left)
        # 边界条件（周期性）
        u[0] = u[-2]
        u[-1] = u[1]
        
        # 更新时间步
        un = u.copy()
    
    return u            

nxc = 40
nx = nxc + 1
L  = 2
dx = L / (nxc)
nt = 25    #nt is the number of timesteps we want to calculate
dt = .025  #dt is the amount of time each timestep covers (delta t)
c = 1      #assume wavespeed of c = 1

T = dt * nt
x = np.linspace(0, 2, nxc + 1)

# 计算初始条件
u = initial_condition(x)

# 计算数值解
u_numerical = rusanov(u, nx, nt, dx, dt, c)

# 计算理论解
u_analytical = analytical_solution(x, T, c)
       
# 可视化
plt.figure(figsize=(10, 6))
plt.scatter(x, u_numerical, facecolor="none", edgecolor="blue", s=20, linewidths=0.5, label='Numerical (Rusanov)')
plt.plot(x, u_analytical, 'r--', label='Analytical')
plt.xlabel('x')
plt.ylabel('u')
plt.title(f'1D Convection Equation at t = {T:.2f}')
plt.legend()
plt.tight_layout()
#plt.grid(True)
plt.show()