import numpy as np
import matplotlib.pyplot as plt

# 参数设置
a = 1.0          # 对流速度
nx = 100         # 空间网格点数
dx = 1.0 / nx    # 空间步长
dt = 0.005       # 时间步长（满足 CFL 条件）
t_max = 0.2      # 总模拟时间
nt = int(t_max / dt)  # 时间步数

# 初始化网格和初值
x = np.linspace(0, 1, nx+1)  # 空间网格（包括边界）
u = np.zeros(nx+1)           # 数值解
u_exact = np.zeros(nx+1)     # 理论解

# 设置初始条件（方波）
for i in range(nx+1):
    if 0.2 <= x[i] <= 0.4:
        u[i] = 1.0
    else:
        u[i] = 0.0

# 复制初始条件作为理论解的参考
u_initial = u.copy()

# Rusanov flux 函数
def rusanov_flux(u_L, u_R, a):
    F_L = a * u_L  # 左状态通量
    F_R = a * u_R  # 右状态通量
    alpha = abs(a) # 最大波速
    flux = 0.5 * (F_L + F_R) - 0.5 * alpha * (u_R - u_L)
    return flux

# 时间推进
u_new = np.zeros(nx+1)
for n in range(nt):
    # 计算每个单元边界的通量
    #for i in range(nx):
    for i in range(1, nx-1):
        f_left = rusanov_flux(u[i-1], u[i], a)
        f_right = rusanov_flux(u[i], u[i+1], a)
        # 更新 u
        u_new[i] = u[i] - (dt / dx) * (f_right - f_left)        
        
    # 更新边界条件（周期性边界）
    u_new[0] = u_new[nx]
    u_new[nx] = u_new[0]
    
    # 更新 u
    u = u_new.copy()

# 计算理论解（初始条件平移）
for i in range(nx+1):
    x_shifted = (x[i] - a * t_max) % 1.0  # 周期性平移
    u_exact[i] = u_initial[int(x_shifted / dx)]

# 绘图
plt.figure(figsize=(10, 6))
plt.plot(x, u, 'b-', label='Numerical Solution (Rusanov)')
plt.plot(x, u_exact, 'r--', label='Exact Solution')
plt.xlabel('x')
plt.ylabel('u')
plt.title(f'1D Linear Convection Equation at t = {t_max}')
plt.legend()
plt.grid(True)
plt.show()

# 计算 L2 误差
error = np.sqrt(np.sum((u - u_exact)**2) * dx)
print(f"L2 Error: {error:.6f}")