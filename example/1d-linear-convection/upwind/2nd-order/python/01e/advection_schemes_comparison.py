import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

# 参数设置（统一）
nx = 101
L = 4.0
x = np.linspace(0, L, nx)
dx = x[1] - x[0]
c = 1.0
dt = 0.01
total_time = 1.0
nt = int(total_time / dt + 1e-8)

# 初始条件
def initial_square_wave(x):
    u = np.zeros_like(x)
    u[(x >= 0.5) & (x <= 1.5)] = 1.0
    return u

u_initial = initial_square_wave(x)

# minmod 限幅器
def minmod(a, b):
    return 0.5 * (np.sign(a) + np.sign(b)) * np.minimum(np.abs(a), np.abs(b))

# 一阶迎风
def first_order_upwind(u_initial, dx, dt, c, nt):
    u = u_initial.copy()
    courant = c * dt / dx
    for n in range(nt):
        u = u - courant * (u - np.roll(u, 1))
    return u

# 二阶迎风 TVD (minmod)
def second_order_tvd_minmod(u_initial, dx, dt, c, nt):
    u = u_initial.copy()
    for n in range(nt):
        delta_left  = (u - np.roll(u, 1)) / dx
        delta_right = (np.roll(u, -1) - u) / dx
        sigma = minmod(delta_left, delta_right)
        u_interface = u + 0.5 * dx * sigma
        flux = c * u_interface
        u = u - (c * dt / dx) * (flux - np.roll(flux, 1))
    return u

# Lax-Wendroff 中心格式
def lax_wendroff(u_initial, dx, dt, c, nt):
    u = u_initial.copy()
    courant = c * dt / dx
    alpha = courant**2 / 2
    beta = courant / 2
    for n in range(nt):
        u_ip1 = np.roll(u, -1)
        u_im1 = np.roll(u, 1)
        u = u - beta * (u_ip1 - u_im1) + alpha * (u_ip1 - 2*u + u_im1)
    return u

# 计算所有解
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

u1 = first_order_upwind(u_initial, dx, dt, c, nt)
u2 = second_order_tvd_minmod(u_initial, dx, dt, c, nt)
u_central = lax_wendroff(u_initial, dx, dt, c, nt)

# 误差输出函数
def print_errors(label, u_num, u_exact, dx):
    error_linf = np.max(np.abs(u_num - u_exact))
    error_l1   = np.sum(np.abs(u_num - u_exact)) * dx
    error_l2   = np.sqrt(np.sum((u_num - u_exact)**2) * dx)
    print(f"=== {label} 误差统计 (t={total_time}) ===")
    print(f"L∞ 误差: {error_linf:.6f}")
    print(f"L1 误差: {error_l1:.6f}")
    print(f"L2 误差: {error_l2:.6f}")
    print(f"最大值: {np.max(u_num):.6f}, 最小值: {np.min(u_num):.6f}\n")

print_errors("一阶迎风数值解", u1, u_exact, dx)
print_errors("二阶迎风（minmod TVD）数值解", u2, u_exact, dx)
print_errors("Lax-Wendroff 中心格式数值解", u_central, u_exact, dx)

# 绘图对比
plt.figure(figsize=(12, 8))
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解 (精确平移)')
plt.plot(x, u1, 'bo-', lw=1, markersize=6, markevery=2, label='一阶迎风')
plt.plot(x, u2, 'rs-', lw=1, markersize=6, markevery=2, label='二阶迎风（minmod TVD）')
plt.plot(x, u_central, 'g^-', lw=1, markersize=6, markevery=2, label='Lax-Wendroff 中心格式')
plt.xlabel('x')
plt.ylabel('u')
plt.title('一维对流方程：中心格式 vs 1阶 vs 2阶迎风格式 方波传播对比 (nx=101, t=1.0)')
plt.legend(fontsize=11)
plt.grid(True)
plt.ylim(-0.3, 1.4)
plt.show()