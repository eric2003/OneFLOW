import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# ========== 参数设置 ==========
nx = 101
L = 4.0
x = np.linspace(0, L, nx)
dx = x[1] - x[0]
c = 1.0
dt = 0.01
total_time = 1.0
nt = int(total_time / dt + 1e-8)

# ========== 初始条件 ==========
def initial_square_wave(x):
    u = np.zeros_like(x)
    for i in range(len(x)):
        if x[i] >= 0.5 and x[i] <= 1.5:
            u[i] = 1.0
    return u

u_initial = initial_square_wave(x)

# ========== minmod 限幅器 ==========
def minmod(a, b):
    if a * b <= 0:
        return 0.0
    return np.sign(a) * min(abs(a), abs(b))

# ========== 一阶迎风 ==========
def first_order_upwind_point(u_initial, dx, dt, c, nt):
    nx = len(u_initial)
    u = u_initial.copy()
    courant = c * dt / dx
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            im1 = nx - 1 if i == 0 else i - 1
            u_new[i] = u[i] - courant * (u[i] - u[im1])
        u = u_new.copy()
    return u

# ========== Lax-Wendroff ==========
def lax_wendroff_point(u_initial, dx, dt, c, nt):
    nx = len(u_initial)
    u = u_initial.copy()
    courant = c * dt / dx
    alpha = courant**2 / 2
    beta = courant / 2
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            im1 = nx - 1 if i == 0 else i - 1
            ip1 = 1 if i == nx - 1 else i + 1
            u_new[i] = u[i] - beta * (u[ip1] - u[im1]) + alpha * (u[ip1] - 2*u[i] + u[im1])
        u = u_new.copy()
    return u

# ========== Rusanov 格式 ==========
def rusanov_scalar_point(u_initial, dx, dt, c, nt):
    nx = len(u_initial)
    u = u_initial.copy()
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            im1 = nx - 1 if i == 0 else i - 1
            ip1 = 1 if i == nx - 1 else i + 1
            
            u_L = u[i]
            u_R = u[ip1]
            f_L = c * u_L
            f_R = c * u_R
            lambda_max = abs(c) # 线性对流
            
            flux_ip12 = 0.5 * (f_L + f_R) - 0.5 * lambda_max * (u_R - u_L)
            
            u_L_im1 = u[im1]
            u_R_im1 = u[i]
            f_L_im1 = c * u_L_im1
            f_R_im1 = c * u_R_im1
            flux_im12 = 0.5 * (f_L_im1 + f_R_im1) - 0.5 * lambda_max * (u_R_im1 - u_L_im1)
            
            flux_diff = flux_ip12 - flux_im12
            u_new[i] = u[i] - (c * dt / dx) * flux_diff
        u = u_new.copy()
    return u

# ========== 计算所有解 ==========
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

u1_upwind = first_order_upwind_point(u_initial, dx, dt, c, nt)
u_central = lax_wendroff_point(u_initial, dx, dt, c, nt)
u_rusanov = rusanov_scalar_point(u_initial, dx, dt, c, nt)

# ========== 误差输出函数 ==========
def print_errors(label, u_num, u_exact, dx):
    error_linf = np.max(np.abs(u_num - u_exact))
    error_l1   = np.sum(np.abs(u_num - u_exact)) * dx
    error_l2   = np.sqrt(np.sum((u_num - u_exact)**2) * dx)
    print(f"=== {label} 误差统计 (t={total_time}) ===")
    print(f"L∞ 误差: {error_linf:.6f}")
    print(f"L1 误差: {error_l1:.6f}")
    print(f"L2 误差: {error_l2:.6f}\n")

print_errors("一阶迎风数值解", u1_upwind, u_exact, dx)
print_errors("Lax-Wendroff 中心格式数值解", u_central, u_exact, dx)
print_errors("Rusanov 格式数值解", u_rusanov, u_exact, dx)

# ========== 绘图对比 ==========
plt.figure(figsize=(12, 10))

# 子图1：所有方法对比
plt.subplot(2, 2, 1)
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解')
plt.plot(x, u1_upwind, 'bo-', lw=1, markersize=4, label='一阶迎风')
plt.plot(x, u_central, 'g^-', lw=1, markersize=4, label='Lax-Wendroff')
plt.plot(x, u_rusanov, 'rs-', lw=1, markersize=4, label='Rusanov格式')
plt.xlabel('x')
plt.ylabel('u')
plt.title('线性对流方程：各方法对比')
plt.legend()
plt.grid(True)

# 子图2：Rusanov与一阶迎风对比（应该非常接近）
plt.subplot(2, 2, 2)
plt.plot(x, u_exact, 'k-', lw=2, label='理论解', alpha=0.5)
plt.plot(x, u1_upwind, 'bo-', lw=1.5, markersize=4, label='一阶迎风', alpha=0.6)
plt.plot(x, u_rusanov, 'rs-', lw=1.5, markersize=4, label='Rusanov格式', alpha=0.8)
plt.xlabel('x')
plt.ylabel('u')
plt.title('Rusanov vs 一阶迎风（线性对流）')
plt.legend()
plt.grid(True)

# 子图3：误差对比（对数坐标）
plt.subplot(2, 2, 3)
error_rusanov = np.abs(u_rusanov - u_exact)
error_upwind = np.abs(u1_upwind - u_exact)
error_central = np.abs(u_central - u_exact)
plt.plot(x, error_upwind, 'b-', lw=1.5, label='一阶迎风', alpha=0.7)
plt.plot(x, error_central, 'g-', lw=1.5, label='Lax-Wendroff', alpha=0.7)
plt.plot(x, error_rusanov, 'r-', lw=2, label='Rusanov格式', alpha=0.8)
plt.xlabel('x')
plt.ylabel('绝对误差')
plt.title('误差对比（Rusanov与迎风几乎重合）')
plt.legend()
plt.grid(True)
plt.yscale('log')

# 子图4：Rusanov的耗散特性（与精确解对比）
plt.subplot(2, 2, 4)
plt.plot(x, u_exact, 'k-', lw=2, label='理论解（方波）')
plt.plot(x, u_rusanov, 'r--', lw=2, label='Rusanov格式（被平滑）')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Rusanov格式的数值耗散效应')
plt.legend()
plt.grid(True)
plt.ylim(-0.1, 1.2)

plt.tight_layout()
plt.savefig('rusanov_comparison.png', dpi=300, bbox_inches='tight')
plt.show()