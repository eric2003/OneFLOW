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
    for i in range(len(x)):
        if x[i] >= 0.5 and x[i] <= 1.5:
            u[i] = 1.0
    return u

u_initial = initial_square_wave(x)

# minmod 限幅器
def minmod(a, b):
    if a * b <= 0:
        return 0.0
    return np.sign(a) * min(abs(a), abs(b))

# 一阶迎风（点索引形式）
def first_order_upwind_point(u_initial, dx, dt, c, nt):
    nx = len(u_initial)
    u = u_initial.copy()
    courant = c * dt / dx
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            # 处理周期性边界条件
            if i == 0:
                im1 = nx - 1  # 左邻居是最后一个点
            else:
                im1 = i - 1
            
            # 一阶迎风格式：u_i^{n+1} = u_i^n - c*Δt/dx * (u_i^n - u_{i-1}^n)
            u_new[i] = u[i] - courant * (u[i] - u[im1])
        u = u_new.copy()
    return u

# 二阶迎风 TVD (minmod) 点索引形式
def second_order_tvd_minmod_point(u_initial, dx, dt, c, nt):
    nx = len(u_initial)
    u = u_initial.copy()
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        sigma = np.zeros_like(u)
        
        # 第一步：计算minmod限幅器
        for i in range(nx):
            # 处理邻居索引（周期性边界）
            if i == 0:
                im1 = nx - 1
                ip1 = 1
            elif i == nx - 1:
                im1 = i - 1
                ip1 = 0
            else:
                im1 = i - 1
                ip1 = i + 1
            
            # 计算左右差分
            delta_left = (u[i] - u[im1]) / dx
            delta_right = (u[ip1] - u[i]) / dx
            
            # minmod限幅器
            sigma[i] = minmod(delta_left, delta_right)
        
        # 第二步：计算界面值和通量
        for i in range(nx):
            if i == 0:
                im1 = nx - 1
            else:
                im1 = i - 1
            
            # 界面值：u_{i+1/2} = u_i + 0.5*dx*σ_i
            u_interface = u[i] + 0.5 * dx * sigma[i]
            flux = c * u_interface
            
            # 通量差：F_{i+1/2} - F_{i-1/2}
            if i == 0:
                flux_im1 = c * (u[im1] + 0.5 * dx * sigma[im1])
            else:
                flux_im1 = c * (u[im1] + 0.5 * dx * sigma[im1])
            
            flux_diff = flux - flux_im1
            
            # 更新解：u_i^{n+1} = u_i^n - Δt/dx * (F_{i+1/2} - F_{i-1/2})
            u_new[i] = u[i] - (c * dt / dx) * flux_diff
        
        u = u_new.copy()
    return u

# Lax-Wendroff 中心格式（点索引形式）
def lax_wendroff_point(u_initial, dx, dt, c, nt):
    nx = len(u_initial)
    u = u_initial.copy()
    courant = c * dt / dx
    alpha = courant**2 / 2
    beta = courant / 2
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            # 处理邻居索引（周期性边界）
            if i == 0:
                im1 = nx - 1
                ip1 = 1
            elif i == nx - 1:
                im1 = i - 1
                ip1 = 0
            else:
                im1 = i - 1
                ip1 = i + 1
            
            # Lax-Wendroff公式：
            # u_i^{n+1} = u_i^n - β*(u_{i+1}^n - u_{i-1}^n) + α*(u_{i+1}^n - 2u_i^n + u_{i-1}^n)
            u_new[i] = u[i] - beta * (u[ip1] - u[im1]) + alpha * (u[ip1] - 2*u[i] + u[im1])
        u = u_new.copy()
    return u

# 计算所有解
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

u1 = first_order_upwind_point(u_initial, dx, dt, c, nt)
u2 = second_order_tvd_minmod_point(u_initial, dx, dt, c, nt)
u_central = lax_wendroff_point(u_initial, dx, dt, c, nt)

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