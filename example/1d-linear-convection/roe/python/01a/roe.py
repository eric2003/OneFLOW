import numpy as np
import matplotlib.pyplot as plt

# 解决中文显示问题
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False  # 正常显示负号

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

# ========== 一阶迎风格式（点索引） ==========
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

# ========== 二阶迎风 TVD (minmod) 点索引 ==========
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

# ========== Lax-Wendroff 中心格式（点索引） ==========
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

# ========== 标量Roe格式（修正版） ==========
def roe_scalar_point(u_initial, dx, dt, c, nt):
    """
    标量对流方程的Roe格式（修正版）
    对于线性对流 f(u) = c*u，应该退化为一阶迎风
    """
    nx = len(u_initial)
    u = u_initial.copy()
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        
        for i in range(nx):
            # 处理周期性边界
            if i == 0:
                im1 = nx - 1  # i-1
                ip1 = 1       # i+1
            elif i == nx - 1:
                im1 = i - 1
                ip1 = 0
            else:
                im1 = i - 1
                ip1 = i + 1
            
            # ========== 关键修正 ==========
            # Roe格式计算界面 i+1/2 的通量
            # 使用 u_L = u[i] 和 u_R = u[ip1]
            u_L = u[i]
            u_R = u[ip1]
            
            # 通量函数 f(u) = c*u
            f_L = c * u_L
            f_R = c * u_R
            
            # Roe平均速度（线性对流情况下为常数c）
            a_roe = c
            
            # Roe通量公式
            F_ip12 = 0.5 * (f_L + f_R) - 0.5 * np.abs(a_roe) * (u_R - u_L)
            
            # 计算界面 i-1/2 的通量
            # 使用 u_L = u[im1] 和 u_R = u[i]
            u_L_im1 = u[im1]
            u_R_im1 = u[i]
            
            f_L_im1 = c * u_L_im1
            f_R_im1 = c * u_R_im1
            F_im12 = 0.5 * (f_L_im1 + f_R_im1) - 0.5 * np.abs(a_roe) * (u_R_im1 - u_L_im1)
            
            # 通量差
            flux_diff = F_ip12 - F_im12
            
            # 更新解
            u_new[i] = u[i] - (c * dt / dx) * flux_diff
        
        u = u_new.copy()
    
    return u

# ========== 非线性Burgers方程的Roe格式 ==========
def roe_burgers_point(u_initial, dx, dt, nt):
    """
    Burgers方程: u_t + (u^2/2)_x = 0
    使用Roe格式
    """
    nx = len(u_initial)
    u = u_initial.copy()
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        
        for i in range(nx):
            # 边界处理
            if i == 0:
                im1 = nx - 1
                ip1 = 1
            elif i == nx - 1:
                im1 = i - 1
                ip1 = 0
            else:
                im1 = i - 1
                ip1 = i + 1
            
            # --- Burgers方程的Roe格式 ---
            # 界面 i+1/2
            u_L = u[i]
            u_R = u[ip1]
            
            # 通量函数 f(u) = u^2/2
            f_L = 0.5 * u_L**2
            f_R = 0.5 * u_R**2
            
            # Roe平均速度 a = (u_L + u_R) / 2
            a_roe = 0.5 * (u_L + u_R)
            
            # Roe通量
            F_ip12 = 0.5 * (f_L + f_R) - 0.5 * np.abs(a_roe) * (u_R - u_L)
            
            # 界面 i-1/2
            u_L_im1 = u[im1]
            u_R_im1 = u[i]
            
            f_L_im1 = 0.5 * u_L_im1**2
            f_R_im1 = 0.5 * u_R_im1**2
            F_im12 = 0.5 * (f_L_im1 + f_R_im1) - 0.5 * np.abs(0.5 * (u_L_im1 + u_R_im1)) * (u_R_im1 - u_L_im1)
            
            # 通量差
            flux_diff = F_ip12 - F_im12
            
            u_new[i] = u[i] - (dt / dx) * flux_diff
        
        u = u_new.copy()
    
    return u

# ========== 计算所有解 ==========
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

u1_upwind = first_order_upwind_point(u_initial, dx, dt, c, nt)
u2_tvd = second_order_tvd_minmod_point(u_initial, dx, dt, c, nt)
u_central = lax_wendroff_point(u_initial, dx, dt, c, nt)
u_roe_linear = roe_scalar_point(u_initial, dx, dt, c, nt)
u_roe_burgers = roe_burgers_point(u_initial, dx, dt, nt)

# ========== 误差输出函数 ==========
def print_errors(label, u_num, u_exact, dx):
    error_linf = np.max(np.abs(u_num - u_exact))
    error_l1   = np.sum(np.abs(u_num - u_exact)) * dx
    error_l2   = np.sqrt(np.sum((u_num - u_exact)**2) * dx)
    print(f"=== {label} 误差统计 (t={total_time}) ===")
    print(f"L∞ 误差: {error_linf:.6f}")
    print(f"L1 误差: {error_l1:.6f}")
    print(f"L2 误差: {error_l2:.6f}")
    print(f"最大值: {np.max(u_num):.6f}, 最小值: {np.min(u_num):.6f}\n")

print_errors("一阶迎风数值解", u1_upwind, u_exact, dx)
print_errors("二阶迎风（minmod TVD）数值解", u2_tvd, u_exact, dx)
print_errors("Lax-Wendroff 中心格式数值解", u_central, u_exact, dx)
print_errors("Roe格式（线性对流）- 修正版", u_roe_linear, u_exact, dx)

print("\n=== 特别说明 ===")
print(f"Roe与一阶迎风的最大差异: {np.max(np.abs(u_roe_linear - u1_upwind)):.10f}")
print(f"Roe与一阶迎风的L1差异: {np.sum(np.abs(u_roe_linear - u1_upwind)) * dx:.10f}")

# ========== 绘图对比 ==========
plt.figure(figsize=(14, 10))

# 子图1：线性对流 - 所有方法对比
plt.subplot(2, 2, 1)
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解')
plt.plot(x, u1_upwind, 'bo-', lw=1, markersize=4, label='一阶迎风', alpha=0.7)
plt.plot(x, u2_tvd, 'rs-', lw=1, markersize=4, label='二阶TVD', alpha=0.7)
plt.plot(x, u_central, 'g^-', lw=1, markersize=4, label='Lax-Wendroff', alpha=0.7)
plt.plot(x, u_roe_linear, 'm*-', lw=1.5, markersize=6, label='Roe格式（修正）')
plt.xlabel('x')
plt.ylabel('u')
plt.title('线性对流方程：各方法对比')
plt.legend()
plt.grid(True)

# 子图2：Roe格式与一阶迎风对比（应该重合）
plt.subplot(2, 2, 2)
plt.plot(x, u_exact, 'k-', lw=2, label='理论解', alpha=0.5)
plt.plot(x, u1_upwind, 'bo-', lw=1.5, markersize=4, label='一阶迎风', alpha=0.6)
plt.plot(x, u_roe_linear, 'm*-', lw=2, markersize=6, label='Roe格式', alpha=0.8)
plt.xlabel('x')
plt.ylabel('u')
plt.title('Roe格式 vs 一阶迎风（线性对流）- 应该重合')
plt.legend()
plt.grid(True)

# 子图3：Burgers方程 - Roe格式
plt.subplot(2, 2, 3)
plt.plot(x, u_initial, 'k--', lw=2, label='初始条件')
plt.plot(x, u_roe_burgers, 'ro-', lw=1, markersize=4, label='Roe格式（Burgers）')
plt.xlabel('x')
plt.ylabel('u')
plt.title('Burgers方程：Roe格式')
plt.legend()
plt.grid(True)

# 子图4：误差对比
plt.subplot(2, 2, 4)
error_roe = np.abs(u_roe_linear - u_exact)
error_upwind = np.abs(u1_upwind - u_exact)
error_tvd = np.abs(u2_tvd - u_exact)
error_central = np.abs(u_central - u_exact)
plt.plot(x, error_upwind, 'b-', lw=1.5, label='一阶迎风', alpha=0.7)
plt.plot(x, error_tvd, 'r-', lw=1.5, label='二阶TVD', alpha=0.7)
plt.plot(x, error_central, 'g-', lw=1.5, label='Lax-Wendroff', alpha=0.7)
plt.plot(x, error_roe, 'm-', lw=2, label='Roe格式', alpha=0.8)
plt.xlabel('x')
plt.ylabel('绝对误差')
plt.title('误差对比（Roe与迎风应该重合）')
plt.legend()
plt.grid(True)
plt.yscale('log')

plt.tight_layout()
plt.savefig('roe_comparison_fixed.png', dpi=300, bbox_inches='tight')
plt.show()