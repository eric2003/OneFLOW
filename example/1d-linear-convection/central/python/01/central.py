import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ========== 参数设置 ==========
nx = 101
L = 4.0
x = np.linspace(0, L, nx)
dx = x[1] - x[0]
c = 1.0
dt = 0.01  # 即使很小的时间步长也会不稳定
total_time = 0.5
nt = int(total_time / dt + 1e-8)

courant = c * dt / dx
print(f"CFL数: {courant:.3f}")

# ========== 初始条件 ==========
def initial_square_wave(x):
    u = np.zeros_like(x)
    for i in range(len(x)):
        if x[i] >= 0.8 and x[i] <= 1.2:
            u[i] = 1.0
    return u

u_initial = initial_square_wave(x)

# ========== 纯中心格式 ==========
def pure_central_point(u_initial, dx, dt, c, nt):
    """纯中心格式（显式）"""
    nx = len(u_initial)
    u = u_initial.copy()
    nu = c * dt / dx
    
    print(f"纯中心格式的放大因子 |G| = sqrt(1 + {nu:.3f}^2 * sin^2) > 1")
    print(f"最小放大因子: {np.sqrt(1 + nu**2 * 0):.6f}")
    print(f"最大放大因子: {np.sqrt(1 + nu**2 * 1):.6f}")
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            # 周期性边界
            if i == 0:
                im1 = nx - 1
                ip1 = 1
            elif i == nx - 1:
                im1 = i - 1
                ip1 = 0
            else:
                im1 = i - 1
                ip1 = i + 1
            
            # 纯中心格式
            u_new[i] = u[i] - 0.5 * nu * (u[ip1] - u[im1])
        
        u = u_new.copy()
        
        # 检测是否发散
        if np.max(np.abs(u)) > 10 * np.max(np.abs(u_initial)):
            print(f"\n⚠️  在第 {n+1} 步发生数值爆炸！")
            print(f"   最大值: {np.max(np.abs(u)):.2e}")
            break
    
    return u

# ========== 计算 ==========
u_central = pure_central_point(u_initial, dx, dt, c, nt)
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

# ========== 绘图 ==========
plt.figure(figsize=(12, 8))

# 子图1：纯中心格式的灾难
plt.subplot(2, 1, 1)
plt.plot(x, u_exact, 'k-', lw=2, label='理论解')
plt.plot(x, u_central, 'r-', lw=1.5, label='纯中心格式', alpha=0.7)
plt.xlabel('x')
plt.ylabel('u')
plt.title('纯中心格式：数值不稳定示例')
plt.legend()
plt.grid(True)
plt.ylim(-2, 3)

# 子图2：误差增长（对数坐标）
plt.subplot(2, 1, 2)
error = np.abs(u_central - u_exact)
plt.semilogy(x, error + 1e-10, 'r-', lw=1.5, label='绝对误差')
plt.xlabel('x')
plt.ylabel('误差（对数坐标）')
plt.title('误差呈指数增长')
plt.legend()
plt.grid(True, which='both')

plt.tight_layout()
plt.savefig('pure_central_instability.png', dpi=300, bbox_inches='tight')
plt.show()

# ========== 误差统计 ==========
error_linf = np.max(np.abs(u_central - u_exact))
error_l1 = np.sum(np.abs(u_central - u_exact)) * dx
print(f"\n=== 纯中心格式误差统计 ===")
print(f"L∞ 误差: {error_linf:.6f}")
print(f"L1 误差: {error_l1:.6f}")