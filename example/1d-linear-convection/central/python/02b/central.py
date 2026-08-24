import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ========== 隐式中心格式（完全修正版）==========
def implicit_central_crank_nicolson(u_initial, dx, dt, c, nt):
    """
    隐式中心格式（Crank-Nicolson）- 完全修正版
    使用正确的矩阵构建方法
    """
    nx = len(u_initial)
    u = np.zeros((nx, nt))
    u[:, 0] = u_initial.copy()
    
    alpha = c * dt / (4 * dx)
    
    print(f"隐式中心格式: α={alpha:.6f}, CFL={c*dt/dx:.3f}")
    
    # 构建周期性三对角矩阵
    main_diag = np.full(nx, 1.0)
    upper_diag = np.full(nx, alpha)
    lower_diag = np.full(nx, -alpha)
    
    A = diags([lower_diag, main_diag, upper_diag], 
              offsets=[-1, 0, 1], 
              shape=(nx, nx), 
              format='csr')
    
    # 周期性边界条件
    A[0, -1] = -alpha
    A[-1, 0] = alpha
    
    for n in range(nt - 1):
        rhs = np.zeros(nx)
        for i in range(nx):
            im1 = (i - 1) % nx
            ip1 = (i + 1) % nx
            rhs[i] = u[i, n] - alpha * (u[ip1, n] - u[im1, n])
        
        try:
            u_new = spsolve(A, rhs)
            u[:, n+1] = u_new
        except Exception as e:
            print(f"  ⚠️  求解失败: {e}")
            return u
        
        if np.max(np.abs(u[:, n+1])) > 1e10:
            print(f"  ⚠️  数值溢出！")
            return u
    
    return u

# ========== Rusanov 格式 ==========
def rusanov_scalar_point(u_initial, dx, dt, c, nt):
    """Rusanov格式（局部Lax-Friedrichs）"""
    nx = len(u_initial)
    u = u_initial.copy()
    lambda_max = abs(c)
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            im1 = (i - 1) % nx
            ip1 = (i + 1) % nx
            
            u_L = u[i]
            u_R = u[ip1]
            f_L = c * u_L
            f_R = c * u_R
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

# ========== 参数设置 ==========
nx = 101
L = 4.0
x = np.linspace(0, L, nx)
dx = x[1] - x[0]
c = 1.0

# 测试小CFL数
cfl_numbers = [0.1, 0.2, 0.3]
total_time = 1.0

# 初始条件（光滑正弦波）
def initial_smooth_wave(x):
    # 在x∈[1,2]区间内生成一个正弦波包
    u = np.zeros_like(x)
    mask = (x >= 1.0) & (x <= 2.0)
    u[mask] = np.sin(np.pi * (x[mask] - 1.0))  # 周期为1的正弦波
    return u

u_initial = initial_smooth_wave(x)
# 精确解：正弦波以速度c平移
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

# ========== 计算所有格式 ==========
print("="*60)
print("小CFL数下隐式中心格式 vs Rusanov（光滑初始条件）")
print("="*60)

results = {}

for cfl in cfl_numbers:
    dt = cfl * dx / c
    nt = int(total_time / dt + 1e-8)
    
    print(f"\nCFL = {cfl:.2f}, dt = {dt:.6f}, nt = {nt}")
    
    # 隐式中心格式
    u_implicit_matrix = implicit_central_crank_nicolson(u_initial, dx, dt, c, nt)
    e_inf_imp = np.max(np.abs(u_implicit_matrix[:, -1] - u_exact))
    e_l1_imp = np.sum(np.abs(u_implicit_matrix[:, -1] - u_exact)) * dx
    
    # Rusanov格式
    u_rusanov = rusanov_scalar_point(u_initial, dx, dt, c, nt)
    e_inf_rus = np.max(np.abs(u_rusanov - u_exact))
    e_l1_rus = np.sum(np.abs(u_rusanov - u_exact)) * dx
    
    results[cfl] = {
        'implicit_u': u_implicit_matrix[:, -1],
        'implicit_einf': e_inf_imp,
        'implicit_el1': e_l1_imp,
        'rusanov_u': u_rusanov,
        'rusanov_einf': e_inf_rus,
        'rusanov_el1': e_l1_rus
    }
    
    print(f"  隐式中心: L∞={e_inf_imp:.6f}, L1={e_l1_imp:.6f}")
    print(f"  Rusanov:  L∞={e_inf_rus:.6f}, L1={e_l1_rus:.6f}")
    if e_inf_imp > 0:
        ratio = e_inf_rus / e_inf_imp
        print(f"  精度比: Rusanov/隐式 = {ratio:.2f}x")

# ========== 绘图对比 ==========
plt.figure(figsize=(12, 8))

# 子图1：解的对比
plt.subplot(2, 2, 1)
for cfl in [0.1, 0.2, 0.3]:
    plt.plot(x, results[cfl]['implicit_u'], '--', lw=1.5, 
             label=f'隐式中心(CFL={cfl})', alpha=0.7)
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解', alpha=0.9)
plt.xlabel('x')
plt.ylabel('u')
plt.title('隐式中心格式解对比')
plt.legend()
plt.grid(True)

# 子图2：Rusanov解对比
plt.subplot(2, 2, 2)
for cfl in [0.1, 0.2, 0.3]:
    plt.plot(x, results[cfl]['rusanov_u'], '-', lw=1.5, 
             label=f'Rusanov(CFL={cfl})', alpha=0.7)
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解', alpha=0.9)
plt.xlabel('x')
plt.ylabel('u')
plt.title('Rusanov格式解对比')
plt.legend()
plt.grid(True)

# 子图3：误差对比（线性坐标）
plt.subplot(2, 2, 3)
for cfl in cfl_numbers:
    error_imp = np.abs(results[cfl]['implicit_u'] - u_exact)
    error_rus = np.abs(results[cfl]['rusanov_u'] - u_exact)
    plt.plot(x, error_imp, '--', lw=1.5, label=f'隐式(CFL={cfl})', alpha=0.6)
    plt.plot(x, error_rus, '-', lw=1.5, label=f'Rusanov(CFL={cfl})', alpha=0.8)
plt.xlabel('x')
plt.ylabel('绝对误差')
plt.title('误差对比（线性坐标）')
plt.legend()
plt.grid(True)

# 子图4：误差对比（对数坐标）
plt.subplot(2, 2, 4)
for cfl in cfl_numbers:
    error_imp = np.abs(results[cfl]['implicit_u'] - u_exact) + 1e-10
    error_rus = np.abs(results[cfl]['rusanov_u'] - u_exact) + 1e-10
    plt.semilogy(x, error_imp, '--', lw=1.5, label=f'隐式(CFL={cfl})', alpha=0.6)
    plt.semilogy(x, error_rus, '-', lw=1.5, label=f'Rusanov(CFL={cfl})', alpha=0.8)
plt.xlabel('x')
plt.ylabel('绝对误差（对数）')
plt.title('误差对比（对数坐标）')
plt.legend()
plt.grid(True, which='both')

plt.tight_layout()
plt.savefig('implicit_vs_rusanov_smooth.png', dpi=300, bbox_inches='tight')
plt.show()

# ========== 误差统计表格 ==========
print("\n" + "="*60)
print("误差统计对比（光滑初始条件）")
print("="*60)
print(f"{'CFL':<8} {'隐式 L∞':<15} {'Rusanov L∞':<15} {'精度比':<10}")
print("-"*60)
for cfl in cfl_numbers:
    r = results[cfl]
    ratio = r['rusanov_einf'] / r['implicit_einf'] if r['implicit_einf'] > 0 else 0
    print(f"{cfl:<8.2f} {r['implicit_einf']:<15.6f} {r['rusanov_einf']:<15.6f} {ratio:<10.2f}x")
print("="*60)