import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
import time

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ========== 隐式中心格式（Crank-Nicolson）==========
def implicit_central_crank_nicolson(u_initial, dx, dt, c, nt):
    """
    隐式中心格式（Crank-Nicolson）
    无条件稳定，二阶精度
    
    参数:
        u_initial: 初始条件
        dx: 空间步长
        dt: 时间步长
        c: 对流速度
        nt: 时间步数
    
    返回:
        u: 最终解（nx × nt 的矩阵，每列是一个时间步）
    """
    nx = len(u_initial)
    u = np.zeros((nx, nt))
    u[:, 0] = u_initial.copy()
    
    alpha = c * dt / (4 * dx)  # 关键系数
    nu = c * dt / dx  # CFL数
    
    print(f"隐式中心格式: α={alpha:.6f}, CFL={nu:.3f}")
    print(f"网格数: {nx}, 时间步数: {nt}")
    
    # 构建三对角矩阵（周期性边界）
    # 使用 scipy 的 banded 格式存储：[上对角线, 主对角线, 下对角线]
    # 对于周期性边界，需要特殊处理
    
    for n in range(nt - 1):
        # 构建右端向量 RHS
        rhs = np.zeros(nx)
        for i in range(nx):
            im1 = (i - 1) % nx
            ip1 = (i + 1) % nx
            
            # RHS_i = α*u_{i-1}^n + u_i^n - α*u_{i+1}^n
            rhs[i] = alpha * u[im1, n] + u[i, n] - alpha * u[ip1, n]
        
        # 求解线性系统 A*u^{n+1} = rhs
        # 方法1：使用 scipy 的循环三对角求解器
        u[:, n+1] = solve_periodic_tridiagonal(nx, alpha, rhs)
        
        # 检查数值稳定性
        if np.max(np.abs(u[:, n+1])) > 1e10:
            print(f"⚠️  数值不稳定！")
            break
    
    return u

# ========== 求解周期性三对角系统 ==========
def solve_periodic_tridiagonal(n, alpha, rhs):
    """
    求解周期性三对角线性系统
    
    矩阵形式:
    [ 1,  α,  0, ..., -α] [u0]   [rhs0]
    [-α,  1,  α, ...,  0] [u1]   [rhs1]
    [ 0, -α,  1, ...,  0] [u2] = [rhs2]
    [..., ..., ..., ..., ...] [...]   [...]
    [ α,  0,  0, ...,  1] [un]   [rhsn]
    """
    # 使用 Thomas 算法的变体（Sherman-Morrison 公式）
    
    # 构建非周期系统的三对角矩阵
    a = np.full(n, -alpha)  # 下对角线
    b = np.full(n, 1.0)     # 主对角线
    c = np.full(n, alpha)   # 上对角线
    
    # 修正周期性边界条件
    # 使用 Sherman-Morrison 公式处理角落元素
    # A = A' + uv^T，其中 u = [α, 0, ..., 0, -α]^T, v = [1, 0, ..., 0, -1]^T
    
    # 步骤1：求解 A'y = rhs
    y = thomas_algorithm(a, b, c, rhs)
    
    # 步骤2：求解 A'z = u
    u_vec = np.zeros(n)
    u_vec[0] = alpha
    u_vec[-1] = -alpha
    z = thomas_algorithm(a, b, c, u_vec)
    
    # 步骤3：求解 A'w = v
    v_vec = np.zeros(n)
    v_vec[0] = 1.0
    v_vec[-1] = -1.0
    w = thomas_algorithm(a, b, c, v_vec)
    
    # 步骤4：Sherman-Morrison 修正
    denom = 1.0 + np.dot(v_vec, z)
    if abs(denom) < 1e-10:
        # 退化情况，使用直接求解
        return solve_direct_periodic(n, alpha, rhs)
    
    factor = np.dot(v_vec, y) / denom
    solution = y - factor * z
    
    return solution

# ========== Thomas 算法（三对角求解器）==========
def thomas_algorithm(a, b, c, d):
    """
    Thomas 算法求解三对角系统
    a: 下对角线 (长度 n，a[0] 未使用)
    b: 主对角线 (长度 n)
    c: 上对角线 (长度 n，c[n-1] 未使用)
    d: 右端向量 (长度 n)
    """
    n = len(d)
    
    # 前向消元
    for i in range(1, n):
        w = a[i] / b[i-1]
        b[i] -= w * c[i-1]
        d[i] -= w * d[i-1]
    
    # 回代
    x = np.zeros(n)
    x[-1] = d[-1] / b[-1]
    for i in range(n-2, -1, -1):
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    
    return x

# ========== 直接求解（备用方法）==========
def solve_direct_periodic(n, alpha, rhs):
    """
    直接构建并求解周期性三对角系统（用于小规模问题）
    """
    # 构建完整矩阵
    A = np.zeros((n, n))
    for i in range(n):
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        A[i, im1] = -alpha
        A[i, i] = 1.0
        A[i, ip1] = alpha
    
    # 直接求解
    solution = np.linalg.solve(A, rhs)
    return solution

# ========== 显式中心格式（用于对比）==========
def explicit_central(u_initial, dx, dt, c, nt):
    """显式中心格式（不稳定）"""
    nx = len(u_initial)
    u = u_initial.copy()
    nu = c * dt / dx
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            im1 = (i - 1) % nx
            ip1 = (i + 1) % nx
            u_new[i] = u[i] - 0.5 * nu * (u[ip1] - u[im1])
        u = u_new.copy()
        
        if np.max(np.abs(u)) > 1e10:
            print(f"⚠️  显式中心格式在第 {n+1} 步爆炸！")
            return u
    
    return u

# ========== 参数设置 ==========
nx = 101
L = 4.0
x = np.linspace(0, L, nx)
dx = x[1] - x[0]
c = 1.0

# 测试不同 CFL 数
cfl_numbers = [0.25, 0.5, 1.0, 2.0, 5.0]
dt_base = 0.01
total_time = 1.0

# 初始条件
def initial_square_wave(x):
    u = np.zeros_like(x)
    for i in range(len(x)):
        if x[i] >= 0.8 and x[i] <= 1.2:
            u[i] = 1.0
    return u

u_initial = initial_square_wave(x)
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

# ========== 测试不同 CFL 数 ==========
print("="*60)
print("测试隐式中心格式在不同 CFL 数下的表现")
print("="*60)

results = []

for cfl in cfl_numbers:
    dt = cfl * dx / c
    nt = int(total_time / dt + 1e-8)
    
    print(f"\nCFL = {cfl:.2f}, dt = {dt:.6f}, nt = {nt}")
    
    # 隐式中心格式
    start_time = time.time()
    u_implicit = implicit_central_crank_nicolson(u_initial, dx, dt, c, nt)
    elapsed_implicit = time.time() - start_time
    
    # 显式中心格式（仅小 CFL 测试）
    if cfl <= 0.5:
        start_time = time.time()
        u_explicit = explicit_central(u_initial, dx, dt, c, nt)
        elapsed_explicit = time.time() - start_time
    else:
        u_explicit = None
        elapsed_explicit = None
    
    # 计算误差
    e_inf_imp = np.max(np.abs(u_implicit[:, -1] - u_exact))
    e_l1_imp = np.sum(np.abs(u_implicit[:, -1] - u_exact)) * dx
    
    results.append({
        'CFL': cfl,
        'dt': dt,
        'nt': nt,
        'implicit_Linf': e_inf_imp,
        'implicit_L1': e_l1_imp,
        'implicit_time': elapsed_implicit,
        'explicit_time': elapsed_explicit
    })
    
    print(f"  隐式: L∞={e_inf_imp:.6f}, 时间={elapsed_implicit:.4f}s")
    if u_explicit is not None:
        e_inf_exp = np.max(np.abs(u_explicit - u_exact))
        print(f"  显式: L∞={e_inf_exp:.6f}, 时间={elapsed_explicit:.4f}s")

# ========== 绘图对比 ==========
plt.figure(figsize=(14, 10))

# 子图1：不同 CFL 数的误差对比
plt.subplot(2, 2, 1)
cfl_vals = [r['CFL'] for r in results]
linf_vals = [r['implicit_Linf'] for r in results]
plt.plot(cfl_vals, linf_vals, 'o-', lw=2, markersize=8, label='隐式中心格式')
plt.xlabel('CFL 数')
plt.ylabel('L∞ 误差')
plt.title('隐式中心格式：误差 vs CFL 数')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')

# 子图2：计算时间对比
plt.subplot(2, 2, 2)
implicit_times = [r['implicit_time'] for r in results]
plt.plot(cfl_vals, implicit_times, 'o-', lw=2, markersize=8, 
         color='red', label='隐式中心格式')
plt.xlabel('CFL 数')
plt.ylabel('计算时间 (秒)')
plt.title('计算时间 vs CFL 数')
plt.grid(True)
plt.xscale('log')

# 子图3：CFL=5.0 时的解对比
plt.subplot(2, 2, 3)
cfl_5_idx = [i for i, r in enumerate(results) if r['CFL'] == 5.0][0]
u_imp_5 = implicit_central_crank_nicolson(
    u_initial, dx, results[cfl_5_idx]['dt'], c, results[cfl_5_idx]['nt']
)
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解')
plt.plot(x, u_imp_5[:, -1], 'r--', lw=1.5, label=f'隐式中心 (CFL=5.0)', alpha=0.8)
plt.xlabel('x')
plt.ylabel('u')
plt.title('大 CFL 数下的解 (CFL=5.0)')
plt.legend()
plt.grid(True)

# 子图4：误差分布（CFL=0.25）
plt.subplot(2, 2, 4)
cfl_025_idx = [i for i, r in enumerate(results) if r['CFL'] == 0.25][0]
u_imp_025 = implicit_central_crank_nicolson(
    u_initial, dx, results[cfl_025_idx]['dt'], c, results[cfl_025_idx]['nt']
)
error_imp = np.abs(u_imp_025[:, -1] - u_exact)
plt.plot(x, error_imp, 'r-', lw=1.5, label='隐式中心', alpha=0.8)
plt.xlabel('x')
plt.ylabel('绝对误差')
plt.title('误差分布 (CFL=0.25)')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('implicit_central_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# ========== 结果表格 ==========
print("\n" + "="*60)
print("隐式中心格式性能总结")
print("="*60)
print(f"{'CFL':<8} {'dt':<10} {'nt':<8} {'L∞误差':<12} {'计算时间(s)':<12}")
print("-"*60)
for r in results:
    print(f"{r['CFL']:<8.2f} {r['dt']:<10.6f} {r['nt']:<8d} "
          f"{r['implicit_Linf']:<12.6f} {r['implicit_time']:<12.4f}")
print("="*60)

# ========== 与显式格式对比 ==========
print("\n" + "="*60)
print("隐式 vs 显式中心格式对比 (CFL=0.25)")
print("="*60)

dt_test = 0.25 * dx / c
nt_test = int(total_time / dt_test + 1e-8)

# 显式
u_exp = explicit_central(u_initial, dx, dt_test, c, nt_test)
e_exp = np.max(np.abs(u_exp - u_exact))

# 隐式
u_imp = implicit_central_crank_nicolson(u_initial, dx, dt_test, c, nt_test)
e_imp = np.max(np.abs(u_imp[:, -1] - u_exact))

print(f"显式中心: L∞={e_exp:.6f} (不稳定，仅小CFL可用)")
print(f"隐式中心: L∞={e_imp:.6f} (稳定，任意CFL可用)")
print(f"隐式精度提升: {e_exp/e_imp:.2f} 倍")
print("="*60)