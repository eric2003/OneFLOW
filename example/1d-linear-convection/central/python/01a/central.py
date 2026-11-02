import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

# ========== 参数设置 ==========
nx = 201
L = 4.0
x = np.linspace(0, L, nx)
dx = x[1] - x[0]
c = 1.0
dt = 0.005  # 较小的时间步长
total_time = 1.0
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

# ========== 人工粘性系数函数 ==========
def compute_artificial_viscosity(u, dx, method='constant', alpha=0.5):
    """
    计算人工粘性系数
    
    参数:
        u: 当前解
        dx: 空间步长
        method: 'constant'（常数）或 'gradient'（基于梯度）
        alpha: 系数（0.1-1.0）
    
    返回:
        epsilon: 人工粘性系数数组
    """
    nx = len(u)
    epsilon = np.zeros(nx)
    
    if method == 'constant':
        # 常数系数
        epsilon[:] = alpha * abs(c) * dx
    
    elif method == 'gradient':
        # 基于梯度的自适应系数
        for i in range(nx):
            im1 = nx - 1 if i == 0 else i - 1
            ip1 = 1 if i == nx - 1 else i + 1
            
            # 计算梯度
            grad_left = abs(u[i] - u[im1]) / dx
            grad_right = abs(u[ip1] - u[i]) / dx
            grad_max = max(grad_left, grad_right)
            
            # 人工粘性与梯度成正比
            epsilon[i] = alpha * dx * grad_max
    
    elif method == 'smoothed':
        # 平滑版本（避免剧烈变化）
        epsilon_temp = np.zeros(nx)
        for i in range(nx):
            im1 = nx - 1 if i == 0 else i - 1
            ip1 = 1 if i == nx - 1 else i + 1
            
            grad_left = abs(u[i] - u[im1]) / dx
            grad_right = abs(u[ip1] - u[i]) / dx
            grad_max = max(grad_left, grad_right)
            
            epsilon_temp[i] = alpha * dx * grad_max
        
        # 平滑处理
        epsilon = alpha * epsilon_temp + (1 - alpha) * np.roll(epsilon_temp, 1)
    
    return epsilon

# ========== 带人工粘性的中心格式 ==========
def central_with_viscosity_point(u_initial, dx, dt, c, nt, epsilon_method='constant', alpha=0.5):
    """
    带人工粘性的中心格式
    
    参数:
        epsilon_method: 'constant', 'gradient', 或 'smoothed'
        alpha: 人工粘性系数（0.1-1.0）
    """
    nx = len(u_initial)
    u = u_initial.copy()
    nu = c * dt / dx
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        
        # 计算人工粘性系数
        epsilon = compute_artificial_viscosity(u, dx, epsilon_method, alpha)
        
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
            
            # 中心对流项
            convection = -0.5 * nu * (u[ip1] - u[im1])
            
            # 人工粘性项（二阶导数）
            # ∂²u/∂x² ≈ (u_{i+1} - 2u_i + u_{i-1}) / Δx²
            viscosity = (epsilon[i] * dt / dx**2) * (u[ip1] - 2*u[i] + u[im1])
            
            # 更新
            u_new[i] = u[i] + convection + viscosity
        
        u = u_new.copy()
        
        # 检测发散
        if np.max(np.abs(u)) > 10 * np.max(np.abs(u_initial)):
            print(f"⚠️  在第 {n+1} 步数值爆炸！")
            break
    
    return u

# ========== 纯中心格式（用于对比） ==========
def pure_central_point(u_initial, dx, dt, c, nt):
    nx = len(u_initial)
    u = u_initial.copy()
    nu = c * dt / dx
    
    for n in range(nt):
        u_new = np.zeros_like(u)
        for i in range(nx):
            if i == 0:
                im1 = nx - 1
                ip1 = 1
            elif i == nx - 1:
                im1 = i - 1
                ip1 = 0
            else:
                im1 = i - 1
                ip1 = i + 1
            
            u_new[i] = u[i] - 0.5 * nu * (u[ip1] - u[im1])
        
        u = u_new.copy()
        
        if np.max(np.abs(u)) > 10 * np.max(np.abs(u_initial)):
            print(f"⚠️  纯中心格式在第 {n+1} 步爆炸！")
            return u
    
    return u

# ========== 计算所有解 ==========
u_exact = np.roll(u_initial, int(round(c * total_time / dx)))

# 不同人工粘性系数的对比
u_const_01 = central_with_viscosity_point(u_initial, dx, dt, c, nt, 'constant', alpha=0.1)
u_const_05 = central_with_viscosity_point(u_initial, dx, dt, c, nt, 'constant', alpha=0.5)
u_const_10 = central_with_viscosity_point(u_initial, dx, dt, c, nt, 'constant', alpha=1.0)
u_grad = central_with_viscosity_point(u_initial, dx, dt, c, nt, 'gradient', alpha=0.5)
u_pure = pure_central_point(u_initial, dx, dt, c, nt)

# ========== 误差计算 ==========
def compute_errors(u_num, u_exact, dx):
    error_linf = np.max(np.abs(u_num - u_exact))
    error_l1 = np.sum(np.abs(u_num - u_exact)) * dx
    error_l2 = np.sqrt(np.sum((u_num - u_exact)**2) * dx)
    return error_linf, error_l1, error_l2

print("\n=== 误差统计 ===")
for name, u_sol in [
    ("纯中心格式", u_pure),
    ("人工粘性(α=0.1)", u_const_01),
    ("人工粘性(α=0.5)", u_const_05),
    ("人工粘性(α=1.0)", u_const_10),
    ("自适应粘性", u_grad)
]:
    if u_sol is not None:
        e_inf, e_l1, e_l2 = compute_errors(u_sol, u_exact, dx)
        print(f"{name:20s}: L∞={e_inf:.6f}, L1={e_l1:.6f}, L2={e_l2:.6f}")

# ========== 绘图对比 ==========
plt.figure(figsize=(14, 12))

# 子图1：不同粘性系数的对比
plt.subplot(3, 2, 1)
plt.plot(x, u_exact, 'k-', lw=2.5, label='理论解')
plt.plot(x, u_const_01, 'b--', lw=1.5, label='α=0.1', alpha=0.7)
plt.plot(x, u_const_05, 'r-', lw=1.5, label='α=0.5', alpha=0.8)
plt.plot(x, u_const_10, 'g-.', lw=1.5, label='α=1.0', alpha=0.7)
plt.xlabel('x')
plt.ylabel('u')
plt.title('不同人工粘性系数的影响')
plt.legend()
plt.grid(True)

# 子图2：自适应粘性 vs 常数粘性
plt.subplot(3, 2, 2)
plt.plot(x, u_exact, 'k-', lw=2, label='理论解', alpha=0.5)
plt.plot(x, u_const_05, 'r-', lw=1.5, label='常数粘性(α=0.5)', alpha=0.7)
plt.plot(x, u_grad, 'm--', lw=1.5, label='自适应粘性', alpha=0.8)
plt.xlabel('x')
plt.ylabel('u')
plt.title('自适应粘性 vs 常数粘性')
plt.legend()
plt.grid(True)

# 子图3：纯中心格式的灾难（如果没爆炸）
plt.subplot(3, 2, 3)
if u_pure is not None and np.max(np.abs(u_pure)) < 100:
    plt.plot(x, u_exact, 'k-', lw=2, label='理论解')
    plt.plot(x, u_pure, 'c-', lw=1.5, label='纯中心格式', alpha=0.7)
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('纯中心格式（不稳定）')
    plt.legend()
    plt.grid(True)
    plt.ylim(-2, 3)
else:
    plt.text(0.5, 0.5, '纯中心格式已爆炸', ha='center', va='center', 
             transform=plt.gca().transAxes, fontsize=14, color='red')
    plt.title('纯中心格式（不稳定）')

# 子图4：误差对比（线性坐标）
plt.subplot(3, 2, 4)
plt.plot(x, np.abs(u_const_01 - u_exact), 'b--', lw=1.5, label='α=0.1', alpha=0.7)
plt.plot(x, np.abs(u_const_05 - u_exact), 'r-', lw=1.5, label='α=0.5', alpha=0.8)
plt.plot(x, np.abs(u_const_10 - u_exact), 'g-.', lw=1.5, label='α=1.0', alpha=0.7)
plt.plot(x, np.abs(u_grad - u_exact), 'm--', lw=1.5, label='自适应', alpha=0.8)
plt.xlabel('x')
plt.ylabel('绝对误差')
plt.title('误差对比')
plt.legend()
plt.grid(True)

# 子图5：误差对比（对数坐标）
plt.subplot(3, 2, 5)
plt.semilogy(x, np.abs(u_const_01 - u_exact) + 1e-10, 'b--', lw=1.5, label='α=0.1', alpha=0.7)
plt.semilogy(x, np.abs(u_const_05 - u_exact) + 1e-10, 'r-', lw=1.5, label='α=0.5', alpha=0.8)
plt.semilogy(x, np.abs(u_const_10 - u_exact) + 1e-10, 'g-.', lw=1.5, label='α=1.0', alpha=0.7)
plt.semilogy(x, np.abs(u_grad - u_exact) + 1e-10, 'm--', lw=1.5, label='自适应', alpha=0.8)
plt.xlabel('x')
plt.ylabel('绝对误差（对数）')
plt.title('误差对比（对数坐标）')
plt.legend()
plt.grid(True, which='both')

# 子图6：人工粘性系数分布
plt.subplot(3, 2, 6)
epsilon_const = np.full(nx, 0.5 * abs(c) * dx)
epsilon_grad = compute_artificial_viscosity(u_initial, dx, 'gradient', 0.5)
plt.plot(x, epsilon_const, 'r-', lw=2, label='常数粘性', alpha=0.7)
plt.plot(x, epsilon_grad, 'm--', lw=2, label='自适应粘性', alpha=0.8)
plt.xlabel('x')
plt.ylabel('人工粘性系数 ε')
plt.title('人工粘性系数分布')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.savefig('artificial_viscosity_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

# ========== 粘性系数选择指南 ==========
print("\n" + "="*60)
print("人工粘性系数选择指南")
print("="*60)
print("""
1. 常数系数法:
   - α = 0.1-0.3: 低耗散，精度高，但可能有轻微振荡
   - α = 0.5:   平衡选择（推荐）
   - α = 1.0:   高耗散，非常稳定，但解被过度平滑
   
2. 自适应系数法:
   - 在梯度大的区域（激波、间断）自动增加粘性
   - 在光滑区域减少粘性，保持高精度
   - 推荐 α = 0.5-1.0
   
3. 经验法则:
   - 从 α = 0.5 开始尝试
   - 如果有振荡，增加 α
   - 如果解过于平滑，减小 α
   - CFL 数越小，可以用越小的 α
""")
print("="*60)