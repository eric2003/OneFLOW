import numpy as np

def initial_condition(x):
    """生成初始条件：在 [0.1, 0.3] 区间内为 1.0，其他位置为 0.0。
    
    参数:
        x (np.ndarray): 空间坐标数组
    
    返回:
        np.ndarray: 初始条件数组
    """
    return np.where((x >= 0.1) & (x <= 0.3), 1.0, 0.0)

def analytical_solution(x, t, a, domain_length):
    """计算理论解：初始波形以速度 a 平移，并应用周期边界条件。
    
    参数:
        x (np.ndarray): 空间坐标数组
        t (float): 时间
        a (float): 波速
        domain_length (float): 区域的周期长度 L
    
    返回:
        np.ndarray: 时间 t 后的波形
    """
    # 计算平移后的坐标（无边界处理）
    shifted_x = x - a * t
    
    # 应用周期边界条件，将坐标限制在 [0, domain_length) 范围内
    shifted_x_periodic = (shifted_x + domain_length) % domain_length
    
    # 返回平移后的初始条件
    return initial_condition(shifted_x_periodic)
    
    
import matplotlib.pyplot as plt

# 定义参数
L = 1.0       # 区域长度
a = 1.0       # 波速
t = 0.4       # 时间
x = np.linspace(0, L, 1000)  # 空间网格

# 计算初始条件和理论解
u0 = initial_condition(x)
u_analytical = analytical_solution(x, t, a, L)

# 可视化
plt.plot(x, u0, label="Initial condition")
plt.plot(x, u_analytical, label=f"Analytical (t={t})")
plt.legend()
plt.xlabel("x")
plt.ylabel("u")
plt.show()    