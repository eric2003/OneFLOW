import sympy as sp

def compute_weno5_beta_k(points_values, k):
    """
    计算 WENO5 子模板 k 的光滑因子 β_k，使用 SymPy 符号积分。
    
    参数:
    - points_values: list of sympy symbols or floats, e.g., [v_im2, v_im1, vi] for k=0
    - k: int, 子模板索引 (0: 左偏, 1: 中心, 2: 右偏)
    
    返回:
    - sympy Expr: 符号表达式 of β_k
    """
    if k not in [0, 1, 2]:
        raise ValueError("k must be 0, 1, or 2 for WENO5")
    
    # 根据 k 定义点坐标 (Δx=1, 参考 x_i=0)
    if k == 0:
        coords = [-2, -1, 0]
    elif k == 1:
        coords = [-1, 0, 1]
    else:  # k=2
        coords = [0, 1, 2]
    
    x = sp.symbols('x')
    v0, v1, v2 = points_values  # 符号或数值
    
    # 拉格朗日基函数
    l0 = ((x - coords[1]) * (x - coords[2])) / ((coords[0] - coords[1]) * (coords[0] - coords[2])) * v0
    l1 = ((x - coords[0]) * (x - coords[2])) / ((coords[1] - coords[0]) * (coords[1] - coords[2])) * v1
    l2 = ((x - coords[0]) * (x - coords[1])) / ((coords[2] - coords[0]) * (coords[2] - coords[1])) * v2
    
    # 多项式 p_k(x)
    p_k = l0 + l1 + l2
    
    # 导数 (m=1 to 3)
    dp1 = sp.diff(p_k, x)
    dp2 = sp.diff(p_k, x, 2)
    dp3 = sp.diff(p_k, x, 3)  # 0 for quadratic
    
    # 积分限
    a, b = -sp.Rational(1, 2), sp.Rational(1, 2)
    
    # 各 m 项
    int_m1 = sp.integrate(dp1**2, (x, a, b))
    int_m2 = sp.integrate(dp2**2, (x, a, b))
    int_m3 = sp.integrate(dp3**2, (x, a, b))
    
    # β_k
    beta_k = int_m1 + int_m2 + int_m3
    return sp.simplify(beta_k)

# 示例使用 (符号)
v_im2, v_im1, vi = sp.symbols('v_{i-2} v_{i-1} v_i')
beta0 = compute_weno5_beta_k([v_im2, v_im1, vi], 0)
print(beta0)  # 输出: (4/3)*v_{i-2}^2 + (25/3)*v_{i-1}^2 + (10/3)*v_i^2 - (19/3)*v_{i-2}*v_{i-1} + ...