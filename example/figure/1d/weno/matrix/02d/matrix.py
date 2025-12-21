import sympy as sp
from sympy import symbols, Rational, Matrix, latex

# ---------- 全局符号（保持一致性） ----------
r, i_sym = symbols('r i', integer=True)
xi = symbols('\\xi')

# ---------- 核心函数 ----------
def alpha(j): return -r + j - Rational(1, 2)
def beta(j):  return -r + j + Rational(1, 2)

def vbar(idx_expr):
    """生成 \overline{v}_{expr} 的 SymPy 符号（用于计算）"""
    # 使用带 LaTeX 名称的 Symbol，确保 latex() 输出正确
    name = r"\overline{v}_{" + latex(idx_expr) + r"}"
    return sp.Symbol(name)

def build_v_vector(k_val):
    """构建 v = [v_{i-r}, v_{i-r+1}, ..., v_{i-r+k-1}]^T"""
    return Matrix([vbar(i_sym - r + j) for j in range(k_val)])

def build_M_matrix(k_val, r_is_symbolic=True, specific_r=None):
    """
    构建 M_{j,m} = ∫_{α_j}^{β_j} ξ^m dξ,  j,m = 0,...,k-1
    
    参数:
        k_val (int): 多项式阶数 k
        r_is_symbolic (bool): 是否保留 r 为符号
        specific_r (int): 若 r_is_symbolic=False，指定具体 r 值
    """
    if not r_is_symbolic and specific_r is None:
        raise ValueError("specific_r must be provided if r_is_symbolic=False")
    
    M = sp.zeros(k_val, k_val)
    r_use = r if r_is_symbolic else specific_r
    
    for j in range(k_val):
        a_j = (-r_use + j - Rational(1, 2))
        b_j = (-r_use + j + Rational(1, 2))
        for m in range(k_val):
            # ∫ ξ^m dξ = (b^{m+1} - a^{m+1}) / (m+1)
            integral = (b_j**(m + 1) - a_j**(m + 1)) / Rational(m + 1)
            M[j, m] = sp.simplify(integral)
    return M

def solve_coefficients(k_val, r_val=None):
    """
    求解 a = M^{-1} v
    
    返回:
        a_vec: 系数向量 [a0, a1, ..., a_{k-1}]^T （SymPy Matrix）
        M: 使用的 M 矩阵
        v: 使用的 v 向量
    """
    if r_val is None:
        # 保留 r 为符号（第 1 层）
        M = build_M_matrix(k_val, r_is_symbolic=True)
        v = build_v_vector(k_val)
    else:
        # 固定 r（第 2 层）
        M = build_M_matrix(k_val, r_is_symbolic=False, specific_r=r_val)
        # 构建 v 向量（此时 r 已固定）
        v = Matrix([vbar(i_sym - r_val + j) for j in range(k_val)])
    
    # 求逆并计算 a = M^{-1} v
    M_inv = M.inv()  # SymPy 会自动简化
    a_vec = M_inv * v
    return a_vec, M, v

# ---------- 工具：美化 LaTeX 输出 ----------
def print_solution_for_k_r(k_val, r_val=None):
    """打印 a = M^{-1} v 的完整 LaTeX（通用或具体）"""
    a_vec, M, v = solve_coefficients(k_val, r_val)
    
    print(f"\n=== WENO Coefficients: k={k_val}" + (f", r={r_val}" if r_val is not None else "(r symbolic)") + " ===")
    print("M =")
    print(latex(M))
    print("\nv =")
    print(latex(v))
    print("\na = M^{-1} v =")
    print(latex(a_vec))
    return a_vec, M, v

# ---------- 主程序：按需生成 ----------
if __name__ == "__main__":
    # 示例 1: k=3, r 保持符号（通用表达式）
    print_solution_for_k_r(k_val=3, r_val=None)
    
    # 示例 2: k=3, r=0,1,2（具体表达式）
    for r_example in [0, 1, 2]:
        print_solution_for_k_r(k_val=3, r_val=r_example)