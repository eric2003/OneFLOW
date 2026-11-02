import sympy as sp

def construct_system(k, r_symbolic=True, specific_r=None):
    """
    构造积分约束系统：∫_{α}^{β} p(r, φ) dφ = v̄_{i−r+j}
    
    参数：
        k (int): 多项式阶数 (p 是 k−1 次多项式)
        r_symbolic (bool): 是否将 r 保留为符号（True）还是使用 specific_r（False）
        specific_r (int): 若 r_symbolic=False，指定具体的 r 值 (0 ≤ r ≤ k−1)
    
    返回：
        A (Matrix): k×k 系数矩阵 A_{j,m} = ∫ φ^m dφ over [α(r,j), β(r,j)]
        rhs_symbols (list): 右端符号 [v̄_{i−r}, v̄_{i−r+1}, ..., v̄_{i−r+k−1}]
        phi (Symbol): 积分变量（内部用 phi，LaTeX 可替换）
        r (Symbol or int): 使用的 r
    """
    # 积分变量（内部用 phi，LaTeX 显示为 \phi）
    phi = sp.Symbol('phi')
    
    if r_symbolic:
        r = sp.Symbol('r', integer=True)
    else:
        if specific_r is None:
            raise ValueError("specific_r must be provided if r_symbolic=False")
        if not (0 <= specific_r <= k - 1):
            raise ValueError(f"specific_r={specific_r} must be in [0, {k-1}]")
        r = specific_r

    # 多项式系数 a_0(r), ..., a_{k-1}(r)
    a = sp.symbols(f'a0:{k}', cls=sp.Function)
    # 注意：这里 a0(r) 是函数形式，保留对 r 的依赖
    
    # 构造多项式 p(r, phi) = Σ_{m=0}^{k-1} a_m(r) * phi^m
    p = sum(a[m](r) * phi**m for m in range(k))
    
    # 右端符号：v̄_{i - r + j}
    v = sp.symbols(f'vbar0:{k}')  # vbar0, vbar1, ..., vbar_{k-1}
    rhs_symbols = [v[j] for j in range(k)]  # 对应 j=0,...,k−1
    
    # 构建系数矩阵 A：A[j, m] = ∫_{α}^{β} φ^m dφ
    A = sp.zeros(k, k)
    
    for j in range(k):
        # α(r,j) = -r + j - 1/2, β(r,j) = -r + j + 1/2
        alpha = -r + j - sp.Rational(1, 2)
        beta  = -r + j + sp.Rational(1, 2)
        
        for m in range(k):
            # 积分 ∫ φ^m dφ = (β^{m+1} - α^{m+1}) / (m+1)
            integral = (beta**(m + 1) - alpha**(m + 1)) / sp.Rational(m + 1)
            A[j, m] = sp.simplify(integral)
    
    return A, rhs_symbols, phi, r, a

# -----------------------------
# 示例 1：通用表达式（k=3，r 为符号）
# -----------------------------
print("=== 通用表达式 (k=3, r symbolic) ===")
A_gen, rhs_gen, phi, r_sym, a_sym = construct_system(k=3, r_symbolic=True)

# 打印矩阵 A（LaTeX，用 \phi 替代 \xi）
# SymPy 默认用 phi，若原用 xi 则需替换，但这里直接定义为 phi
print("系数矩阵 A (k=3):")
sp.pprint(A_gen)
print("\nLaTeX (A):")
latex_A = sp.latex(A_gen)
latex_A_phi = latex_A.replace(r'\xi', r'\phi')  # 虽然这里没用 xi，但保险起见
print(latex_A_phi)

# -----------------------------
# 示例 2：具体 r=0, k=3
# -----------------------------
print("\n=== 具体情况 (k=3, r=0) ===")
A_r0, rhs_r0, phi, r_val, a_val = construct_system(k=3, r_symbolic=False, specific_r=0)
sp.pprint(A_r0)
print("\nLaTeX (A for r=0):")
print(sp.latex(A_r0))

# -----------------------------
# 示例 3：保留接口以供将来代入 r 值
# -----------------------------
print("\n=== 从通用表达式代入 r=1 ===")
A_r1 = A_gen.subs(r_sym, 1)
sp.pprint(A_r1)
print("\nLaTeX (A for r=1 via subs):")
print(sp.latex(A_r1))