import sympy as sp

# 定义符号
i = sp.symbols('i', integer=True)
v = sp.Function('v')
# 定义四个连续点的函数值
v_i = v(i)
v_ip1 = v(i + 1)
v_ip2 = v(i + 2)
v_ip3 = v(i + 3)

# 定义展开式（按照给定的系数）
beta0_expanded = (
    sp.Rational(2107, 240) * v_i**2
    - sp.Rational(1567, 40) * v_i * v_ip1
    + sp.Rational(3521, 120) * v_i * v_ip2
    - sp.Rational(309, 40) * v_i * v_ip3
    + sp.Rational(11003, 240) * v_ip1**2
    - sp.Rational(8623, 120) * v_ip1 * v_ip2
    + sp.Rational(2321, 120) * v_ip1 * v_ip3
    + sp.Rational(7043, 240) * v_ip2**2
    - sp.Rational(647, 40) * v_ip2 * v_ip3
    + sp.Rational(547, 240) * v_ip3**2
)

# 变量列表（为了方便提取系数）
vars_list = [v_i, v_ip1, v_ip2, v_ip3]

def complete_the_square(poly, var, vars_list):
    """
    对多项式 poly 相对于变量 var 进行配方。
    返回：(平方项, 剩余多项式)
    """
    # 提取 A (var^2 系数), D (var^1 系数), C (无 var 项)
    A = poly.coeff(var, 2)
    D = poly.coeff(var, 1)
    C = poly.as_poly(vars_list).as_expr() - A * var**2 - D * var  # 剩余部分（简化版）
    
    if A == 0:
        return 0, poly
    
    # delta = D / (2 * A)  # 注意：标准配方是 -D/(2A)，但根据二次形式 A x^2 + B x + C，delta = -B/(2A)
    # 修正：D 是 B，这里 delta = -D / (2*A)
    delta = -D / (2 * A)
    square_part = A * (var + delta)**2
    remaining = C - D**2 / (4 * A)
    
    return sp.simplify(square_part), sp.simplify(remaining)

def successive_completion(poly, vars_list, max_steps=10):
    """
    逐次配方法：反复配方直到剩余为常数或简单形式。
    选择当前二次系数绝对值最大的变量配方。
    返回：SOS 列表（平方项们）
    """
    current_poly = poly
    sos_terms = []
    steps = 0
    
    while steps < max_steps:
        # 找当前二次系数最大的变量
        max_coeff = 0
        best_var = None
        for var in vars_list:
            coeff = current_poly.coeff(var, 2)
            if abs(coeff) > max_coeff:
                max_coeff = abs(coeff)
                best_var = var
        
        if max_coeff == 0:
            # 无二次项，剩余为线性/常数（理想 SOS 为 0）
            break
        
        # 配方
        square_part, remaining = complete_the_square(current_poly, best_var, vars_list)
        sos_terms.append(square_part)
        current_poly = remaining
        steps += 1
        print(f"Step {steps}: Distributed {best_var}, square: {square_part}, remaining degree: {sp.degree(current_poly)}")
    
    if current_poly != 0:
        sos_terms.append(current_poly)  # 剩余作为最后一项（可能需进一步因子）
    
    return sos_terms

# 运行逐次配方法
print("Original beta0_expanded:")
sp.pprint(beta0_expanded)
print("\nDegree:", sp.degree(beta0_expanded))

sos_terms = successive_completion(beta0_expanded, vars_list)

print("\n=== SOS 分解（逐次配方） ===")
sos = sum(sos_terms)
sp.pprint(sos)
print("\n完整 SOS 形式：")
sp.pprint(sp.expand(sos))

# 验证：展开 SOS 与原始差值
difference = sp.simplify(sp.expand(sos) - beta0_expanded)
print("\n验证差值（应为 0）:", difference)

# 如果剩余复杂，可尝试因子化或手动调整
for term in sos_terms:
    print("\nTerm:")
    sp.pprint(sp.factor(term))