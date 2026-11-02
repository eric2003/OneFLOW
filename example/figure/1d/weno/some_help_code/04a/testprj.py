from sympy import symbols, Rational, simplify, expand, latex, factor

# 定义符号（v0 = v(i), v1 = v(i+1), v2 = v(i+2)）
v0, v1, v2 = symbols('v0 v1 v2')

# beta0 的展开形式
beta0 = (Rational(10,3)*v0**2 - Rational(31,3)*v0*v1 + Rational(11,3)*v0*v2 +
         Rational(25,3)*v1**2 - Rational(19,3)*v1*v2 + Rational(4,3)*v2**2)

def complete_the_square(poly, var):
    """
    对多项式 poly 相对于变量 var 进行配方。
    返回：(平方项, 剩余多项式)
    """
    A = poly.coeff(var, 2)
    D = poly.coeff(var, 1)
    C = poly - A * var**2 - D * var
    
    if A == 0:
        return 0, poly
    
    # delta = D / (2 * A)
    delta = D / (2 * A)
    square_part = A * (var + delta)**2
    remaining = C - (D**2) / (4 * A)
    
    return simplify(square_part), simplify(remaining)

# 应用：选择 var = v1 (系数最大)
var = v1
square_part, remaining = complete_the_square(beta0, var)

# 简化剩余为平方形式
remaining_squared = factor(remaining)

# 修正：SOS 使用简化平方
sos = square_part + remaining_squared

# 验证
expanded_sos = expand(sos)
difference = simplify(beta0 - expanded_sos)
print("验证：展开 SOS 与原始差值 =", difference)  # 应为 0

# LaTeX 输出
print("\n原始公式 LaTeX:")
print(latex(beta0))
print("\n平方项 LaTeX:")
print(latex(square_part))
print("\n剩余项 LaTeX (简化平方):")
print(latex(remaining_squared))
print("\n完整 SOS LaTeX:")
print(latex(sos))

# 目标形式验证
target = Rational(13,12) * (v0 - 2*v1 + v2)**2 + Rational(1,4) * (3*v0 - 4*v1 + v2)**2
target_expanded = expand(target)
target_diff = simplify(beta0 - target_expanded)
print("\n目标形式与原始差值 =", target_diff)  # 应为 0