from sympy import symbols, Rational, simplify, pprint, expand

# 定义符号（v0 = v(i), v1 = v(i+1), v2 = v(i+2)）
v0, v1, v2 = symbols('v0 v1 v2')

# beta0 的展开形式
beta0 = (Rational(10,3)*v0**2 - Rational(31,3)*v0*v1 + Rational(11,3)*v0*v2 +
         Rational(25,3)*v1**2 - Rational(19,3)*v1*v2 + Rational(4,3)*v2**2)

def complete_the_square(poly, var):
    """
    对多项式poly相对于变量var进行配方。
    返回：(平方项, 剩余多项式)
    """
    # 提取系数：A (var^2), D (var^1), C (常数项)
    A = poly.coeff(var, 2)
    D = poly.coeff(var, 1)
    C = poly - A * var**2 - D * var
    
    if A == 0:
        return 0, poly  # 无二次项
    
    # 移位 delta = -D / (2 A)
    delta = -D / (2 * A)
    # 平方项
    square_part = A * (var + delta)**2
    # 剩余
    remaining = C - (D**2) / (4 * A)
    
    return simplify(square_part), simplify(remaining)

# 应用：选择 var = v1 (v(i+1)，系数最大)
var = v1
square_part, remaining = complete_the_square(beta0, var)

print("beta0 的 SOS 分解：")
print("平方项 1：")
pprint(square_part)
print("\n剩余（平方项 2）：")
pprint(remaining)
print("\n完整形式：beta0 = 平方项1 + 剩余")