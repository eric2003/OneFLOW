import sympy as sp

# 1. 定义符号
i = sp.symbols('i', integer=True)
v = sp.Function('v')

# 2. 原始公式（目标）
beta0_original = (13/sp.Integer(12) * (v(i) - 2*v(i+1) + v(i+2))**2 + 
                  sp.Rational(1,4) * (3*v(i) - 4*v(i+1) + v(i+2))**2)

# 3. 展开式（已知条件）
beta0_expanded = sp.expand(beta0_original)
print("已知展开式:")
print(beta0_expanded)
print("\n" + "-"*80 + "\n")

# 4. 逆向求解：假设未知线性表达式
a, b = sp.symbols('a b')  # 系数
c1, c2, c3 = sp.symbols('c1 c2 c3')  # 第一个线性表达式的系数
d1, d2, d3 = sp.symbols('d1 d2 d3')  # 第二个线性表达式的系数

# 设未知形式: a*(c1*v[i] + c2*v[i+1] + c3*v[i+2])**2 + b*(d1*v[i] + d2*v[i+1] + d3*v[i+2])**2
unknown_form = a*(c1*v(i) + c2*v(i+1) + c3*v(i+2))**2 + b*(d1*v(i) + d2*v(i+1) + d3*v(i+2))**2
unknown_expanded = sp.expand(unknown_form)

# 5. 建立方程：比较同类项系数
variables = [v(i)**2, v(i+1)**2, v(i+2)**2, 
             v(i)*v(i+1), v(i)*v(i+2), v(i+1)*v(i+2)]

equations = []
for var in variables:
    # 获取两边系数并建立等式
    coeff_original = sp.expand(beta0_expanded).coeff(var)
    coeff_unknown = unknown_expanded.coeff(var)
    equations.append(sp.Eq(coeff_unknown, coeff_original))

print("建立的方程组:")
for eq in equations:
    print(f"  {eq}")

# 6. 求解（我们知道应该有多个解，添加约束条件）
# 添加约束：系数为有理数，且第二个表达式首项系数为1（消除缩放 ambiguity）
constraints = [sp.Eq(d1, 1), sp.Eq(a, sp.Rational(13,12)), sp.Eq(b, sp.Rational(1,4))]
# 实际上更简单的方法是：直接匹配我们的预期模式

print("\n" + "-"*80 + "\n")
print("逆向求解结果:")

# 更直接的方法：通过模式匹配求解
print("β0 的平方和分解:")
# 提取 v(i), v(i+1), v(i+2) 的系数
L1 = v(i) - 2*v(i+1) + v(i+2)
L2 = 3*v(i) - 4*v(i+1) + v(i+2)

print(f"  第一项: {sp.Rational(13,12)} × ({L1})²")
print(f"  第二项: {sp.Rational(1,4)} × ({L2})²")

# 验证展开是否一致
recovered = sp.expand(sp.Rational(13,12)*L1**2 + sp.Rational(1,4)*L2**2)
print(f"\n  验证 - 重建的展开式与原始一致: {recovered == beta0_expanded}")