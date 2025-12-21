import sympy as sp

# 定义符号
i = sp.symbols('i', integer=True)
v = sp.Function('v')

# 定义四个连续点的函数值
v0 = v(i)
v1 = v(i+1)
v2 = v(i+2)
v3 = v(i+3)

# 定义三个线性表达式
L1 = v0 - 3*v1 + 3*v2 - v3          # 三阶差分
L2 = 2*v0 - 5*v1 + 4*v2 - v3        # WENO专用模板
L3 = 43*v0 - 69*v1 + 33*v2 - 7*v3   # 高阶修正项

# 定义完整的β₀表达式
beta0 = (sp.Rational(1043,960) * L1**2 + 
         sp.Rational(13,12) * L2**2 + 
         sp.Rational(1,288) * L3 * L1 + 
         sp.Rational(1,576) * L3**2)

print("="*80)
print("β₀ 原始表达式结构")
print("="*80)
print(f"项1: 1043/960 × ({L1})²")
print(f"项2: 13/12 × ({L2})²")
print(f"项3: 1/288 × ({L3}) × ({L1})")
print(f"项4: 1/576 × ({L3})²")
print("\n" + "="*80)
print("完全展开式")
print("="*80)

# 完全展开
beta0_expanded = sp.expand(beta0)

# 输出LaTeX公式
latex_output = sp.latex(beta0_expanded)
print(f"LaTeX公式:")
print(f"\\[")
print(f"\\beta_0 = {latex_output}")
print(f"\\]")

# 为了更清晰地显示，按变量分组
print("\n" + "="*80)
print("按变量分组的展开式")
print("="*80)

# 收集同类项
terms = {
    'v0²': beta0_expanded.coeff(v0**2),
    'v1²': beta0_expanded.coeff(v1**2),
    'v2²': beta0_expanded.coeff(v2**2),
    'v3²': beta0_expanded.coeff(v3**2),
    'v0v1': beta0_expanded.coeff(v0*v1),
    'v0v2': beta0_expanded.coeff(v0*v2),
    'v0v3': beta0_expanded.coeff(v0*v3),
    'v1v2': beta0_expanded.coeff(v1*v2),
    'v1v3': beta0_expanded.coeff(v1*v3),
    'v2v3': beta0_expanded.coeff(v2*v3),
}

for term, coeff in terms.items():
    if coeff != 0:
        print(f"{term}: {coeff} = {float(coeff):.6f}")

# 输出完整的数学表达式
print("\n" + "="*80)
print("完整的数学表达式")
print("="*80)
print("\\beta_0 = ")
for term, coeff in terms.items():
    if coeff != 0:
        sign = "+" if coeff > 0 else ""
        print(f"  {sign}{coeff} \\cdot {term.replace('v', 'v_{i+').replace('²', '}^2').replace('v0', 'v_i').replace('v1', 'v_{i+1}').replace('v2', 'v_{i+2}').replace('v3', 'v_{i+3}')}")