import sympy as sp

# 定义符号变量
i = sp.symbols('i', integer=True)
v = sp.Function('v')

# 定义三个beta表达式
beta0 = (13/sp.Integer(12) * (v(i) - 2*v(i+1) + v(i+2))**2 + 
         sp.Rational(1,4) * (3*v(i) - 4*v(i+1) + v(i+2))**2)

beta1 = (13/sp.Integer(12) * (v(i-1) - 2*v(i) + v(i+1))**2 + 
         sp.Rational(1,4) * (v(i-1) - v(i+1))**2)

beta2 = (13/sp.Integer(12) * (v(i-2) - 2*v(i-1) + v(i))**2 + 
         sp.Rational(1,4) * (v(i-2) - 4*v(i-1) + 3*v(i))**2)

# 展开表达式
beta0_expanded = sp.expand(beta0)
beta1_expanded = sp.expand(beta1)
beta2_expanded = sp.expand(beta2)

print(f"beta0_expanded={beta0_expanded}")
print(f"beta1_expanded={beta1_expanded}")
print(f"beta2_expanded={beta2_expanded}")

print("β0 的展开式:")
print(sp.latex(beta0_expanded))
print("\nβ1 的展开式:")
print(sp.latex(beta1_expanded))
print("\nβ2 的展开式:")
print(sp.latex(beta2_expanded))