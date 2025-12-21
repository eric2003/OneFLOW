from sympy import symbols, Rational, expand, latex

# 定义符号变量
v_im2, v_im1, v_i, v_ip1, v_ip2 = symbols('v_{i-2} v_{i-1} v_i v_{i+1} v_{i+2}')

# β0 的表达式
beta0_expr = Rational(13, 12) * (v_i - 2 * v_ip1 + v_ip2)**2 + Rational(1, 4) * (3 * v_i - 4 * v_ip1 + v_ip2)**2
beta0_expanded = expand(beta0_expr)

# β1 的表达式
beta1_expr = Rational(13, 12) * (v_im1 - 2 * v_i + v_ip1)**2 + Rational(1, 4) * (v_im1 - v_ip1)**2
beta1_expanded = expand(beta1_expr)

# β2 的表达式
beta2_expr = Rational(13, 12) * (v_im2 - 2 * v_im1 + v_i)**2 + Rational(1, 4) * (v_im2 - 4 * v_im1 + 3 * v_i)**2
beta2_expanded = expand(beta2_expr)

# 输出 LaTeX 公式
print("β0 = " + latex(beta0_expanded))
print("β1 = " + latex(beta1_expanded))
print("β2 = " + latex(beta2_expanded))