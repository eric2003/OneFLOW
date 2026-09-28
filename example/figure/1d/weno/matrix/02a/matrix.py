import sympy as sp
from sympy import symbols, Function, latex
from sympy.printing.latex import LatexPrinter

# ---- Step 1: 定义基础符号（可全局配置） ----
k = sp.Symbol('k', integer=True, positive=True)
r, j, xi = sp.symbols('r j \\xi')  # 注意：xi 命名为 \xi 以便 LaTeX

# 多项式系数：a_m(r) 作为函数
a = Function('a')
# 为简化，我们用 a(0)(r), a(1)(r), ... 表示 a_0(r), a_1(r), ...
# 但 LaTeX 中希望显示为 a_0(r)，需自定义打印

# 右端符号：定义一个带下标的符号生成器
def vbar(index_expr):
    # 使用 Dummy symbol with name containing LaTeX
    name = f"\\overline{{v}}_{{{sp.latex(index_expr)}}}"
    return sp.Symbol(name)

# ---- Step 2: 构造多项式 p(r, xi) ----
def build_polynomial(order, use_r_in_coeff=True):
    terms = []
    for m in range(order):
        if use_r_in_coeff:
            coeff = a(m)(r)  # a_m(r)
        else:
            coeff = sp.Symbol(f'a_{m}')  # a_m（无 r 依赖）
        terms.append(coeff * xi**m)
    return sum(terms)

# ---- Step 3: 构造积分限 α, β ----
alpha = -r + j - sp.Rational(1, 2)
beta  = -r + j + sp.Rational(1, 2)

# ---- Step 4: 构造各个表达式（保持符号） ----
# 多项式（带 r 依赖）
p_r_xi = build_polynomial(k, use_r_in_coeff=True)

# 积分等式（抽象）
integral_eq_abstract = sp.Eq(
    sp.Integral(p_r_xi, (xi, alpha, beta)),
    vbar(sp.Symbol('i') - r + j)
)

# 积分等式（显式展开，带 a_m(r)）
p_expanded = p_r_xi
integral_eq_expanded = sp.Eq(
    sp.Integral(p_expanded, (xi, alpha, beta)),
    vbar(sp.Symbol('i') - r + j)
)

# 积分等式（a_m 不带 r）
p_no_r = build_polynomial(k, use_r_in_coeff=False)
integral_eq_no_r = sp.Eq(
    sp.Integral(p_no_r, (xi, alpha, beta)),
    vbar(sp.Symbol('i') - r + j)
)

# ---- Step 5: 生成 LaTeX 输出（array 环境） ----
def generate_latex_array(*exprs, alpha_expr, beta_expr):
    lines = []
    for expr in exprs:
        lines.append(latex(expr))
    
    # 添加 α, β 定义（j 范围）
    j_range = r'j=0,1,\ldots,k-1'
    alpha_line = f"{latex(alpha_expr)}={latex(alpha)},\\quad {j_range}"
    beta_line  = f"{latex(beta_expr)}={latex(beta)},\\quad {j_range}"
    
    full_latex = r"\begin{array}{l}" + "\n"
    full_latex += " \\\\\n".join(lines + [alpha_line, beta_line])
    full_latex += "\n\\end{array}"
    return full_latex

# ---- 执行生成 ----
# 定义 α(r,j), β(r,j) 作为符号（用于等式左边）
alpha_sym = sp.Symbol(r'\alpha(r,j)')
beta_sym  = sp.Symbol(r'\beta(r,j)')

latex_output = generate_latex_array(
    sp.Eq(sp.Symbol('p(r,\\xi)'), p_r_xi),
    sp.Eq(sp.Integral(sp.Symbol('p(r,\\xi)'), (xi, alpha_sym, beta_sym)), vbar(sp.Symbol('i') - r + j)),
    integral_eq_expanded,
    integral_eq_no_r,
    alpha_expr=alpha_sym,
    beta_expr=beta_sym
)

print(latex_output)