import sympy as sp
from sympy import symbols, Function, latex, Rational

# ---------- 配置 ----------
# 当前仅用于 LaTeX 显示，不用于求和
k = sp.Symbol('k', integer=True, positive=True)
r, j, xi = symbols('r j \\xi')
i_sym = symbols('i')  # 用于 vbar 下标

# 多项式系数函数（用于具体 k 时）
a = Function('a')

# ---------- 工具函数 ----------
def vbar_expr():
    """返回 \overline{v}_{i - r + j} 的 LaTeX 表示"""
    index = i_sym - r + j
    return r"\overline{v}_{" + latex(index) + r"}"

def polynomial_latex_with_cdots(use_r=True):
    """生成 a_0 + a_1 \\xi + \\cdots + a_{k-1} \\xi^{k-1} 的 LaTeX"""
    if use_r:
        term0 = r"a_{0}(r)"
        term1 = r"a_{1}(r) \xi"
        term_last = r"a_{k-1}(r) \xi^{k-1}"
    else:
        term0 = r"a_{0}"
        term1 = r"a_{1} \xi"
        term_last = r"a_{k-1} \xi^{k-1}"
    return f"{term0} + {term1} + \\cdots + {term_last}"

def integral_latex_with_cdots(use_r=True):
    """生成带省略号的积分表达式 LaTeX"""
    poly_str = polynomial_latex_with_cdots(use_r)
    alpha_str = latex(-r + j - Rational(1, 2))
    beta_str = latex(-r + j + Rational(1, 2))
    vbar_str = vbar_expr()
    return f"\\displaystyle \\int_{{{alpha_str}}}^{{{beta_str}}} ({poly_str}) \\, d\\xi = {vbar_str}"

# ---------- 构造具体 k 的表达式（用于将来矩阵生成） ----------
def build_full_polynomial(k_val, use_r=True):
    """当 k 是具体整数时，构建完整多项式（用于计算）"""
    terms = []
    for m in range(k_val):
        if use_r:
            coeff = a(m)(r)
        else:
            coeff = sp.Symbol(f'a_{m}')
        terms.append(coeff * xi**m)
    return sum(terms)

# ---------- 生成最终 LaTeX array ----------
def generate_weno_latex_array():
    lines = []

    # Line 1: p(r,xi) = a0(r) + a1(r) xi + ... + a_{k-1}(r) xi^{k-1}
    p_def = r"p(r,\xi) = " + polynomial_latex_with_cdots(use_r=True)
    lines.append(p_def)

    # Line 2: ∫ p(r,xi) dxi = vbar
    alpha_sym = r"\alpha(r,j)"
    beta_sym = r"\beta(r,j)"
    vbar_str = vbar_expr()
    lines.append(rf"\displaystyle\int_{{{alpha_sym}}}^{{{beta_sym}}} p(r,\xi) \, d\xi = {vbar_str}")

    # Line 3: 展开积分（带 a_m(r)）
    lines.append(integral_latex_with_cdots(use_r=True))

    # Line 4: 展开积分（不带 r）
    lines.append(integral_latex_with_cdots(use_r=False))

    # Line 5-6: α, β 定义
    alpha_def = latex(-r + j - Rational(1, 2))
    beta_def = latex(-r + j + Rational(1, 2))
    j_range = r"j=0,1,\ldots,k-1"
    lines.append(rf"\alpha(r,j)={alpha_def},\quad {j_range}")
    lines.append(rf"\beta(r,j)={beta_def},\quad {j_range}")

    # 组合成 array
    latex_output = "\\begin{array}{l}\n" + " \\\\\n".join(lines) + "\n\\end{array}"
    return latex_output

# ---------- 主程序 ----------
if __name__ == "__main__":
    latex_code = generate_weno_latex_array()
    print(latex_code)

    # 未来：当你需要 k=3 的具体系统时
    # p3 = build_full_polynomial(3, use_r=True)
    # print("p3 =", p3)