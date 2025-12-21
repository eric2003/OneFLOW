import sympy as sp
import re
from sympy import symbols, Rational, Matrix, latex, simplify

# ---------- 全局符号 ----------
r, i_sym = symbols('r i', integer=True)
xi = symbols(r'\xi')

# ---------- 积分限命名策略 ----------
# 提供基名（不含下标），如 \alpha, \beta
DEFAULT_BOUNDARY_NAMES = {
    'long_left': r'\alpha(r,j)',    # 用于原始公式
    'long_right': r'\beta(r,j)',
    # M 矩阵中：\alpha_j 表示第 j 个单元左边界
    'short_left_template': r'\alpha_{{{j}}}',   # 用于 M 矩阵：\alpha_0, \alpha_1, ...
    'short_right_template': r'\beta_{{{j}}}'
}

# ---------- 辅助函数 ----------
def vbar(idx_expr):
    name = r"\overline{v}_{" + latex(idx_expr) + r"}"
    return sp.Symbol(name)

def lower_limit(j): return -r + j - Rational(1, 2)
def upper_limit(j):  return -r + j + Rational(1, 2)

# ---------- 向量与矩阵构造 ----------
def build_v_vector(k_val, r_val=None):
    if r_val is None:
        return Matrix([vbar(i_sym - r + j) for j in range(k_val)])
    else:
        return Matrix([vbar(i_sym - r_val + j) for j in range(k_val)])

def build_M_matrix(k_val, r_val=None):
    M = sp.zeros(k_val, k_val)
    for j in range(k_val):
        if r_val is None:
            a_j = lower_limit(j)
            b_j = upper_limit(j)
        else:
            a_j = (-r_val + j - Rational(1, 2))
            b_j = (-r_val + j + Rational(1, 2))
        for m in range(k_val):
            integral = (b_j**(m + 1) - a_j**(m + 1)) / Rational(m + 1)
            M[j, m] = simplify(integral)
    return M

def solve_coefficients(k_val, r_val=None):
    M = build_M_matrix(k_val, r_val)
    v = build_v_vector(k_val, r_val)
    a_vec = M.inv() * v
    return a_vec, M, v

# ---------- LaTeX: 原始公式 ----------
def generate_original_latex(boundary_names=None):
    if boundary_names is None:
        boundary_names = DEFAULT_BOUNDARY_NAMES
    left = boundary_names['long_left']
    right = boundary_names['long_right']
    lines = [
        r"p(r,\xi) = a_{0}(r)+a_{1}(r)\xi +a_{2}(r)\xi^2+\cdots+a_{k-1}(r)\xi^{k-1}",
        rf"\displaystyle\int_{{{left}}}^{{{right}}} p(r,\xi) d\xi=\overline{{v}}_{{i-r+j}}",
        rf"\displaystyle \int_{{{left}}}^{{{right}}} (a_{{0}}(r)+a_{{1}}(r)\xi +\cdots+a_{{k-1}}(r)\xi^{{k-1}}) d\xi=\overline{{v}}_{{i-r+j}}",
        rf"\displaystyle \int_{{{left}}}^{{{right}}} (a_{{0}} + a_{{1}} \xi + \cdots + a_{{k-1}} \xi^{{k-1}}) d\xi=\overline{{v}}_{{i-r+j}}",
        rf"{left}=-r+j-1/2,\quad j=0,1,\cdots,k-1",
        rf"{right}=-r+j+1/2,\quad j=0,1,\cdots,k-1"
    ]
    return r"\begin{array}{l}" + "\n" + " \\\\\n".join(lines) + "\n\\end{array}"

# ---------- LaTeX: 扩展系统 ----------
def generate_extended_latexOld(boundary_names=None):
    if boundary_names is None:
        boundary_names = DEFAULT_BOUNDARY_NAMES
    base_l = boundary_names['short_base_left']   # e.g., \alpha
    base_r = boundary_names['short_base_right']  # e.g., \beta

    lines = [
        r"\mathbf{a}=[a_{0},a_{1},\dots,a_{k-1}]^T,\quad \mathbf{\phi}(\xi) = [\xi^0, \xi^1, \dots, \xi^{k-1}]^{T}",
        r"\mathbf{v}=[\overline{v}_{i - r},\overline{v}_{i - r + 1},\dots,\overline{v}_{i + k - r - 1}]^T",
        r"\mathbf{a}=M^{-1}\mathbf{v}",
        r"\begin{bmatrix} a_0 \\ a_1 \\ \vdots \\ a_{k-1} \end{bmatrix} = M^{-1} \begin{bmatrix} \overline{v}_{i - r} \\ \overline{v}_{i - r + 1} \\ \vdots \\ \overline{v}_{i + k - r - 1} \end{bmatrix}"
    ]
    
    # 构建 M 矩阵 LaTeX：使用 \alpha_0, \alpha_1, \alpha_{k-1}
    def make_row(j_str):
        return " & ".join([
            rf"\int_{{{base_l}_{{{j_str}}}}}^{{{base_r}_{{{j_str}}}}} d\xi",
            rf"\int_{{{base_l}_{{{j_str}}}}}^{{{base_r}_{{{j_str}}}}} \xi^{{1}} d\xi",
            r"\cdots",
            rf"\int_{{{base_l}_{{{j_str}}}}}^{{{base_r}_{{{j_str}}}}} \xi^{{k-1}} d\xi"
        ])
    
    M_body = " \\\\\n".join([
        make_row("0"),
        make_row("1"),
        r"\vdots & \vdots & \ddots & \vdots",
        make_row("k-1")
    ])
    M_latex = f"M = \\begin{{bmatrix}}\n{M_body}\n\\end{{bmatrix}}"
    lines.append(M_latex)
    return r"\begin{array}{l}" + "\n" + " \\\\\n".join(lines) + "\n\\end{array}"
    
def generate_extended_latex(boundary_names=None):
    if boundary_names is None:
        boundary_names = DEFAULT_BOUNDARY_NAMES
    left_template = boundary_names['short_left_template']   # e.g., r'\alpha_{{{j}}}'
    right_template = boundary_names['short_right_template'] # e.g., r'\beta_{{{j}}}'

    lines = [
        r"\mathbf{a}=[a_{0},a_{1},\dots,a_{k-1}]^T,\quad \mathbf{\phi}(\xi) = [\xi^0, \xi^1, \dots, \xi^{k-1}]^{T}",
        r"\mathbf{v}=[\overline{v}_{i - r},\overline{v}_{i - r + 1},\dots,\overline{v}_{i + k - r - 1}]^T",
        r"\mathbf{a}=M^{-1}\mathbf{v}",
        r"\begin{bmatrix} a_0 \\ a_1 \\ \vdots \\ a_{k-1} \end{bmatrix} = M^{-1} \begin{bmatrix} \overline{v}_{i - r} \\ \overline{v}_{i - r + 1} \\ \vdots \\ \overline{v}_{i + k - r - 1} \end{bmatrix}"
    ]
    
    def make_row(j_str):
        left = left_template.format(j=j_str)
        right = right_template.format(j=j_str)
        return " & ".join([
            rf"\int_{{{left}}}^{{{right}}} d\xi",
            rf"\int_{{{left}}}^{{{right}}} \xi^{{1}} d\xi",
            r"\cdots",
            rf"\int_{{{left}}}^{{{right}}} \xi^{{k-1}} d\xi"
        ])
    
    M_body = " \\\\\n".join([
        make_row("0"),
        make_row("1"),
        r"\vdots & \vdots & \ddots & \vdots",
        make_row("k-1")
    ])
    M_latex = f"M = \\begin{{bmatrix}}\n{M_body}\n\\end{{bmatrix}}"
    lines.append(M_latex)
    return r"\begin{array}{l}" + "\n" + " \\\\\n".join(lines) + "\n\\end{array}"    

# ---------- 排序与格式化（无变动） ----------
def _is_vbar_symbol(sym):
    return isinstance(sym, sp.Symbol) and sym.name.startswith(r"\overline{v}_{")

def _extract_index_from_vbar(v_symbol):
    match = re.search(r'\\overline\{v\}_\{(.+)\}', v_symbol.name)
    if not match:
        return v_symbol.name
    index_latex = match.group(1)
    try:
        return simplify(eval(index_latex, {'i': i_sym, 'r': r}))
    except:
        return index_latex

def _sort_key_from_index(index_expr):
    if isinstance(index_expr, str):
        return (1, index_expr)
    try:
        offset = simplify(index_expr - i_sym)
        if offset.is_number:
            return (0, float(offset))
        else:
            return (0.5, str(offset))
    except:
        return (1, str(index_expr))

def format_a_expression(expr, factored=True):
    if expr.is_Number:
        return latex(expr)
    expr = simplify(expr)
    has_vbar = any(_is_vbar_symbol(s) for s in expr.free_symbols)
    if not has_vbar:
        return latex(expr)

    terms = expr.as_ordered_terms() if expr.is_Add else [expr]
    parsed_terms = []
    for term in terms:
        v_syms = [s for s in term.free_symbols if _is_vbar_symbol(s)]
        if not v_syms:
            parsed_terms.append((term, None, None))
        else:
            if len(v_syms) != 1:
                raise ValueError(f"Unexpected term: {term}")
            v = v_syms[0]
            coeff = simplify(term / v)
            index_expr = _extract_index_from_vbar(v)
            parsed_terms.append((coeff, v, index_expr))
    
    def term_sort_key(item):
        _, v, idx = item
        if v is None:
            return (-1, 0)
        return _sort_key_from_index(idx)
    
    parsed_terms.sort(key=term_sort_key)

    latex_parts = []
    for coeff, v, _ in parsed_terms:
        if v is None:
            latex_parts.append(latex(coeff))
        else:
            if factored:
                if coeff == 1:
                    s = latex(v)
                elif coeff == -1:
                    s = "-" + latex(v)
                else:
                    c_latex = latex(coeff)
                    if any(op in c_latex for op in ['+', '-']) and not c_latex.startswith('-'):
                        c_latex = f"({c_latex})"
                    s = f"{c_latex} {latex(v)}"
            else:
                if coeff.is_Rational and abs(coeff.p) == 1:
                    sign = "-" if coeff.p < 0 else ""
                    den = str(coeff.q)
                    s = f"{sign}\\frac{{{latex(v)}}}{{{den}}}"
                else:
                    s = latex(coeff * v)
            latex_parts.append(s)
    
    result = latex_parts[0]
    for part in latex_parts[1:]:
        if part.startswith("-"):
            result += " " + part
        else:
            result += " + " + part
    return result

def latex_a_vector_formatted(a_vec, factored=True):
    entries = [format_a_expression(a_vec[i], factored=factored) for i in range(len(a_vec))]
    body = " \\\\ ".join(entries)
    return f"\\begin{{bmatrix}} {body} \\end{{bmatrix}}"

# ---------- 主输出函数 ----------
def print_solution_for_k_r(k_val, r_val=None, factored=True):
    a_vec, M, v = solve_coefficients(k_val, r_val)
    desc = f"k={k_val}" + (f", r={r_val}" if r_val is not None else " (r symbolic)")
    print(f"\n=== WENO Coefficients: {desc} ===")
    print(latex_a_vector_formatted(a_vec, factored=factored))
    return a_vec, M, v

# ---------- 示例：不同命名策略 ----------
X_BOUNDARY_NAMES = {
    'long_left': r'x_{j-\frac{1}{2}}',
    'long_right': r'x_{j+\frac{1}{2}}',
    # 对 M 矩阵，提供 j 的占位符，如 {j} - 1/2
    'short_left_template': r'x_{{{j}-\frac{{1}}{{2}}}}',
    'short_right_template': r'x_{{{j}+\frac{{1}}{{2}}}}'
}

# ---------- 主程序 ----------
if __name__ == "__main__":
    # 1. 默认命名 (\alpha, \beta)
    print("=== Original LaTeX (default: \\alpha, \\beta) ===")
    print(generate_original_latex())
    
    print("\n\n=== Extended System (default) ===")
    print(generate_extended_latex())
    
    # 2. 使用 x_{j±1/2} 命名
    print("\n\n=== Original LaTeX (x_{j±1/2}) ===")
    print(generate_original_latex(boundary_names=X_BOUNDARY_NAMES))
    
    print("\n\n=== Extended System (x_{j±1/2}) ===")
    print(generate_extended_latex(boundary_names=X_BOUNDARY_NAMES))
    
    # 3. 系数输出 (不受命名影响)
    for r_val in [0, 1, 2]:
        print_solution_for_k_r(k_val=3, r_val=r_val, factored=True)