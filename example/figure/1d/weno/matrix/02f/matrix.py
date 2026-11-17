import sympy as sp
import re
from sympy import symbols, Rational, Matrix, latex, simplify

# ---------- 全局符号 ----------
r, i_sym = symbols('r i', integer=True)
xi = symbols(r'\xi')

# ---------- 辅助函数 ----------
def vbar(idx_expr):
    name = r"\overline{v}_{" + latex(idx_expr) + r"}"
    return sp.Symbol(name)

def alpha(j): return -r + j - Rational(1, 2)
def beta(j):  return -r + j + Rational(1, 2)

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
            a_j = alpha(j)
            b_j = beta(j)
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
def generate_original_latex():
    lines = [
        r"p(r,\xi) = a_{0}(r)+a_{1}(r)\xi +a_{2}(r)\xi^2+\cdots+a_{k-1}(r)\xi^{k-1}",
        r"\displaystyle\int_{\alpha(r,j)}^{\beta(r,j)} p(r,\xi) d\xi=\overline{v}_{i-r+j}",
        r"\displaystyle \int_{\alpha(r,j)}^{\beta(r,j)} (a_{0}(r)+a_{1}(r)\xi +a_{2}(r)\xi^2+\cdots+a_{k-1}(r)\xi^{k-1}) d\xi=\overline{v}_{i-r+j}",
        r"\displaystyle \int_{\alpha(r,j)}^{\beta(r,j)} (a_{0} + a_{1} \xi + a_{2} \xi^2 + \cdots + a_{k-1} \xi^{k-1}) d\xi=\overline{v}_{i-r+j}",
        r"\alpha(r,j)=-r+j-1/2,\quad j=0,1,\cdots,k-1",
        r"\beta(r,j)=-r+j+1/2,\quad j=0,1,\cdots,k-1"
    ]
    return r"\begin{array}{l}" + "\n" + " \\\\\n".join(lines) + "\n\\end{array}"

# ---------- LaTeX: 扩展系统 ----------
def generate_extended_latex():
    lines = [
        r"\mathbf{a}=[a_{0},a_{1},\dots,a_{k-1}]^T,\quad \mathbf{\phi}(\xi) = [\xi^0, \xi^1, \dots, \xi^{k-1}]^{T}",
        r"\mathbf{v}=[\overline{v}_{i - r},\overline{v}_{i - r + 1},\dots,\overline{v}_{i + k - r - 1}]^T",
        r"\mathbf{a}=M^{-1}\mathbf{v}",
        r"\begin{bmatrix} a_0 \\ a_1 \\ \vdots \\ a_{k-1} \end{bmatrix} = M^{-1} \begin{bmatrix} \overline{v}_{i - r} \\ \overline{v}_{i - r + 1} \\ \vdots \\ \overline{v}_{i + k - r - 1} \end{bmatrix}"
    ]
    M_latex = r"M = \begin{bmatrix}" + "\n"
    M_latex += r"\int_{\alpha_{0}}^{\beta_{0}}d\xi & \int_{\alpha_{0}}^{\beta_{0}}{\xi}^{1}d\xi & \cdots & \int_{\alpha_{0}}^{\beta_{0}}{\xi}^{k-1}d\xi \\" + "\n"
    M_latex += r"\int_{\alpha_{1}}^{\beta_{1}}d\xi & \int_{\alpha_{1}}^{\beta_{1}}{\xi}^{1}d\xi & \cdots & \int_{\alpha_{1}}^{\beta_{1}}{\xi}^{k-1}d\xi \\" + "\n"
    M_latex += r"\vdots & \vdots & \ddots & \vdots \\" + "\n"
    M_latex += r"\int_{\alpha_{k-1}}^{\beta_{k-1}}d\xi & \int_{\alpha_{k-1}}^{\beta_{k-1}}{\xi}^{1}d\xi & \cdots & \int_{\alpha_{k-1}}^{\beta_{k-1}}{\xi}^{k-1}d\xi" + "\n"
    M_latex += r"\end{bmatrix}"
    lines.append(M_latex)
    return r"\begin{array}{l}" + "\n" + " \\\\\n".join(lines) + "\n\\end{array}"

# ---------- 排序与格式化 ----------
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
    
    # 正确检测是否包含 vbar 符号
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
    
    # 排序：vbar 项按 index 排序，常数项放前（WENO 中通常无常数项）
    def term_sort_key(item):
        _, v, idx = item
        if v is None:
            return (-1, 0)
        return _sort_key_from_index(idx)
    
    parsed_terms.sort(key=term_sort_key)

    # 生成 LaTeX
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

# ---------- 主程序 ----------
if __name__ == "__main__":
    # 1. 原始公式
    print("=== Original LaTeX Block ===")
    print(generate_original_latex())
    
    # 2. 扩展系统
    print("\n\n=== Extended System LaTeX ===")
    print(generate_extended_latex())
    
    # 3. 具体系数 (k=3, r=0,1,2)
    for r_val in [0, 1, 2]:
        print_solution_for_k_r(k_val=3, r_val=r_val, factored=True)