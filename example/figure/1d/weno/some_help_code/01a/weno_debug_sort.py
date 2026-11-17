import sympy as sp
import re
from sympy import symbols, Rational, Matrix, latex, simplify

# ---------- 全局符号 ----------
r, i_sym = symbols('r i', integer=True)
xi = symbols(r'\xi')

# ---------- vbar 构造 ----------
def vbar(idx_expr):
    name = r"\overline{v}_{" + latex(idx_expr) + r"}"
    return sp.Symbol(name)

# ---------- 积分限 ----------
def alpha(j): return -r + j - Rational(1, 2)
def beta(j):  return -r + j + Rational(1, 2)

# ---------- 矩阵与向量 ----------
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

# ---------- 排序辅助函数 ----------
def _is_vbar_symbol(sym):
    return isinstance(sym, sp.Symbol) and sym.name.startswith(r"\overline{v}_{")

def _extract_index_from_vbar(v_symbol):
    match = re.search(r'\\overline\{v\}_\{(.+)\}', v_symbol.name)
    if not match:
        return v_symbol.name
    index_latex = match.group(1)
    try:
        # 使用全局 i_sym, r
        index_expr = eval(index_latex, {'i': i_sym, 'r': r})
        return simplify(index_expr)
    except Exception as e:
        print(f"[EXTRACT FAILED] name={v_symbol.name}, index_latex={index_latex}, error={e}")
        return index_latex

def _sort_key_from_index(index_expr):
    if isinstance(index_expr, str):
        print(f"  [SORT KEY] string fallback: {index_expr}")
        return (1, index_expr)
    try:
        offset = simplify(index_expr - i_sym)
        if offset.is_number:
            key = (0, float(offset))
            print(f"  [SORT KEY] {index_expr} → offset={offset} → key={key}")
            return key
        else:
            key = (0.5, str(offset))
            print(f"  [SORT KEY] non-numeric offset: {offset} → key={key}")
            return key
    except Exception as e:
        key = (1, str(index_expr))
        print(f"  [SORT KEY] exception: {e} → key={key}")
        return key

# ---------- 格式化表达式（带详细调试） ----------
def format_a_expression(expr, factored=True):
    if expr.is_Number:
        return latex(expr)
    expr = simplify(expr)
    if not expr.has(_is_vbar_symbol):
        return latex(expr)

    terms = expr.as_ordered_terms() if expr.is_Add else [expr]
    parsed_terms = []
    for term in terms:
        v_syms = [s for s in term.free_symbols if _is_vbar_symbol(s)]
        if not v_syms:
            parsed_terms.append((term, None, None))
        else:
            if len(v_syms) != 1:
                raise ValueError(f"Term has multiple vbar: {term}")
            v = v_syms[0]
            coeff = simplify(term / v)
            index_expr = _extract_index_from_vbar(v)
            parsed_terms.append((coeff, v, index_expr))
    
    # === DEBUG: 打印解析结果 ===
    print(f"\n[DEBUG] Expression: {expr}")
    for coeff, v, idx in parsed_terms:
        if v is None:
            print(f"  CONSTANT: {coeff}")
        else:
            print(f"  TERM: coeff={coeff}, v={v.name}, index_expr={idx} (type={type(idx)})")

    # === 排序 ===
    def term_sort_key(item):
        _, v, idx = item
        if v is None:
            return (-1, 0)  # 常数项最先
        return _sort_key_from_index(idx)

    print("[DEBUG] Sorting terms...")
    parsed_terms.sort(key=term_sort_key)

    # === 生成 LaTeX ===
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
    
    # 合并
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

# ---------- 主程序 ----------
if __name__ == "__main__":
    print("=== Testing k=3, r=0 ===")
    a_vec, M, v = solve_coefficients(k_val=3, r_val=0)
    print("\nRaw a_vec from SymPy:")
    for i, ai in enumerate(a_vec):
        print(f"a{i} = {ai}")

    print("\nFormatted output (factored=True):")
    formatted = latex_a_vector_formatted(a_vec, factored=True)
    print(formatted)