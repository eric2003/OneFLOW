import sympy as sp
import re
from sympy import symbols, Rational, Matrix, latex, simplify

# ---------- Global symbols ----------
r, i_sym = symbols('r i', integer=True)
xi = symbols(r'\xi')

# ---------- Boundary naming strategy ----------
DEFAULT_BOUNDARY_NAMES = {
    'long_left': r'\alpha(r,j)',    
    'long_right': r'\beta(r,j)',
    'short_left_template': r'\alpha_{{{j}}}',   
    'short_right_template': r'\beta_{{{j}}}'
}

X_BOUNDARY_NAMES = {
    'long_left': r'x_{j-\frac{1}{2}}',
    'long_right': r'x_{j+\frac{1}{2}}',
    'short_left_template': r'x_{{{j}-\frac{{1}}{{2}}}}',
    'short_right_template': r'x_{{{j}+\frac{{1}}{{2}}}}'
}

# ---------- Helper functions ----------
def vbar(idx_expr):
    name = r"\overline{v}_{" + latex(idx_expr) + r"}"
    return sp.Symbol(name)

def lower_limit(j): 
    return -r + j - Rational(1, 2)

def upper_limit(j):  
    return -r + j + Rational(1, 2)

# ---------- Vector and matrix construction ----------
def build_v_vector(k_val, r_val=None):
    if r_val is None:
        return Matrix([vbar(i_sym - r + j) for j in range(k_val)])
    else:
        return Matrix([vbar(i_sym - r_val + j) for j in range(k_val)])

def build_M_matrix(k_val, r_val=None):
    M = sp.zeros(k_val, k_val)
    for j in range(k_val):
        if r_val is not None:
            a_j = lower_limit(j).subs(r, r_val)
            b_j = upper_limit(j).subs(r, r_val)
        else:
            a_j = lower_limit(j)
            b_j = upper_limit(j)
        
        for m in range(k_val):
            integral = (b_j**(m + 1) - a_j**(m + 1)) / Rational(m + 1)
            M[j, m] = simplify(integral)
    return M

def solve_coefficients(k_val, r_val=None):
    M = build_M_matrix(k_val, r_val)
    v = build_v_vector(k_val, r_val)
    a_vec = M.inv() * v
    return a_vec, M, v

# ---------- NEW: Smoothness indicator Sr computation ----------
def compute_Sr(k_val, r_val=None):
    """
    Compute the WENO smoothness indicator Sr = sum_{l=1}^{k-1} ∫_{-1/2}^{1/2} (d^l p/dξ^l)^2 dξ
    
    Parameters:
        k_val (int): Polynomial order (degree k-1)
        r_val (int, optional): Specific r value. If None, keep r symbolic.
    
    Returns:
        Sr_expr: Symbolic expression for Sr in terms of vbar values
        p_expr: The polynomial p(r,ξ)
        a_vec: Coefficient vector [a0, a1, ..., a_{k-1}]
    """
    # Step 1: Get coefficients a = [a0, a1, ..., a_{k-1}]
    a_vec, M, v = solve_coefficients(k_val, r_val)
    
    # Step 2: Construct polynomial p(r,ξ) = sum_{m=0}^{k-1} a_m * ξ^m
    p_expr = sum(a_vec[m] * xi**m for m in range(k_val))
    
    # Step 3: Compute Sr = sum_{l=1}^{k-1} ∫_{-1/2}^{1/2} (p^{(l)}(ξ))^2 dξ
    Sr_expr = 0
    for l in range(1, k_val):  # l = 1, 2, ..., k-1
        # Compute l-th derivative
        p_l_deriv = sp.diff(p_expr, xi, l)
        # Square it
        p_l_squared = p_l_deriv**2
        # Integrate over [-1/2, 1/2]
        integral_l = sp.integrate(p_l_squared, (xi, -Rational(1, 2), Rational(1, 2)))
        Sr_expr += simplify(integral_l)
    
    # Final simplification
    Sr_expr = simplify(Sr_expr)
    return Sr_expr, p_expr, a_vec

def format_Sr_expression(Sr_expr, factored=True):
    """
    Format Sr expression with vbar terms sorted by index (similar to coefficient formatting).
    """
    if Sr_expr.is_Number:
        return latex(Sr_expr)
    Sr_expr = simplify(Sr_expr)
    
    # Reuse the existing vbar detection and sorting logic
    has_vbar = any(_is_vbar_symbol(s) for s in Sr_expr.free_symbols)
    if not has_vbar:
        return latex(Sr_expr)

    terms = Sr_expr.as_ordered_terms() if Sr_expr.is_Add else [Sr_expr]
    parsed_terms = []
    for term in terms:
        v_syms = [s for s in term.free_symbols if _is_vbar_symbol(s)]
        if not v_syms:
            parsed_terms.append((term, None, None))
        else:
            # Sr may have products of vbar terms (like v_i * v_{i+1})
            # For now, handle single vbar terms (which is the case for k=3)
            # For higher k, we might have cross terms
            if len(v_syms) == 1:
                v = v_syms[0]
                coeff = simplify(term / v)
                index_expr = _extract_index_from_vbar(v)
                parsed_terms.append((coeff, v, index_expr))
            else:
                # Handle cross terms or complex expressions
                parsed_terms.append((term, None, None))
    
    # Sort single vbar terms, keep others as-is
    single_vbar_terms = []
    other_terms = []
    for item in parsed_terms:
        if item[1] is not None:
            single_vbar_terms.append(item)
        else:
            other_terms.append(item)
    
    # Sort single vbar terms
    def term_sort_key(item):
        _, v, idx = item
        if v is None:
            return (-1, 0)
        return _sort_key_from_index(idx)
    
    single_vbar_terms.sort(key=term_sort_key)
    all_terms = other_terms + single_vbar_terms

    # Generate LaTeX (simplified for Sr - usually no factoring needed)
    latex_parts = []
    for coeff, v, _ in all_terms:
        if v is None:
            latex_parts.append(latex(coeff))
        else:
            # For Sr, we typically want expanded form
            latex_parts.append(latex(coeff * v))
    
    # Combine terms
    if not latex_parts:
        return "0"
    result = latex_parts[0]
    for part in latex_parts[1:]:
        if part.startswith("-"):
            result += " " + part
        else:
            result += " + " + part
    return result

# ---------- Existing functions (unchanged) ----------
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

def generate_extended_latex(boundary_names=None):
    if boundary_names is None:
        boundary_names = DEFAULT_BOUNDARY_NAMES
    left_template = boundary_names['short_left_template']
    right_template = boundary_names['short_right_template']

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

def print_solution_for_k_r(k_val, r_val=None, factored=True):
    a_vec, M, v = solve_coefficients(k_val, r_val)
    desc = f"k={k_val}" + (f", r={r_val}" if r_val is not None else " (r symbolic)")
    print(f"\n=== WENO Coefficients: {desc} ===")
    print(latex_a_vector_formatted(a_vec, factored=factored))
    return a_vec, M, v

# ---------- Main execution with Sr computation ----------
if __name__ == "__main__":
    k_val = 3
    
    # Print original and extended systems
    print("=== Original LaTeX (default) ===")
    print(generate_original_latex())
    
    print("\n\n=== Extended System (default) ===")
    print(generate_extended_latex())
    
    # Compute coefficients and Sr for r = 0, 1, 2
    for r_val in [0, 1, 2]:
        print(f"\n" + "="*60)
        print(f"=== k={k_val}, r={r_val} ===")
        
        # Coefficients
        a_vec, M, v = solve_coefficients(k_val, r_val)
        print("\nCoefficients a = M^{-1} v:")
        print(latex_a_vector_formatted(a_vec, factored=True))
        
        # Smoothness indicator Sr
        Sr_expr, p_expr, a_vec_sr = compute_Sr(k_val, r_val)
        print(f"\nSmoothness indicator S_{r_val}:")
        print(format_Sr_expression(Sr_expr))
        
        # Optional: Show the polynomial
        print(f"\nPolynomial p({r_val},ξ):")
        print(latex(p_expr))
    
    # Also show symbolic Sr (r as symbol) for completeness
    print(f"\n" + "="*60)
    print(f"=== Symbolic Sr (k={k_val}, r symbolic) ===")
    Sr_sym, p_sym, a_sym = compute_Sr(k_val, r_val=None)
    print("S_r =")
    print(format_Sr_expression(Sr_sym))