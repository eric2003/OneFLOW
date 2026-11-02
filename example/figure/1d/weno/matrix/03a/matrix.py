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

# ---------- Enhanced: Smoothness indicator Sr computation ----------
def compute_Sr(k_val, r_val=None):
    """
    Compute the WENO smoothness indicator Sr = sum_{l=1}^{k-1} ∫_{-1/2}^{1/2} (d^l p/dξ^l)^2 dξ
    
    This works for any k >= 2 and handles cross terms correctly.
    """
    # Get coefficients a = [a0, a1, ..., a_{k-1}]
    a_vec, M, v = solve_coefficients(k_val, r_val)
    
    # Construct polynomial p(r,ξ) = sum_{m=0}^{k-1} a_m * ξ^m
    p_expr = sum(a_vec[m] * xi**m for m in range(k_val))
    
    # Compute Sr = sum_{l=1}^{k-1} ∫_{-1/2}^{1/2} (p^{(l)}(ξ))^2 dξ
    Sr_expr = 0
    for l in range(1, k_val):  # l = 1, 2, ..., k-1
        p_l_deriv = sp.diff(p_expr, xi, l)
        p_l_squared = p_l_deriv**2
        integral_l = sp.integrate(p_l_squared, (xi, -Rational(1, 2), Rational(1, 2)))
        Sr_expr += simplify(integral_l)
    
    Sr_expr = simplify(Sr_expr)
    return Sr_expr, p_expr, a_vec

def format_Sr_expression(Sr_expr):
    """
    Format Sr expression properly, handling both single terms and cross terms.
    
    For WENO smoothness indicators, the standard format is expanded form
    showing all quadratic terms explicitly.
    """
    if Sr_expr.is_Number:
        return latex(Sr_expr)
    
    Sr_expr = simplify(Sr_expr)
    
    # For Sr, we want to expand everything to show the quadratic form
    # This is the standard way smoothness indicators are presented in literature
    expanded_expr = sp.expand(Sr_expr)
    
    # If it's just a number after expansion
    if expanded_expr.is_Number:
        return latex(expanded_expr)
    
    # Handle the general case by converting directly to LaTeX
    # The expanded form will show all cross terms properly
    return latex(expanded_expr)

# ---------- Utility function to get vbar symbols from expression ----------
def get_vbar_symbols(expr):
    """Extract all vbar symbols from an expression."""
    return [s for s in expr.free_symbols if _is_vbar_symbol(s)]

# ---------- Existing helper functions ----------
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

# ---------- Main execution for both k=3 and k=4 ----------
if __name__ == "__main__":
    # Test both k=3 and k=4
    for k_val in [3, 4]:
        print(f"\n{'='*80}")
        print(f"COMPUTING FOR k = {k_val}")
        print(f"{'='*80}")
        
        print("\n=== Original LaTeX ===")
        print(generate_original_latex())
        
        # Compute for r = 0, 1, ..., k-1
        for r_val in range(k_val):
            print(f"\n" + "-"*60)
            print(f"r = {r_val} (k = {k_val})")
            
            # Coefficients
            a_vec, M, v = solve_coefficients(k_val, r_val)
            print("\nCoefficients:")
            print(latex_a_vector_formatted(a_vec, factored=True))
            
            # Smoothness indicator
            Sr_expr, p_expr, _ = compute_Sr(k_val, r_val)
            print(f"\nSmoothness indicator S_{r_val}:")
            print(format_Sr_expression(Sr_expr))
            
            # Show the stencil size
            v_symbols = get_vbar_symbols(Sr_expr)
            if v_symbols:
                print(f"\nStencil involves {len(v_symbols)} points: {[s.name for s in v_symbols]}")
    
    # Show symbolic version for k=4
    print(f"\n{'='*80}")
    print("SYMBOLIC S_r FOR k=4 (r as symbol)")
    print(f"{'='*80}")
    Sr_sym, p_sym, _ = compute_Sr(4, r_val=None)
    print("S_r =")
    print(format_Sr_expression(Sr_sym))