import sympy as sp
import re
from sympy import symbols, Rational, Matrix, latex, simplify

# ---------- Global symbols ----------
r, i_sym = symbols('r i', integer=True)
xi = symbols(r'\xi')

# ---------- Boundary naming strategy ----------
# Template-based naming to avoid double subscripts and allow flexible notation
DEFAULT_BOUNDARY_NAMES = {
    'long_left': r'\alpha(r,j)',    # Used in original formulas
    'long_right': r'\beta(r,j)',
    'short_left_template': r'\alpha_{{{j}}}',   # Template for M matrix: \alpha_{j}
    'short_right_template': r'\beta_{{{j}}}'
}

# Alternative naming for cell-centered finite volume methods
X_BOUNDARY_NAMES = {
    'long_left': r'x_{j-\frac{1}{2}}',
    'long_right': r'x_{j+\frac{1}{2}}',
    'short_left_template': r'x_{{{j}-\frac{{1}}{{2}}}}',
    'short_right_template': r'x_{{{j}+\frac{{1}}{{2}}}}'
}

# ---------- Helper functions ----------
def vbar(idx_expr):
    """Create a Symbol with LaTeX name \overline{v}_{idx_expr} for proper rendering."""
    name = r"\overline{v}_{" + latex(idx_expr) + r"}"
    return sp.Symbol(name)

def lower_limit(j): 
    """Lower integration limit: α(r,j) = -r + j - 1/2"""
    return -r + j - Rational(1, 2)

def upper_limit(j):  
    """Upper integration limit: β(r,j) = -r + j + 1/2"""
    return -r + j + Rational(1, 2)

# ---------- Vector and matrix construction ----------
def build_v_vector(k_val, r_val=None):
    """Construct the right-hand side vector v = [v̄_{i-r+j}]_{j=0}^{k-1}."""
    if r_val is None:
        return Matrix([vbar(i_sym - r + j) for j in range(k_val)])
    else:
        return Matrix([vbar(i_sym - r_val + j) for j in range(k_val)])

def build_M_matrix(k_val, r_val=None):
    """
    Construct the moment matrix M where M[j,m] = ∫_{α_j}^{β_j} ξ^m dξ.
    
    Parameters:
        k_val (int): Polynomial order (degree k-1)
        r_val (int, optional): Specific r value. If None, keep r symbolic.
    """
    M = sp.zeros(k_val, k_val)
    for j in range(k_val):
        if r_val is not None:
            # Substitute specific r value into symbolic limits
            a_j = lower_limit(j).subs(r, r_val)
            b_j = upper_limit(j).subs(r, r_val)
        else:
            a_j = lower_limit(j)
            b_j = upper_limit(j)
        
        for m in range(k_val):
            # Analytical integration: ∫ ξ^m dξ = (β^{m+1} - α^{m+1}) / (m+1)
            integral = (b_j**(m + 1) - a_j**(m + 1)) / Rational(m + 1)
            M[j, m] = simplify(integral)
    return M

def solve_coefficients(k_val, r_val=None):
    """Solve a = M^{-1} v for polynomial coefficients."""
    M = build_M_matrix(k_val, r_val)
    v = build_v_vector(k_val, r_val)
    a_vec = M.inv() * v
    return a_vec, M, v

# ---------- LaTeX generation: Original formulas ----------
def generate_original_latex(boundary_names=None):
    """
    Generate the original LaTeX block with integration constraints.
    
    Uses boundary_names['long_left'] and boundary_names['long_right'] 
    for integration limits.
    """
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

# ---------- LaTeX generation: Extended system ----------
def generate_extended_latex(boundary_names=None):
    """
    Generate the extended system LaTeX including a = M^{-1}v and matrix M.
    
    Uses template-based naming for M matrix entries to avoid double subscripts.
    Templates should contain {j} placeholder for row index.
    """
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
        """Generate a matrix row with proper boundary notation."""
        # Safe template substitution
        left = left_template.format(j=j_str)
        right = right_template.format(j=j_str)
        return " & ".join([
            rf"\int_{{{left}}}^{{{right}}} d\xi",
            rf"\int_{{{left}}}^{{{right}}} \xi^{{1}} d\xi",
            r"\cdots",
            rf"\int_{{{left}}}^{{{right}}} \xi^{{k-1}} d\xi"
        ])
    
    # Construct matrix with representative rows (0, 1, ..., k-1)
    M_body = " \\\\\n".join([
        make_row("0"),
        make_row("1"),
        r"\vdots & \vdots & \ddots & \vdots",
        make_row("k-1")
    ])
    M_latex = f"M = \\begin{{bmatrix}}\n{M_body}\n\\end{{bmatrix}}"
    lines.append(M_latex)
    return r"\begin{array}{l}" + "\n" + " \\\\\n".join(lines) + "\n\\end{array}"

# ---------- Expression formatting with sorted vbar terms ----------
def _is_vbar_symbol(sym):
    """Check if symbol represents \overline{v}_{...}."""
    return isinstance(sym, sp.Symbol) and sym.name.startswith(r"\overline{v}_{")

def _extract_index_from_vbar(v_symbol):
    """Extract the index expression from \overline{v}_{index} symbol name."""
    match = re.search(r'\\overline\{v\}_\{(.+)\}', v_symbol.name)
    if not match:
        return v_symbol.name
    index_latex = match.group(1)
    try:
        # Evaluate in context of global symbols i and r
        return simplify(eval(index_latex, {'i': i_sym, 'r': r}))
    except:
        return index_latex

def _sort_key_from_index(index_expr):
    """
    Generate sort key for vbar indices.
    Priority: numeric offsets from i (i-1, i, i+1, ...) > symbolic > strings.
    """
    if isinstance(index_expr, str):
        return (1, index_expr)
    try:
        offset = simplify(index_expr - i_sym)
        if offset.is_number:
            return (0, float(offset))  # Sort by numeric offset
        else:
            return (0.5, str(offset))  # Symbolic offsets
    except:
        return (1, str(index_expr))

def format_a_expression(expr, factored=True):
    """
    Format polynomial coefficient expression with vbar terms sorted by index.
    
    Parameters:
        expr: SymPy expression for a coefficient
        factored: If True, output as (1/24)*v̄; if False, output as v̄/24
    """
    if expr.is_Number:
        return latex(expr)
    expr = simplify(expr)
    
    # Check if expression contains any vbar symbols
    has_vbar = any(_is_vbar_symbol(s) for s in expr.free_symbols)
    if not has_vbar:
        return latex(expr)

    # Parse terms into (coefficient, vbar_symbol, index) tuples
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
    
    # Sort terms: vbar terms by index offset, constants first
    def term_sort_key(item):
        _, v, idx = item
        if v is None:
            return (-1, 0)  # Constants first
        return _sort_key_from_index(idx)
    
    parsed_terms.sort(key=term_sort_key)

    # Generate LaTeX for each term
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
                    # Add parentheses for compound coefficients
                    if any(op in c_latex for op in ['+', '-']) and not c_latex.startswith('-'):
                        c_latex = f"({c_latex})"
                    s = f"{c_latex} {latex(v)}"
            else:
                # Fraction format: v̄/denominator
                if coeff.is_Rational and abs(coeff.p) == 1:
                    sign = "-" if coeff.p < 0 else ""
                    den = str(coeff.q)
                    s = f"{sign}\\frac{{{latex(v)}}}{{{den}}}"
                else:
                    s = latex(coeff * v)
            latex_parts.append(s)
    
    # Combine terms with proper sign handling
    result = latex_parts[0]
    for part in latex_parts[1:]:
        if part.startswith("-"):
            result += " " + part
        else:
            result += " + " + part
    return result

def latex_a_vector_formatted(a_vec, factored=True):
    """Format coefficient vector as bmatrix with sorted terms."""
    entries = [format_a_expression(a_vec[i], factored=factored) for i in range(len(a_vec))]
    body = " \\\\ ".join(entries)
    return f"\\begin{{bmatrix}} {body} \\end{{bmatrix}}"

# ---------- Main output functions ----------
def print_solution_for_k_r(k_val, r_val=None, factored=True):
    """Print formatted coefficient vector for given k and r values."""
    a_vec, M, v = solve_coefficients(k_val, r_val)
    desc = f"k={k_val}" + (f", r={r_val}" if r_val is not None else " (r symbolic)")
    print(f"\n=== WENO Coefficients: {desc} ===")
    print(latex_a_vector_formatted(a_vec, factored=factored))
    return a_vec, M, v

# ---------- Main execution ----------
if __name__ == "__main__":
    # 1. Default naming (\alpha, \beta)
    print("=== Original LaTeX (default: \\alpha, \\beta) ===")
    print(generate_original_latex())
    
    print("\n\n=== Extended System (default) ===")
    print(generate_extended_latex())
    
    # 2. Cell-centered naming (x_{j±1/2})
    print("\n\n=== Original LaTeX (x_{j±1/2}) ===")
    print(generate_original_latex(boundary_names=X_BOUNDARY_NAMES))
    
    print("\n\n=== Extended System (x_{j±1/2}) ===")
    print(generate_extended_latex(boundary_names=X_BOUNDARY_NAMES))
    
    # 3. Coefficient solutions for k=3, r=0,1,2
    for r_val in [0, 1, 2]:
        print_solution_for_k_r(k_val=3, r_val=r_val, factored=True)