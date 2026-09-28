import sympy as sp

# 作为符号处理（避免 Function 问题）
v_i, v_ip1, v_ip2, v_ip3 = sp.symbols('v_i v_{i+1} v_{i+2} v_{i+3}')

beta0_expanded = (
    sp.Rational(2107, 240) * v_i**2 - sp.Rational(1567, 40) * v_i * v_ip1 + sp.Rational(3521, 120) * v_i * v_ip2 - sp.Rational(309, 40) * v_i * v_ip3 +
    sp.Rational(11003, 240) * v_ip1**2 - sp.Rational(8623, 120) * v_ip1 * v_ip2 + sp.Rational(2321, 120) * v_ip1 * v_ip3 +
    sp.Rational(7043, 240) * v_ip2**2 - sp.Rational(647, 40) * v_ip2 * v_ip3 + sp.Rational(547, 240) * v_ip3**2
)

def complete_the_square(poly, var):
    poly_exp = sp.expand(poly)
    A = poly_exp.coeff(var, 2)
    B = poly_exp.coeff(var, 1)
    C = poly_exp - A * var**2 - B * var
    
    if A == 0:
        return 0, poly
    
    delta = -B / (2 * A)
    square_part = A * (var + delta)**2
    remaining = C - B**2 / (4 * A)
    
    return sp.simplify(square_part), sp.simplify(remaining)

def successive_completion(poly, order):
    current_poly = poly
    sos_terms = []
    
    for var in order:
        A = sp.expand(current_poly).coeff(var, 2)
        if A == 0:
            break
        square_part, remaining = complete_the_square(current_poly, var)
        sos_terms.append(square_part)
        current_poly = remaining
    
    if current_poly != 0:
        sos_terms.append(current_poly)
    
    return sos_terms

# 尝试顺序直到 diff = 0
orders = [
    [v_ip1, v_ip2, v_i, v_ip3],
    [v_ip2, v_ip1, v_i, v_ip3],
    [v_ip2, v_ip1, v_ip3, v_i],
    [v_ip1, v_ip3, v_ip2, v_i]
]

best_sos = None
for idx, order in enumerate(orders, 1):
    print(f"--- Order {idx}: {order} ---")
    sos_terms = successive_completion(beta0_expanded, order)
    sos = sum(sos_terms)
    diff = sp.simplify(sp.expand(sos) - beta0_expanded)
    print("Diff:", diff)
    if diff == 0:
        best_sos = sos
        print("Success! LaTeX SOS:")
        print(sp.latex(best_sos))
        break

if best_sos is None:
    print("No exact SOS found with these orders. Try larger steps or SDP.")

# PSD 检查
A = sp.Matrix([
    [sp.Rational(2107,240), sp.Rational(-1567,80), sp.Rational(3521,240), sp.Rational(-309,80)],
    [sp.Rational(-1567,80), sp.Rational(11003,240), sp.Rational(-8623,240), sp.Rational(2321,240)],
    [sp.Rational(3521,240), sp.Rational(-8623,240), sp.Rational(7043,240), sp.Rational(-647,80)],
    [sp.Rational(-309,80), sp.Rational(2321,240), sp.Rational(-647,80), sp.Rational(547,240)]
])
eigs_num = [float(e) for e in A.applyfunc(sp.N).eigenvals().values()]
print("Numerical eigenvalues (all >=0 for PSD):", eigs_num)