from sympy import Rational
from itertools import product
from math import gcd as math_gcd
from collections import defaultdict  # 添加以防

# 修正矩阵（off-diagonal 非 0）
A_beta0 = [
    [Rational(10,3), Rational(-31,6), Rational(11,6)],
    [Rational(-31,6), Rational(25,3), Rational(-19,6)],
    [Rational(11,6), Rational(-19,6), Rational(4,3)]
]

A_beta1 = [
    [Rational(4,3), Rational(-13,6), Rational(5,6)],
    [Rational(-13,6), Rational(13,3), Rational(-13,6)],
    [Rational(5,6), Rational(-13,6), Rational(4,3)]
]

A_beta2 = [
    [Rational(4,3), Rational(-19,6), Rational(11,6)],
    [Rational(-19,6), Rational(25,3), Rational(-31,6)],
    [Rational(11,6), Rational(-31,6), Rational(10,3)]
]

ALL_BETAS = {'beta0': A_beta0, 'beta1': A_beta1, 'beta2': A_beta2}

def gcd_list(lst):
    g = 0
    for x in lst:
        g = math_gcd(g, abs(x))
    return g

def normalize_vector(vec):
    non_zero = [x for x in vec if x != 0]
    if not non_zero:
        return tuple(vec)
    g = gcd_list(non_zero)
    vec_norm = tuple(int(x) // g if x != 0 else 0 for x in vec)
    for i in range(len(vec_norm)):
        if vec_norm[i] != 0:
            if vec_norm[i] < 0:
                vec_norm = tuple(-int(x) if x != 0 else 0 for x in vec_norm)
            break
    return vec_norm

def check_sos_2_terms(l1, l2, A):
    a, b, c = l1
    d, e, f = l2
    p = a**2
    q = d**2
    r = b**2
    s = e**2
    det = p * s - q * r
    if det == 0:
        return None
    lambda1 = (A[0][0] * s - A[1][1] * q) / det
    lambda2 = (A[1][1] * p - A[0][0] * r) / det
    pred_22 = lambda1 * c**2 + lambda2 * f**2
    pred_01 = lambda1 * a * b + lambda2 * d * e
    pred_02 = lambda1 * a * c + lambda2 * d * f
    pred_12 = lambda1 * b * c + lambda2 * e * f
    if (pred_22 == A[2][2] and pred_01 == A[0][1] and pred_02 == A[0][2] and pred_12 == A[1][2] and
        lambda1 > 0 and lambda2 > 0):
        return (l1, lambda1, l2, lambda2)
    return None

def find_all_sos_2(A, max_c=4, name='beta'):
    print(f"\n--- Searching for {name} with max_c={max_c}... ---")
    print(f"Total iterations: {(2*max_c + 1)**6}")
    sos_list = []
    processed = set()
    total = (2 * max_c + 1) ** 6
    count = 0
    for coeffs in product(range(-max_c, max_c + 1), repeat=6):
        count += 1
        if count % 50000 == 0:
            print(f"Progress: {count}/{total} ({count / total * 100:.1f}%) | Processed pairs: {len(processed)}")
        l1_raw = coeffs[0:3]
        l2_raw = coeffs[3:6]
        if all(x == 0 for x in l1_raw) or all(x == 0 for x in l2_raw):
            continue
        l1 = normalize_vector(l1_raw)
        l2 = normalize_vector(l2_raw)
        pair = tuple(sorted([l1, l2]))
        if pair in processed:
            continue
        processed.add(pair)
        res = check_sos_2_terms(l1, l2, A)
        if res:
            sos_list.append(res)
            print(f"Found one for {name}: {res}")
    print(f"Search for {name} complete. Found {len(sos_list)} representations.")
    return sos_list

def classify_supports(sos_list):
    groups = defaultdict(list)
    for l1, lam1, l2, lam2 in sos_list:
        supp1 = sum(1 for x in l1 if x != 0)
        supp2 = sum(1 for x in l2 if x != 0)
        key = f"{max(supp1, supp2)}+{min(supp1, supp2)}"
        groups[key].append((l1, lam1, l2, lam2))
    return groups

if __name__ == "__main__":
    max_c = 4  # 可增大到 6
    all_results = {}
    for name, A in ALL_BETAS.items():
        sos = find_all_sos_2(A, max_c, name)
        all_results[name] = classify_supports(sos)
    
    print("\n=== 所有 β 的 SOS 表示（按支持集分类） ===")
    for name, groups in all_results.items():
        print(f"\n{name.upper()}:")
        for combo, forms in sorted(groups.items()):
            print(f"  {combo} 组合 ({len(forms)} 个):")
            for i, (l1, lam1, l2, lam2) in enumerate(forms, 1):
                print(f"    {i}. \\beta_{name} = {lam1} ({l1[0]} v_0 + {l1[1]} v_1 + {l1[2]} v_2)^2 + {lam2} ({l2[0]} v_0 + {l2[1]} v_1 + {l2[2]} v_2)^2")
    
    print(f"\n注: max_c={max_c} 下结果；若无稀疏形式，增大 max_c。总运行时间 ~2min。")