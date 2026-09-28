from sympy import Rational
from itertools import product
from math import gcd as math_gcd

# beta0 的对称矩阵（二次型 A，off-diag 是跨项系数的一半）
A_beta = [
    [Rational(10,3), Rational(-31,6), Rational(11,6)],
    [Rational(-31,6), Rational(25,3), Rational(-19,6)],
    [Rational(11,6), Rational(-19,6), Rational(4,3)]
]

def gcd_list(lst):
    """计算列表整数的最大公约数（处理 0）。"""
    g = 0
    for x in lst:
        g = math_gcd(g, abs(x))
    return g

def normalize_vector(vec):
    """标准化向量：gcd=1，第一非零系数正（忽略零系数）。"""
    non_zero = [x for x in vec if x != 0]
    if not non_zero:
        return tuple(vec)
    g = gcd_list(non_zero)
    vec_norm = tuple(int(x) // g if x != 0 else 0 for x in vec)
    # 翻转符号使第一非零正
    for i in range(len(vec_norm)):
        if vec_norm[i] != 0:
            if vec_norm[i] < 0:
                vec_norm = tuple(-int(x) if x != 0 else 0 for x in vec_norm)
            break
    return vec_norm

def check_sos_2_terms(l1, l2, A):
    """检查是否为有效 2 平方 SOS：解 lambda 并验证所有系数。"""
    a, b, c = l1
    d, e, f = l2
    p = a**2
    q = d**2
    r = b**2
    s = e**2
    det = p * s - q * r
    if det == 0:
        return None
    # Cramer 法则解 lambda1, lambda2 (基于 v0² 和 v1² 系数)
    lambda1 = (A[0][0] * s - A[1][1] * q) / det
    lambda2 = (A[1][1] * p - A[0][0] * r) / det
    # 预测剩余系数
    pred_22 = lambda1 * c**2 + lambda2 * f**2
    pred_01 = lambda1 * a * b + lambda2 * d * e
    pred_02 = lambda1 * a * c + lambda2 * d * f
    pred_12 = lambda1 * b * c + lambda2 * e * f
    # 验证（对称矩阵，只查上三角）
    if (pred_22 == A[2][2] and pred_01 == A[0][1] and pred_02 == A[0][2] and pred_12 == A[1][2] and
        lambda1 > 0 and lambda2 > 0):
        return (l1, lambda1, l2, lambda2)
    return None

def find_all_sos_2(A, max_c=4):
    """枚举所有 2 平方 SOS（包括稀疏支持集，如 3+2, 2+2）。"""
    print(f"Starting search with max_c={max_c}... Total iterations: {(2*max_c + 1)**6}")
    sos_list = []
    processed = set()  # 去重已标准化对
    total = (2 * max_c + 1) ** 6
    count = 0
    for coeffs in product(range(-max_c, max_c + 1), repeat=6):
        count += 1
        if count % 50000 == 0:  # 每 5 万迭代打印进度
            print(f"Progress: {count}/{total} ({count / total * 100:.1f}%) | Processed pairs: {len(processed)}")
        l1_raw = coeffs[0:3]
        l2_raw = coeffs[3:6]
        if all(x == 0 for x in l1_raw) or all(x == 0 for x in l2_raw):  # 跳过零向量
            continue
        l1 = normalize_vector(l1_raw)
        l2 = normalize_vector(l2_raw)
        pair = tuple(sorted([l1, l2]))  # 排序避免 (l1,l2) 和 (l2,l1) 重复
        if pair in processed:
            continue
        processed.add(pair)
        res = check_sos_2_terms(l1, l2, A)
        if res:
            sos_list.append(res)
            print(f"Found one: {res}")  # 立即打印发现
    print(f"Search complete. Found {len(sos_list)} representations.")
    return sos_list

def classify_supports(sos_list):
    """分类找到的形式按支持集大小 (e.g., 3+3, 3+2)。"""
    from collections import defaultdict
    groups = defaultdict(list)
    for l1, lam1, l2, lam2 in sos_list:
        supp1 = sum(1 for x in l1 if x != 0)
        supp2 = sum(1 for x in l2 if x != 0)
        key = f"{max(supp1, supp2)}+{min(supp1, supp2)}"  # 排序如 3+2
        groups[key].append((l1, lam1, l2, lam2))
    return groups

# 运行搜索与分类
if __name__ == "__main__":
    sos_2 = find_all_sos_2(A_beta, max_c=4)
    print("\n=== 所有找到的 2 平方 SOS 表示（按支持集分类） ===")
    groups = classify_supports(sos_2)
    for combo, forms in sorted(groups.items()):
        print(f"\n{combo} 组合 ({len(forms)} 个):")
        for i, (l1, lam1, l2, lam2) in enumerate(forms, 1):
            print(f"  {i}. \\beta_0 = {lam1} ({l1[0]} v_0 + {l1[1]} v_1 + {l1[2]} v_2)^2 + {lam2} ({l2[0]} v_0 + {l2[1]} v_1 + {l2[2]} v_2)^2")
    
    print("\n注: 无找到 3+2 或 2+2（max_c=4 下）；增大 max_c=6 可能发现更多稀疏形式，但计算 ~10x 更长。")