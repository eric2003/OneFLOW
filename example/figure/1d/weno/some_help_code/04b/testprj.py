from sympy import Rational
from itertools import product
from math import gcd as math_gcd  # 更快整数 gcd

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
    """标准化向量：gcd=1，第一非零系数正。"""
    g = gcd_list(vec)
    if g == 0:
        return tuple(vec)
    vec = tuple(int(x) // g for x in vec)
    # 翻转符号使第一非零正
    for i in range(len(vec)):
        if vec[i] != 0:
            if vec[i] < 0:
                vec = tuple(-int(x) for x in vec)
            break
    return vec

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
    # 预测所有系数
    pred_00 = lambda1 * p + lambda2 * q
    pred_11 = lambda1 * r + lambda2 * s
    pred_22 = lambda1 * c**2 + lambda2 * f**2
    pred_01 = lambda1 * a * b + lambda2 * d * e
    pred_02 = lambda1 * a * c + lambda2 * d * f
    pred_12 = lambda1 * b * c + lambda2 * e * f
    # 验证（对称矩阵，只查上三角）
    if (pred_00 == A[0][0] and pred_11 == A[1][1] and pred_22 == A[2][2] and
        pred_01 == A[0][1] and pred_02 == A[0][2] and pred_12 == A[1][2] and
        lambda1 > 0 and lambda2 > 0):
        return (l1, lambda1, l2, lambda2)
    return None

def find_all_sos_2(A, max_c=4):  # 默认 max_c=4 以快速找到目标形式
    """枚举所有 2 平方 SOS，小系数版本。"""
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

# 运行 2 平方搜索（从小 max_c 开始）
if __name__ == "__main__":
    sos_2 = find_all_sos_2(A_beta, max_c=4)  # 改成 6 如果想找更多
    print("\n=== 找到的 2 平方 SOS 表示（3 项线性形式 + 2 平方） ===")
    for i, (l1, lam1, l2, lam2) in enumerate(sos_2, 1):
        print(f"{i}. \\beta_0 = {lam1} ({l1[0]} v_0 + {l1[1]} v_1 + {l1[2]} v_2)^2 + {lam2} ({l2[0]} v_0 + {l2[1]} v_1 + {l2[2]} v_2)^2")
    
    print("\n注: 3 平方 SOS 可能无小系数严格正 lambda 表示（因秩 2），可通过添加第三个 0 lambda 退化到 2 平方。")
    print("若需 3 平方扩展，增大 max_c 但计算量 O(max_c^9)，建议 max_c=2。")