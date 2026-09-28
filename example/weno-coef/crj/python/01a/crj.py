import math
from fractions import Fraction  

def divide_as_simplified_fraction(a, b):
    if b == 0:
        raise ValueError("分母不能为零！")
    
    # 计算最大公约数（GCD）
    gcd_value = math.gcd(a, b)
    
    # 化简分数
    simplified_num = a // gcd_value
    simplified_den = b // gcd_value
    
    # 返回最简分数形式
    return f"{simplified_num}/{simplified_den}"

def calculate_crj(r, j, k):
    result = Fraction(0, 1)

    # 外层求和：m 从 j+1 到 k
    for m in range(j + 1, k + 1):
        numerator = 0

        # 内层求和：l 从 0 到 k，且 l ≠ m
        for l in range(0, k + 1):
            if l == m:
                continue  # 跳过 l = m 的情况

            product = 1

            # 内层连乘：q 从 0 到 k，且 q ≠ m, l
            for q in range(0, k + 1):
                if q == m or q == l:
                    continue  # 跳过 q = m 或 q = l 的情况
                product *= (r - q + 1)

            numerator += product

        denominator = 1

        # 分母连乘：l 从 0 到 k，且 l ≠ m
        for l in range(0, k + 1):
            if l == m:
                continue  # 跳过 l = m 的情况
            denominator *= (m - l)

        #result += numerator / denominator
        result += Fraction(numerator, denominator)
    return result

# 计算并输出 k=1, 2, 3 时的所有 c_rj
for k in [1, 2, 3]:
    print(f"=== k = {k} ===")
    for r in range(0, k):
        for j in range(0, k):
            c_rj = calculate_crj(r, j, k)
            print(f"c_{r}{j} = {c_rj}")
    print()