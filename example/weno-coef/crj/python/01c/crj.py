import math
from fractions import Fraction  

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

        result += Fraction(numerator, denominator)
    return result
    

for k in range(1,8):
    print(f"=== k = {k} ===")
    # 计算矩阵并存储到列表中
    mat = []
    for r in range(k):
        row = []
        for j in range(k):
            c_rj = calculate_crj(r, j, k)
            row.append(c_rj)
        mat.append(row)
    # 计算每个分数的最大宽度
    max_width = max(len(str(item)) for row in mat for item in row)        
    for row in mat:
        for item in row:
            print(f"{str(item):^{max_width}}", end=' ')
        print()
    print()