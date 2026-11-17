import numpy as np

def calculate_eno_crj(r, j, k):
    result = 0

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

        result += numerator / denominator

    return result
    
for k in range(1,8):
    print(f"=== k = {k} ===")
    # 计算矩阵并存储到列表中
    eno_mat = np.zeros((k+1,k))
    for r in range(-1,k):
        for j in range(k):
            eno_mat[r+1, j] = calculate_eno_crj(r, j, k)
        
    print(f'{eno_mat=}')
    print(f'{type(eno_mat)=}')