from fractions import Fraction

def calculate_crj(r, j, k):
    result = Fraction(0, 1)
    for m in range(j + 1, k + 1):
        numerator = 0
        for l in range(0, k + 1):
            if l == m:
                continue
            product = 1
            for q in range(0, k + 1):
                if q == m or q == l:
                    continue
                product *= (r - q + 1)
            numerator += product
        denominator = 1
        for l in range(0, k + 1):
            if l == m:
                continue
            denominator *= (m - l)
        result += Fraction(numerator, denominator)
    return result

for k in range(1,8):
    print(f"=== k = {k} ===")
    mat = []
    r_values = range(-1, k)
    for r in r_values:
        row = []
        for j in range(k):
            row.append(calculate_crj(r, j, k))
        mat.append(row)

    # 计算宽度（自动处理分数/负数长度）
    max_width = max(len(str(item)) for row in mat for item in row)
    r_width = max(len(str(r)) for r in r_values)
    r_width = max(r_width, len("r"))  # 表头“r”也占位

    # 打印表头
    header = " ".join(f"{j:^{max_width}}" for j in range(k))
    print(f"{'r':>{r_width}} | {header}")

    # 打印数据行
    for r, row in zip(r_values, mat):
        r_str = f"{r:>{r_width}}"
        cells = " ".join(f"{str(item):^{max_width}}" for item in row)
        print(f"{r_str} | {cells}")
    print()  # 每个 k 后空一行分隔