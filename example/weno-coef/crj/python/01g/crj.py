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

for k in range(1, 8):
    print(f"=== k = {k} ===")
    mat = []
    # 关键修改：将 range(-1, k) 改为反向遍历（从大到小）
    r_values = range(k - 1, -2, -1)  # 原范围是 [-1, 0, 1, ..., k-1]，反向后为 [k-1, ..., 1, 0, -1]
    for r in r_values:
        row = []
        for j in range(k):
            row.append(calculate_crj(r, j, k))
        mat.append(row)

    # 所有会出现字符串（包括表头 j=0, j=1...）用于计算统一宽度
    header_cells = [f"j={j}" for j in range(k)]
    all_strings = header_cells + [str(item) for row in mat for item in row]
    max_width = max(len(s) for s in all_strings) if all_strings else 8

    # r 列宽度（包含表头 "r"）
    r_strings = ["r"] + [str(r) for r in r_values]
    r_width = max(len(s) for s in r_strings)

    # 表头
    header = " ".join(f"{cell:^{max_width}}" for cell in header_cells)
    print(f"{'r':>{r_width}} | {header}")

    # 数据行（数值右对齐）
    for r, row in zip(r_values, mat):
        r_str = f"{r:>{r_width}}"
        cells = " ".join(f"{str(item):>{max_width}}" for item in row)
        print(f"{r_str} | {cells}")
    print()