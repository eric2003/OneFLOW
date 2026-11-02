import numpy as np
from fractions import Fraction

def print_matrix_fraction(matrix):
    # 将矩阵转换为Fraction数组
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in matrix])
    
    # 转换为字符串矩阵并计算最大宽度
    str_matrix = []
    max_width = 0
    for row in fraction_matrix:
        str_row = []
        for f in row:
            if f.denominator == 1:
                s = f"{f.numerator}"
            else:
                s = f"{f.numerator}/{f.denominator}"
            str_row.append(s)
            current_length = len(s)
            if current_length > max_width:
                max_width = current_length
        str_matrix.append(str_row)
    
    # 打印矩阵，每个元素等宽对齐
    print("Matrix in Fraction Form:")
    for row in str_matrix:
        formatted_elements = [f"{s:>{max_width}}" for s in row]
        formatted_row = "  ".join(formatted_elements)
        print(f"[ {formatted_row} ]")

# 示例用法
matrix = [
    [0.5, 3, 0.714285714],
    [10, 1.333333333, 6],
    [0.875, 11.111111111, 2]
]
print_matrix_fraction(matrix)