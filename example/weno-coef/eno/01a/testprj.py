import numpy as np
from fractions import Fraction

# 定义原始矩阵 original_matrix
original_matrix = np.array([[1, 0, 1/12], [1, 1, 13/12], [1, 2, 49/12]])

# 定义一个函数，将矩阵中的每个元素转换为 Fraction 类型
def print_matrix_fraction(matrix):
    # 使用列表推导式和 Fraction 转换每个元素
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in matrix])
    
    # 打印转换后的矩阵
    print("Matrix in Fraction Form:")
    for row in fraction_matrix:
        print(row)

# 调用函数打印矩阵
print_matrix_fraction(original_matrix)