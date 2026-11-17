import numpy as np
from fractions import Fraction

def inverse_matrix(matrix):
    inverse = np.linalg.inv(matrix)
    inverse_fraction = [[Fraction(inverse[i, j]).limit_denominator() for j in range(len(inverse))] for i in range(len(inverse))]
    return inverse_fraction

# 定义一个函数，将矩阵中的每个元素转换为 Fraction 类型并打印
def print_matrix_fraction(matrix):
    # 使用列表推导式和 Fraction 转换每个元素
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in matrix])
    
    # 打印转换后的矩阵
    print("Matrix in Fraction Form:")
    for row in fraction_matrix:
        # 将每个 Fraction 对象转换为字符串，并以指定格式打印
        print("[" + ", ".join(f"{f.numerator}/{f.denominator}" if f.denominator > 1 else str(f.numerator) for f in row) + "]")
        
matrix = np.array([[1, -2, 49/12], [1, -1, 13/12], [1, 0, 1/12]])
inverse = inverse_matrix(matrix)
        
# 调用函数打印矩阵
print_matrix_fraction(matrix)
print_matrix_fraction(inverse)