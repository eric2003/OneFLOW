import numpy as np
from fractions import Fraction

def inverse_matrix(matrix):
    # 将矩阵元素转换为浮点数以计算逆矩阵
    matrix_float = matrix.astype(float)
    inverse = np.linalg.inv(matrix_float)
    # 将逆矩阵元素转换为分数
    inverse_fraction = [[Fraction(inverse[i, j]).limit_denominator() for j in range(len(inverse))] for i in range(len(inverse))]
    return inverse_fraction

# 定义一个函数，将矩阵中的每个元素转换为 Fraction 类型并打印，行列等宽对齐
def print_matrix_fraction(matrix):
    # 找到每列最大宽度
    col_widths = [0] * matrix.shape[1]
    for row in matrix:
        for i, elem in enumerate(row):
            # 将元素转换为字符串形式，包括分子和分母
            elem_str = f"{elem.numerator}/{elem.denominator}" if hasattr(elem, 'denominator') and elem.denominator > 1 else str(elem)
            col_widths[i] = max(col_widths[i], len(elem_str))
    
    # 打印转换后的矩阵
    print("Matrix in Fraction Form:")
    for row in matrix:
        formatted_row = ", ".join(f"{{:^{col_widths[i]}}}".format(f"{elem.numerator}/{elem.denominator}" if hasattr(elem, 'denominator') and elem.denominator > 1 else elem) for i, elem in enumerate(row))
        print(f"[{formatted_row}]")


# 定义原始矩阵 matrix
matrix = np.array([[1, -2, 49/12], [1, -1, 13/12], [1, 0, 1/12]])

# 计算逆矩阵
inverse = inverse_matrix(matrix)

# 调用函数打印矩阵和逆矩阵
print_matrix_fraction(matrix)
print("\nInverse Matrix in Fraction Form:")
print_matrix_fraction(inverse)

# 计算两个矩阵的乘积
product = np.dot(matrix, inverse)

print("\nProduct of Matrix and Inverse Matrix:")
print_matrix_fraction(product)

