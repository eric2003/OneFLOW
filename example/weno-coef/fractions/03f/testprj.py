import numpy as np
from fractions import Fraction

def inverse_matrix(matrix):
    # 将矩阵元素转换为浮点数以计算逆矩阵
    matrix_float = matrix.astype(float)
    inverse = np.linalg.inv(matrix_float)
    # 将逆矩阵元素转换为分数
    inverse_fraction = [[Fraction(inverse[i, j]).limit_denominator() for j in range(len(inverse))] for i in range(len(inverse))]
    return inverse_fraction

# 定义一个函数，将矩阵中的每个元素转换为 Fraction 类型并打印
def print_matrix_fractionOLD(matrix):
    # 使用列表推导式和 Fraction 转换每个元素
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in matrix])
    
    # 打印转换后的矩阵
    print("Matrix in Fraction Form:")
    for row in fraction_matrix:
        # 将每个 Fraction 对象转换为字符串，并以指定格式打印
        print("[" + ", ".join(f"{f.numerator}/{f.denominator}" if f.denominator > 1 else str(f.numerator) for f in row) + "]")
        
# 定义一个函数，将矩阵中的每个元素转换为 Fraction 类型并打印
def print_matrix_fraction(matrix):
    # 使用列表推导式和 Fraction 转换每个元素
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in matrix])
    
    # 打印转换后的矩阵
    print("Matrix in Fraction Form:")
    for row in fraction_matrix:
        # 将每个 Fraction 对象转换为字符串，并以指定格式打印
        # 使用格式化字符串确保每个元素占用相同的宽度
        formatted_row = ", ".join(f"{f.numerator:2}/{f.denominator:2}" if f.denominator > 1 else f"{f.numerator:2}" for f in row)
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

