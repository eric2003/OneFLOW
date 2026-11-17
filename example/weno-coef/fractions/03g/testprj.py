import numpy as np
from fractions import Fraction

def inverse_matrix(matrix):
    # 将矩阵元素转换为浮点数以计算逆矩阵
    matrix_float = matrix.astype(float)
    inverse = np.linalg.inv(matrix_float)
    # 将逆矩阵元素转换为分数
    inverse_fraction = [[Fraction(inverse[i, j]).limit_denominator() for j in range(len(inverse))] for i in range(len(inverse))]
    return inverse_fraction

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
    
    # 打印矩阵，每个元素等宽右对齐，添加逗号
    print("Matrix in Fraction Form:")
    rows = len(str_matrix)
    cols = len(str_matrix[0])
    
    for i in range(rows):
        row_elements = []
        for j in range(cols):
            element = str_matrix[i][j]
            # 右对齐
            formatted_element = f"{element:>{max_width}}"
            # 除最后一列外添加逗号和空格
            if j < cols - 1:
                formatted_element += ", "
            else:
                formatted_element += " "
            row_elements.append(formatted_element)
        # 拼接一行并打印
        formatted_row = "".join(row_elements)
        print(f"[ {formatted_row}]")

 

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

