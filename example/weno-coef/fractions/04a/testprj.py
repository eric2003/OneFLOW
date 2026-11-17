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

# 示例用法
matrix = [
    [0.5, 3, 0.714285714],
    [10, 1.333333333, 6],
    [0.875, 11.111111111, 2]
]
print_matrix_fraction(matrix)