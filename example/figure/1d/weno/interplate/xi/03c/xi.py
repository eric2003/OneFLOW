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
    
    # 转换为字符串矩阵并计算每列的最大宽度
    str_matrix = []
    rows = len(fraction_matrix)
    cols = len(fraction_matrix[0])
    col_widths = [0] * cols  # 每列的最大宽度
    
    # 将数字转换为字符串，并记录每列最大宽度
    for row in fraction_matrix:
        str_row = []
        for j, f in enumerate(row):
            if f.denominator == 1:
                s = f"{f.numerator}"
            else:
                s = f"{f.numerator}/{f.denominator}"
            str_row.append(s)
            current_length = len(s)
            if current_length > col_widths[j]:
                col_widths[j] = current_length
        str_matrix.append(str_row)
    
    # 打印矩阵，每列等宽右对齐，添加逗号
    for i in range(rows):
        row_elements = []
        for j in range(cols):
            element = str_matrix[i][j]
            # 右对齐，使用该列的最大宽度
            formatted_element = f"{element:>{col_widths[j]}}"
            # 除最后一列外添加逗号和空格
            if j < cols - 1:
                formatted_element += ", "
            else:
                formatted_element += " "
            row_elements.append(formatted_element)
        # 拼接一行并打印
        formatted_row = "".join(row_elements)
        print(f"[ {formatted_row}]")

def print_matrix_latex(matrix, matrix_name="A"):
    """
    输出矩阵的 LaTeX 格式
    :param matrix: 输入矩阵（列表嵌套或numpy数组）
    :param matrix_name: 矩阵名称（用于注释）
    """
    # 转换为Fraction数组
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in matrix])
    rows = len(fraction_matrix)
    cols = len(fraction_matrix[0])
    
    # 构建LaTeX字符串
    latex_lines = []
    latex_lines.append(f"% LaTeX 格式：{matrix_name}")
    latex_lines.append("\\begin{bmatrix}")
    
    for i in range(rows):
        row_elements = []
        for j in range(cols):
            f = fraction_matrix[i][j]
            if f.denominator == 1:
                # 整数直接输出分子
                row_elements.append(f"{f.numerator}")
            else:
                # 分数输出为 \frac{分子}{分母}
                row_elements.append(f"\\frac{{{f.numerator}}}{{{f.denominator}}}")
        # 每行元素用 & 分隔，行尾加 \\
        latex_lines.append(" & ".join(row_elements) + " \\\\")
    
    latex_lines.append("\\end{bmatrix}")
    latex_str = "\n".join(latex_lines)
    
    # 打印LaTeX格式（用 === 分隔，方便复制）
    print(f"\n{'='*50}")
    print(f"LaTeX 格式 - {matrix_name}:")
    print(latex_str)
    print(f"{'='*50}\n")

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)

vv = [-1,0,1]
isize = len(vv)
arrays_list = []
for j in vv:
    xia = Fraction(j) - Fraction(1,2)
    xib = Fraction(j) + Fraction(1,2)
    a_list = []
    for i in range(isize):
        val = integral_xi(xib, i) - integral_xi(xia, i)
        a_list.append(val)
    arrays_list.append(a_list)

# 使用 vstack 函数将列表中的数组堆叠成一个矩阵
matrix = np.vstack(arrays_list)

# 打印原始矩阵（分数字符串格式 + LaTeX格式）
print("Original Matrix in Fraction Form:")
print_matrix_fraction(matrix)
print_matrix_latex(matrix, "Original Matrix")

# 计算逆矩阵
inverse = inverse_matrix(matrix)

# 打印逆矩阵（分数字符串格式 + LaTeX格式）
print("Inverse Matrix in Fraction Form:")
print_matrix_fraction(inverse)
print_matrix_latex(inverse, "Inverse Matrix")

# 计算两个矩阵的乘积
product = np.dot(matrix, inverse)

# 打印乘积矩阵（分数字符串格式 + LaTeX格式）
print("Product of Matrix and Inverse Matrix:")
print_matrix_fraction(product)
print_matrix_latex(product, "Matrix Product (Identity Matrix)")