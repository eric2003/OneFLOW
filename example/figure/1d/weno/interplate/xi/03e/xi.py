import numpy as np
from fractions import Fraction

def inverse_matrix(matrix):
    # 将矩阵元素转换为浮点数以计算逆矩阵
    matrix_float = matrix.astype(float)
    inverse = np.linalg.inv(matrix_float)
    # 将逆矩阵元素转换为分数
    inverse_fraction = [[Fraction(inverse[i, j]).limit_denominator() for j in range(len(inverse))] for i in range(len(inverse))]
    return inverse_fraction

def print_matrix_fraction(matrix, is_column_vector=False):
    """
    支持一维向量和二维矩阵的分数字符串打印
    :param matrix: 一维列表（向量）或二维列表/数组（矩阵）
    :param is_column_vector: 一维向量是否按列向量格式打印（默认False：行向量）
    """
    # 步骤1：统一转换为二维矩阵格式
    if isinstance(matrix, (list, np.ndarray)):
        # 若为一维，转为二维（行向量：1×N 或 列向量：N×1）
        if np.ndim(matrix) == 1:
            if is_column_vector:
                # 列向量：N行1列
                two_d_matrix = [[x] for x in matrix]
            else:
                # 行向量：1行N列
                two_d_matrix = [matrix]
        else:
            # 若为二维，直接使用
            two_d_matrix = matrix
    else:
        raise TypeError("输入必须是列表或numpy数组")
    
    # 步骤2：转换为Fraction数组
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in two_d_matrix])
    rows = len(fraction_matrix)
    cols = len(fraction_matrix[0])
    
    # 步骤3：转换为字符串矩阵并计算每列最大宽度
    str_matrix = []
    col_widths = [0] * cols  # 每列的最大宽度
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
    
    # 步骤4：打印（保持原有对齐风格）
    for i in range(rows):
        row_elements = []
        for j in range(cols):
            element = str_matrix[i][j]
            formatted_element = f"{element:>{col_widths[j]}}"  # 右对齐
            if j < cols - 1:
                formatted_element += ", "
            else:
                formatted_element += " "
            row_elements.append(formatted_element)
        formatted_row = "".join(row_elements)
        print(f"[ {formatted_row}]")

def print_matrix_latex(matrix, matrix_name="A", is_column_vector=False):
    """
    支持一维向量和二维矩阵的LaTeX格式打印
    :param matrix: 一维列表（向量）或二维列表/数组（矩阵）
    :param matrix_name: 矩阵名称（注释用）
    :param is_column_vector: 一维向量是否按列向量格式打印（默认False：行向量）
    """
    # 步骤1：统一转换为二维矩阵格式
    if isinstance(matrix, (list, np.ndarray)):
        if np.ndim(matrix) == 1:
            if is_column_vector:
                two_d_matrix = [[x] for x in matrix]  # 列向量：N×1
            else:
                two_d_matrix = [matrix]  # 行向量：1×N
        else:
            two_d_matrix = matrix
    else:
        raise TypeError("输入必须是列表或numpy数组")
    
    # 步骤2：转换为Fraction数组
    fraction_matrix = np.array([[Fraction(x).limit_denominator() for x in row] for row in two_d_matrix])
    rows = len(fraction_matrix)
    cols = len(fraction_matrix[0])
    
    # 步骤3：构建LaTeX字符串
    latex_lines = []
    latex_lines.append(f"% LaTeX 格式：{matrix_name}")
    # 向量用bmatrix，矩阵也用bmatrix（统一风格，可改为pmatrix等）
    latex_lines.append("\\begin{bmatrix}")
    
    for i in range(rows):
        row_elements = []
        for j in range(cols):
            f = fraction_matrix[i][j]
            if f.denominator == 1:
                row_elements.append(f"{f.numerator}")
            else:
                row_elements.append(f"\\frac{{{f.numerator}}}{{{f.denominator}}}")
        # 行内元素用&分隔，行尾加\\
        latex_lines.append(" & ".join(row_elements) + " \\\\")
    
    latex_lines.append("\\end{bmatrix}")
    latex_str = "\n".join(latex_lines)
    
    # 步骤4：打印LaTeX格式
    print(f"\n{'='*50}")
    print(f"LaTeX 格式 - {matrix_name}:")
    print(latex_str)
    print(f"{'='*50}\n")

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)
    
def compute_coef(x,k):
    y = []
    for j in range(k):
        var = x ** j
        y.append(var)
    return y

vv = [-1,0,1]
k = len(vv)

# 测试一维向量xx的打印（支持行向量/列向量两种格式）
xx = compute_coef(Fraction(1,2), k)
print(f"xx（一维行向量）:")
print_matrix_fraction(xx)  # 默认行向量格式
print_matrix_latex(xx, "Vector xx (Row Vector)")  # LaTeX行向量

# 可选：按列向量格式打印xx
print(f"\nxx（一维列向量）:")
print_matrix_fraction(xx, is_column_vector=True)
print_matrix_latex(xx, "Vector xx (Column Vector)", is_column_vector=True)

# 构建二维矩阵并打印
arrays_list = []
for j in vv:
    xia = Fraction(j) - Fraction(1,2)
    xib = Fraction(j) + Fraction(1,2)
    a_list = []
    for i in range(k):
        val = integral_xi(xib, i) - integral_xi(xia, i)
        a_list.append(val)
    arrays_list.append(a_list)

matrix = np.vstack(arrays_list)
print("\nOriginal Matrix in Fraction Form:")
print_matrix_fraction(matrix)
print_matrix_latex(matrix, "Original Matrix")

# 计算并打印逆矩阵
inverse = inverse_matrix(matrix)
print(f"inverse（二维矩阵）:")
print_matrix_fraction(inverse)
print_matrix_latex(inverse, "Inverse Matrix")

# 计算并打印乘积矩阵
product = np.dot(matrix, inverse)
print("Product of Matrix and Inverse Matrix:")
print_matrix_fraction(product)
print_matrix_latex(product, "Matrix Product (Identity Matrix)")

yy = np.dot(xx, inverse)
print(f"yy:")
print_matrix_fraction(yy) 