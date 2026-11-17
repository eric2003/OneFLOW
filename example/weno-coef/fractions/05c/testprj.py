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
    print("Matrix in Fraction Form:")
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
        
def compute_diff_coef(x):
    x0 = x - 1/2
    x1 = x + 1/2
    y = np.zeros(3)
    y[0] = ( np.power(x1, 1) - np.power(x0, 1) ) / float(1)
    y[1] = ( np.power(x1, 2) - np.power(x0, 2) ) / float(2)
    y[2] = ( np.power(x1, 3) - np.power(x0, 3) ) / float(3)
    return y
    
def compute_coef(x):
    y = np.zeros(3)
    y[0] = np.power(x, 0)
    y[1] = np.power(x, 1)
    y[2] = np.power(x, 2)
    return y    
    
def create_matrix():
    x0 = compute_diff_coef(-2)
    x1 = compute_diff_coef(-1)
    x2 = compute_diff_coef(0)
    print(f"{x0=}")
    print(f"{x1=}")
    print(f"{x2=}")

    arrays_list = []

    # 动态地生成一些一维数组并添加到列表中
    arrays_list.append(x0)
    arrays_list.append(x1)
    arrays_list.append(x2)

    # 使用 vstack 函数将列表中的数组堆叠成一个矩阵
    matrix = np.vstack(arrays_list)
    return matrix
    
def create_matrix_new(k):
    k1 = k-1
    print(f'{k1=}')
    
    for s in range(k):
        r = k1 - s
        print(f'{r=},{s=},{k=}')
    r = 2
    s = 0
    arrays_list = []
    for ii in range(-r,s+1):
        print(f'{ii=},{r=},{s=}')
        x = compute_diff_coef(ii)
        print(f"{x=}")
        # 动态地生成一些一维数组并添加到列表中
        arrays_list.append(x)
    #x0 = compute_diff_coef(-2)
    #x1 = compute_diff_coef(-1)
    #x2 = compute_diff_coef(0)
    #print(f"{x0=}")
    #print(f"{x1=}")
    #print(f"{x2=}")


    # 动态地生成一些一维数组并添加到列表中
    #arrays_list.append(x0)
    #arrays_list.append(x1)
    #arrays_list.append(x2)

    # 使用 vstack 函数将列表中的数组堆叠成一个矩阵
    matrix = np.vstack(arrays_list)
    return matrix    
    
#matrix = create_matrix() 
k=3
matrix = create_matrix_new(k)

print_matrix_fraction(matrix)

# 计算逆矩阵
inverse = inverse_matrix(matrix)

print("\nInverse Matrix in Fraction Form:")
print_matrix_fraction(inverse)

# 计算两个矩阵的乘积
product = np.dot(matrix, inverse)

print("\nProduct of Matrix and Inverse Matrix:")
print_matrix_fraction(product)

m1 = compute_coef(0.5).reshape(1, -1)
print_matrix_fraction(m1)

mc = np.dot(m1, inverse)

print("\n Matrix C:")
print_matrix_fraction(mc)


