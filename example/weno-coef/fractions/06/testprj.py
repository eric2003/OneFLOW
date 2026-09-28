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
    #print("Matrix in Fraction Form:")
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
        
def compute_diff_coef(x,k):
    x0 = x - 1/2
    x1 = x + 1/2
    y = np.zeros(k)
    for j in range(k):
        j1 = j + 1
        y[j] = ( np.power(x1, j1) - np.power(x0, j1) ) / float(j1)
    return y
    
def compute_coef(x,k):
    y = np.zeros(k)
    for j in range(k):    
        y[j] = np.power(x, j)
    return y    
    
def create_matrix(r,s,k):
    arrays_list = []
    #k = r + s + 1
    for ii in range(-r,s+1):
        print(f'{ii=},{r=},{s=},{r+s=}')
        x = compute_diff_coef(ii,k)
        print(f"{x=}")
        # 动态地生成一些一维数组并添加到列表中
        arrays_list.append(x)
        
    # 使用 vstack 函数将列表中的数组堆叠成一个矩阵
    matrix = np.vstack(arrays_list)
    return matrix
    

def calc_eno_coef(r,s,k):    
    matrix = create_matrix(r,s,k)
    
    print_matrix_fraction(matrix)
    
    # 计算逆矩阵
    inverse = inverse_matrix(matrix)
    
    #print("\nInverse Matrix in Fraction Form:")
    #print_matrix_fraction(inverse)
    
    # 计算两个矩阵的乘积
    product = np.dot(matrix, inverse)
    
    #print("\nProduct of Matrix and Inverse Matrix:")
    #print_matrix_fraction(product)
    
    m1 = compute_coef(0.5,k).reshape(1, -1)
    #print_matrix_fraction(m1)
    
    mc = np.dot(m1, inverse)
    
    return mc
    
#k=3
k=1
k1 = k-1
print(f'{k=}')
#print(f'{k1=}')
arrays_list = []

for r in range(-1,k):
    s = k1 - r
    print(f'{r=},{s=},{k=}')
    mc = calc_eno_coef(r,s,k)
    arrays_list.append(mc)

# 使用 vstack 函数将列表中的所有数组垂直堆叠成一个矩阵
matrix = np.vstack(arrays_list)

print_matrix_fraction(matrix)
