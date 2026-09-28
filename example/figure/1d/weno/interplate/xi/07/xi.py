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
            #if f.denominator == 1:
            #    s = f"{f.numerator}"
            #else:
            #    s = f"{f.numerator}/{f.denominator}"
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

def integral_xi(x, j):
    return (x ** (j + 1)) / (j + 1)
    
def compute_coef(x,k):
    y = []
    for j in range(k):
        var = x ** j
        y.append(var)
    return y

def id_tostring(rj):
    mystr = str(rj)
    if rj == 0:
        mystr = '  '
    if rj > 0:
        mystr = '+' + str(rj)
    return mystr   
    
def coef_tostring(coef,i):
    mystr = str(coef)
    if coef >= 0:
        if i == 0:
            mystr = ' ' + mystr
        else:
            mystr = '+' + mystr
    return mystr
    
def coef_toabsstring(coef):
    abs_str = str(abs(coef))
    s = '+'
    if coef < 0:
        s = '-'
    return abs_str, s
    
def print_coef_formula(alist,xfrac,ishift=0):
    rows, cols = alist.shape
    #print(f'rows,cols={rows},{cols}')

    widths = np.empty(rows, dtype=int)

    for j in range(rows):
        w = 0
        for i in range(cols):
            absv,ss = coef_toabsstring(alist[i][j])
            ww = len( absv ) + 1
            w = max(w, ww)
        widths[j] = w
        
    #print(f'widths={widths}')
    
    for i in range(rows):
        mystr = ''
        r = i
        for j in range(cols):
            absv,ss = coef_toabsstring(alist[i][j])
            if j == 0 and ss == '+':
                ss = ' '
            tt = f"{absv:>{widths[j]-1}}"
            rj = ishift - r + j
            var_rj = id_tostring(rj)
            mystr += ss + tt + f"*v[i{var_rj}]"
            sxf = ''
            xfrac_new = xfrac + ishift
            if xfrac_new >= 0:
                sxf='+'
            slr='-'
            if xfrac < 0:
                slr='+'
        #print(f'vi{sxf}{xfrac_new}({slr}),r={r}={mystr}')
        print(f'vi{sxf}{xfrac_new}({slr}),{r}={mystr}')
    print()
    
def cal_polynomial_matrix(r, kk):
    arrays_list = []
    for m in range(kk):
        j = -r + m
        xia = Fraction(j) - Fraction(1,2)
        xib = Fraction(j) + Fraction(1,2)
        a_list = []
        for i in range(kk):
            val = integral_xi(xib, i) - integral_xi(xia, i)
            a_list.append(val)
        arrays_list.append(a_list)
    matrix = np.vstack(arrays_list)
    return matrix
    
def cal_polynomial_coefficients(r, kk, xfrac):
    matrix = cal_polynomial_matrix(r, kk)
    inverse = inverse_matrix(matrix)
    xv = compute_coef(xfrac, kk)
    yv = np.dot(xv, inverse)
    return yv
        
def calc_coef_formula(kk, xfrac):
    #xm = compute_coef(xfrac, kk)
    rows_list = []
    for r in range(kk):
        #-r+l
        ym = cal_polynomial_coefficients(r, kk, xfrac)
        rows_list.append(ym)
        
    return np.vstack(rows_list)
    
# i-2 i-1, i i+1 i+2
# vi+1/2(-),r=sum crl*vi-r+l,l=1,kk-1;
# kk=5
# r=0: -r+l=0,1,2,3,4:i,i+1,i+2,i+3,i+4
# r=1: -r+l=-1,0,1,2,3:i-1,i,i+1,i+2,i+3
# r=2: -r+l=-2,-1,0,1,2:i-2,i-1,i,i+1,i+2
        
    
kk=3
xfrac = Fraction(1,2)
xp = compute_coef(Fraction(1,2), kk)
xn = compute_coef(-Fraction(1,2), kk)

mymat1 = calc_coef_formula(kk, xfrac)
mymat2 = calc_coef_formula(kk, -xfrac)

print_matrix_fraction(mymat1)
print_coef_formula(mymat1,xfrac)

print_matrix_fraction(mymat2)
#print_coef_formula(mymat2,-xfrac,1)
print_coef_formula(mymat2,-xfrac)

