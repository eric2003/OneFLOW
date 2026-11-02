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
    rows_list = []
    for r in range(kk):
        #-r+l
        ym = cal_polynomial_coefficients(r, kk, xfrac)
        rows_list.append(ym)
        
    return np.vstack(rows_list)
    
def cal_iplus_index(jstart,cols):
    var_rj_list=[]
    for j in range(cols):
        rj = jstart + j
        var_rj = id_tostring(rj)
        var_rj_list.append(var_rj)
    return var_rj_list
    
def print_coef_formula(alist,xfrac,ishift=0,rin=0):
    rows, cols = alist.shape
    #print(f'rows,cols={rows},{cols}')

    widths = np.empty(cols, dtype=int)

    for j in range(cols):
        w = 0
        for i in range(rows):
            absv,ss = coef_toabsstring(alist[i][j])
            ww = len( absv ) + 1
            w = max(w, ww)
        widths[j] = w
        
    #print(f'widths={widths}')
    
    sxf = ''
    xfrac_new = xfrac + ishift
    if xfrac_new >= 0:
        sxf='+'
    slr='-'
    if xfrac < 0:
        slr='+'
    vstr = f'vi{sxf}{xfrac_new}({slr})'
    for i in range(rows):
        r = i
        row = alist[i]
        if rows == 1:
            r = rin
        mystr = ''
        jstart = ishift - r
        var_rj_list = cal_iplus_index(jstart,cols)
        for j in range(cols):
            absv,ss = coef_toabsstring(row[j])
            if j == 0 and ss == '+':
                ss = ' '
            ss = f"{ss}{absv:>{widths[j]-1}}"
            mystr += ss + f"*v[i{var_rj_list[j]}]"
        print(f'{vstr},{r}={mystr}')
    print()
    
def printhhh(row_matrix, mymat):
    rows_ref, cols_ref = row_matrix.shape
    print(f'rows_ref,cols_ref={rows_ref},{cols_ref}')
    rows, cols = mymat.shape
    print(f'rows,cols={rows},{cols}')
    """
    [ 1/30, -13/60, 47/60, 9/20, -1/20 ]
    vi+1/2(-),2= 1/30*v[i-2]-13/60*v[i-1]+47/60*v[i  ]+9/20*v[i+1]-1/20*v[i+2]
    
    [  1/3,  5/6, -1/6 ]
    [ -1/6,  5/6,  1/3 ]
    [  1/3, -7/6, 11/6 ]
    vi+1/2(-),0= 1/3*v[i  ]+5/6*v[i+1]- 1/6*v[i+2]
    vi+1/2(-),1=-1/6*v[i-1]+5/6*v[i  ]+ 1/3*v[i+1]
    vi+1/2(-),2= 1/3*v[i-2]-7/6*v[i-1]+11/6*v[i  ]
    rows_ref,cols_ref=1,5
    rows,cols=3,3
    r=0:0 1 2
    r=1:-1 0 1
    r=2:-2 -1 0
    -2 [(2, Fraction(1, 3))]
    -1 [(1, Fraction(-1, 6)), (2, Fraction(-7, 6))]
    0 [(0, Fraction(1, 3)), (1, Fraction(5, 6)), (2, Fraction(11, 6))]
    1 [(0, Fraction(5, 6)), (1, Fraction(1, 3))]
    2 [(0, Fraction(-1, 6))]    
    -2 [(2, 1/3)]
    """
    rj_set = set()
    rj_dict = {}
    for i in range(rows):
        r = i
        print(f'r={r}',end=':')
        for j in range(cols):
            rj = - r + j
            rj_set.add(rj)
            print(f'{rj}',end=' ')
            rj_dict.setdefault(rj, []).append((r,mymat[i][j]))
        print()
    print(f'rj_set={rj_set}')
    sorted_rj_set = sorted(rj_set)
    print(f'sorted_rj_set={sorted_rj_set}')
    #print(f'rj_dict={rj_dict}')
    #print(sorted(rj_dict.items()))
       
    for rj, pairs in sorted(rj_dict.items()):
        # 将每个元组格式化为 (idx, 分数) 字符串
        pair_str = ", ".join(f"({idx}, {frac})" for idx, frac in pairs)
        print(f"{rj} [{pair_str}]")
        

    
# i-2 i-1, i i+1 i+2
# vi+1/2(-),r=sum crl*vi-r+l,l=1,kk-1;
# kk=5
# r=0: -r+l=0,1,2,3,4:i,i+1,i+2,i+3,i+4
# r=1: -r+l=-1,0,1,2,3:i-1,i,i+1,i+2,i+3
# r=2: -r+l=-2,-1,0,1,2:i-2,i-1,i,i+1,i+2
        
    
xfrac = Fraction(1,2)
k5=5
mymat5L = calc_coef_formula(k5, xfrac)
print_matrix_fraction(mymat5L)
print_coef_formula(mymat5L,xfrac)

r=2
row_matrixL = mymat5L[r, :].reshape(1, -1)
print_matrix_fraction(row_matrixL)
print_coef_formula(row_matrixL,xfrac,0,r)

mymat5R = calc_coef_formula(k5, -xfrac)
print_matrix_fraction(mymat5R)
print_coef_formula(mymat5R,-xfrac)

row_matrixR = mymat5R[r, :].reshape(1, -1)
print_matrix_fraction(row_matrixR)
print_coef_formula(row_matrixR,-xfrac,0,r)

k3=3
mymat3L = calc_coef_formula(k3, xfrac)
mymat3R = calc_coef_formula(k3, -xfrac)

print_matrix_fraction(mymat3L)
print_coef_formula(mymat3L,xfrac)

print_matrix_fraction(mymat3R)
print_coef_formula(mymat3R,-xfrac)

print_matrix_fraction(row_matrixL)
print_matrix_fraction(mymat3L)

printhhh(row_matrixL, mymat3L)
