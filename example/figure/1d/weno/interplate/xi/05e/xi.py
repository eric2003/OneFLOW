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
    
    
def print_formula(xpcoef,vv):
    #print(f'print_formula vv={vv}')
    ii = 0
    strr = ''
    for v in vv:
        absv = abs(v)
        s = ''
        t = str(absv)
        if v > 0:
            s='+'
        elif v < 0:
            s='-'
        else:
            t = ''
        var1 = xpcoef[ii]
        absv1 = abs(var1)
        ff = r'\cfrac{{{absv1.numerator}}}{{{absv1.denominator}}}'
        if var1 >= 0:
            sf = '+'
        else:
            sf = '-'
        
        if ii == 0 and var1 >= 0:
            fv1 = ff
        else:
            fv1 = sf + ' ' + str(ff)

        var = r'\overline{{v}}_{{i{s}{t}}}'
        strr += f'{fv1}' + var
        ii += 1
            
    print(strr)

        
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
    
#vv = [-1,0,1]
vv = [-2,-1,0,1,2]
k = len(vv)

# 测试一维向量xx的打印（支持行向量/列向量两种格式）
xp = compute_coef(Fraction(1,2), k)
xn = compute_coef(-Fraction(1,2), k)

"""
 r | j=0  j=1  j=2
-1 | 11/6 -7/6  1/3 (i+1,i+2,i+3)
 0 |  1/3  5/6 -1/6 (i  ,i+1,i+2)
 1 | -1/6  5/6  1/3 (i-1,i  ,i+1)
 2 |  1/3 -7/6 11/6 (i-2,i-1,i  )
"""
 
a_str  = ["1/30", "-13/60", "47/60", "9/20", "-1/20"]
a0_str = ["1/3", "5/6", "-1/6"]
a1_str = ["-1/6", "5/6", "1/3"]
a2_str = ["1/3", "-7/6", "11/6"]
a  = [Fraction(s) for s in a_str]
a0 = [Fraction(s) for s in a0_str]
a1 = [Fraction(s) for s in a1_str]
a2 = [Fraction(s) for s in a2_str]
print_matrix_fraction(a)
print_matrix_fraction(a0)
print_matrix_fraction(a1)
print_matrix_fraction(a2)

alist = []
alist.append(a0)
alist.append(a1)
alist.append(a2)

N = len(a)
M = len(a0)
frac_list = [Fraction(0) for _ in range(N)]
print_matrix_fraction(frac_list)

im = 2
for r in range(M):
    print( f'r={r}:', end='')
    for j in range(M):
        id = im - r + j
        print( id, end=' ')
    print()

for r in range(M):
    print( f'r={r}:', end='')
    for j in range(M):
        id = im - r + j
        print( alist[r][j], end=' ')
    print()

dd = [Fraction(0) for _ in range(M)]
dd[0] = Fraction(3,10)
dd[1] = Fraction(6,10)
dd[2] = Fraction(1,10)
    
for r in range(M):
    for j in range(M):
        id = im - r + j
        frac_list[id] += dd[r]*alist[r][j]


print_matrix_fraction(frac_list)

rows=M
cols=M
print_matrix_fraction(alist)

widths = np.empty(M, dtype=int)

for j in range(rows):
    w = 0
    for i in range(cols):
        ww = len( coef_tostring(alist[i][j],j) )
        w = max(w, ww)
    widths[j] = w
    
print(f'widths={widths}')

for i in range(rows):
    mystr = ''
    r = i
    for j in range(cols):
        rj = -r + j
        var_rj = id_tostring(rj)
        mystr += coef_tostring(alist[i][j],j)+ f"*v[i{var_rj}]"
    print(f'mystr={mystr}')
    
for i in range(rows):
    mystr = ''
    r = i
    for j in range(cols):
        rj = -r + j
        var_rj = id_tostring(rj)
        ss = f"{coef_tostring(alist[i][j],j):>{widths[j]}}"
        mystr += ss+ f"*v[i{var_rj}]"
    print(f'mystr={mystr}')    
