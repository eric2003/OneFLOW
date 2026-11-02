import numpy as np
from fractions import Fraction
from collections import defaultdict

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
    
def format_signed_coef(coef,isFirstElement,width):
    """
    vi+1/2(-),0= 1/3*v[i  ]+5/6*v[i+1]- 1/6*v[i+2]
    vi+1/2(-),1=-1/6*v[i-1]+5/6*v[i  ]+ 1/3*v[i+1]
    vi+1/2(-),2= 1/3*v[i-2]-7/6*v[i-1]+11/6*v[i  ]
    """
    abscoef,sign = coef_toabsstring( coef )
    if isFirstElement and sign == '+':
        sign = ' '
    signed_coef_str = f"{sign}{abscoef:>{width-1}}"
    return signed_coef_str
    
def get_sign_and_abs_str(coef):
    """将系数拆分为符号字符（'+', '-', ' '）和绝对值字符串。"""
    if coef >= 0:
        return '+', str(coef)
    else:
        return '-', str(-coef)    
    
def compute_column_widths(coef_matrix):
    """计算每列系数显示所需的最大宽度（含符号位）"""
    rows, cols = coef_matrix.shape
    widths = np.empty(cols, dtype=int)
    for j in range(cols):
        max_width = 0
        for i in range(rows):
            _, abs_str = get_sign_and_abs_str(coef_matrix[i, j])
            width = len(abs_str) + 1  # +1 for sign or space
            max_width = max(max_width, width)
        widths[j] = max_width
    return widths    

def build_variable_index_string(offset):
    """将偏移量转为 '[i+2]', '[i-1]', '[i  ]' 等字符串"""
    if offset == 0:
        return "[i  ]"
    elif offset > 0:
        return f"[i+{offset}]"
    else:
        return f"[i{offset}]"  # offset already includes minus, e.g., -2 → [i-2]
        
def build_variable_indices(start_offset, num_cols):
    """生成每列对应的变量索引字符串列表"""
    return [build_variable_index_string(start_offset + j) for j in range(num_cols)]
    
def build_lhs_label(x_frac: Fraction, shift: int = 0) -> str:
    # 1. 计算总偏移 = shift + x_frac
    total_offset = shift + x_frac

    # 2. 格式化总偏移字符串（带符号）
    if total_offset >= 0:
        offset_str = f"+{total_offset}"  # e.g., +1/2, +3/2
    else:
        offset_str = f"{total_offset}"   # e.g., -1/2（Fraction 会自动带负号）

    # 3. 确定方向标志：由原始 x_frac 的符号决定（非 total_offset！）
    direction = '-' if x_frac >= 0 else '+'

    # 4. 拼接
    return f"vi{offset_str}({direction})"
    
def format_stencil_row(row, jstart, cols, widths):
    terms = []
    ioffset_strs = build_variable_indices(jstart, cols)
    
    for j in range(cols):
        term_str = format_signed_coef(row[j],j==0,widths[j])
        terms.append(f"{term_str}*v{ioffset_strs[j]}")
        
    rhs_label = ''.join(terms)
    return rhs_label   

def print_stencil_formula(coef_matrix,xfrac,ishift=0,base_row=0):
    rows, cols = coef_matrix.shape
    widths = compute_column_widths(coef_matrix)
    
    lhs_label = build_lhs_label(xfrac, ishift)
    
    for i in range(rows):
        r = base_row if rows == 1 else i
        row = coef_matrix[i]
        jstart = ishift - r
        rhs_label = format_stencil_row(row, jstart, cols, widths)
        print(f'{lhs_label},{r}={rhs_label}')
        
    print()
    
    
def printhhh(row_matrix, mymat, rbase, xfrac):
    rows_ref, cols_ref = row_matrix.shape
    print(f'rows_ref,cols_ref={rows_ref},{cols_ref}')
    rows, cols = mymat.shape
    print(f'rows,cols={rows},{cols}')
    
    print_matrix_fraction(row_matrix)
    print_matrix_fraction(mymat)
    
    print(f'rbase={rbase}')
    
    print_stencil_formula(row_matrix,xfrac,0,rbase)
   
    print_stencil_formula(mymat,xfrac)    
    
    """
    rows_ref,cols_ref=1,5
    rows,cols=3,3
    [ 1/30, -13/60, 47/60, 9/20, -1/20 ]
    [  1/3,  5/6, -1/6 ]
    [ -1/6,  5/6,  1/3 ]
    [  1/3, -7/6, 11/6 ]
    rbase=2
    vi+1/2(-),2= 1/30*v[i-2]-13/60*v[i-1]+47/60*v[i  ]+9/20*v[i+1]-1/20*v[i+2]

    vi+1/2(-),0= 1/3*v[i  ]+5/6*v[i+1]- 1/6*v[i+2]
    vi+1/2(-),1=-1/6*v[i-1]+5/6*v[i  ]+ 1/3*v[i+1]
    vi+1/2(-),2= 1/3*v[i-2]-7/6*v[i-1]+11/6*v[i  ]

    r=0:0 1 2
    r=1:-1 0 1
    r=2:-2 -1 0
    rj_set={0, 1, 2, -2, -1}
    sorted_rj_set=[-2, -1, 0, 1, 2]
    -2 [(2, 1/3)]
    -1 [(1, -1/6), (2, -7/6)]
    0 [(0, 1/3), (1, 5/6), (2, 11/6)]
    1 [(0, 5/6), (1, 1/3)]
    2 [(0, -1/6)]
    vi+1/2(-)=d[0]*( 1/3*v[i  ]+5/6*v[i+1]- 1/6*v[i+2])
             +d[1]*(-1/6*v[i-1]+5/6*v[i  ]+ 1/3*v[i+1])
             +d[2]*( 1/3*v[i-2]-7/6*v[i-1]+11/6*v[i  ])
    -2 [(2, 1/3)] 表示下标-2也就是v[i-2]对应的系数为d[2]*1/3
    这个系数应该和vi+1/2(-),2= 1/30*v[i-2]-13/60*v[i-1]+47/60*v[i  ]+9/20*v[i+1]-1/20*v[i+2]
    里面的v[i-2]系数一致，也就是d[2]*1/3v[i-2]=1/30v[i-2]，从而求出d[2]，同理可以求出d[0]，
    最后求出d[1]，我的想法是一步步搞清楚到底该怎么做，在思维不清晰之前先使用代码将各项表达式对应打印出来，
    第一步-2 [(2, 1/3)]已经做到了，但是不直观。第二步应该将其变化为
    v[i-2]: [(d[2]*1/3)=1/30]
    v[i-1]: [(d[1]*(-1/6)+d[2]*( -7/6)=-13/60],以此类推。你有什么建议，并给出代码以便于进一步讨论。
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
       
    for rj, pairs in sorted(rj_dict.items()):
        # 将每个元组格式化为 (idx, 分数) 字符串
        pair_str = ", ".join(f"({idx}, {frac})" for idx, frac in pairs)
        print(f"{rj} [{pair_str}]")
        
def build_offset_map(sub_stencils):
    """为每个空间偏移 k，记录 (模板索引 r, 系数)"""
    rows, cols = sub_stencils.shape
    offset_map = defaultdict(list)
    for r in range(rows):
        for j in range(cols):
            k = j - r  # spatial offset: v[i + k]
            coef = sub_stencils[r, j]
            offset_map[k].append((r, coef))
    return offset_map        
        
def print_weno_equations(sub_stencils, target_dict):
    offset_map = build_offset_map(sub_stencils)
    all_offsets = sorted(target_dict.keys())

    print("WENO linear system (for weights d[0], d[1], d[2]):\n")
    for k in all_offsets:
        terms = []
        for r, coef in offset_map.get(k, []):
            terms.append(f"d[{r}] * ({coef})")
        
        lhs = " + ".join(terms) if terms else "0"
        rhs = target_dict[k]
        print(f"v[i{k:+}] : {lhs} = {rhs}")
    print()
    
def testfun(mymat):
    sub_stencils = mymat
    target_dict = {
        -2: Fraction(1, 30),
        -1: Fraction(-13, 60),
         0: Fraction(47, 60),
         1: Fraction(9, 20),
         2: Fraction(-1, 20)
    }
    print_weno_equations(sub_stencils, target_dict)
 
    
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
print_stencil_formula(mymat5L,xfrac)

r=2
row_matrixL = mymat5L[r, :].reshape(1, -1)
print_matrix_fraction(row_matrixL)
print_stencil_formula(row_matrixL,xfrac,0,r)

mymat5R = calc_coef_formula(k5, -xfrac)
print_matrix_fraction(mymat5R)
print_stencil_formula(mymat5R,-xfrac)

row_matrixR = mymat5R[r, :].reshape(1, -1)
print_matrix_fraction(row_matrixR)
print_stencil_formula(row_matrixR,-xfrac,0,r)

k3=3
mymat3L = calc_coef_formula(k3, xfrac)
mymat3R = calc_coef_formula(k3, -xfrac)

print_matrix_fraction(mymat3L)
print_stencil_formula(mymat3L,xfrac)

print_matrix_fraction(mymat3R)
print_stencil_formula(mymat3R,-xfrac)

#printhhh(row_matrixL, mymat3L, r, xfrac)

testfun(mymat3L)
