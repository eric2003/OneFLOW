import numpy as np
import math
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
    print()

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
    
    
def build_moment_matrix(template_index: int, stencil_width: int) -> np.ndarray:
    r"""
    Build the moment matrix M for a given substencil, where
    
        M @ poly_coeffs = cell_averages
    
    The substencil corresponding to `template_index = r` uses the cells:
    
        $$
        I_{i - r},\ I_{i - r + 1},\ \dots,\ I_{i - r + k - 1}
        $$
    
    with $k = \text{stencil\_width}$. Each cell $I_j$ is the interval $[j - 1/2, j + 1/2]$.

    The matrix entry M[m, i] is the integral of the monomial $\xi^i$ over the m-th cell 
    in the substencil (i.e., over $I_{j_m}$ where $j_m = i - r + m$):

        $$
        M[m, i] = \int_{j_m - 1/2}^{j_m + 1/2} \xi^i \, d\xi
        $$

    Parameters
    ----------
    template_index : int
        Index of the substencil (r = 0, 1, ..., k-1). Larger values shift the stencil left.
    stencil_width : int
        Number of cells in the substencil (k).

    Returns
    -------
    M : np.ndarray of shape (k, k)
        Moment matrix with exact fractional entries.
    """
    rows = []
    for m in range(stencil_width):
        # Spatial index of the m-th cell in the substencil: j = i - r + m
        j = -template_index + m
        left = Fraction(j) - Fraction(1, 2)
        right = Fraction(j) + Fraction(1, 2)
        row = []
        for i in range(stencil_width):
            val = integral_xi(right, i) - integral_xi(left, i)
            row.append(val)
        rows.append(row)
    return np.array(rows, dtype=object)
    
def compute_stencil_coefficients_for_point(
    template_index: int,
    stencil_width: int,
    x_point: Fraction
) -> np.ndarray:    
    r"""
    Compute the reconstruction coefficients for a single substencil used to approximate 
    the point value at `x_point` (e.g., $x = i + 1/2$) from cell averages.

    The substencil corresponding to `template_index = r` (where $r = 0, 1, ..., k-1$)
    uses the following $k = \text{stencil\_width}$ consecutive cells:
    
        $$
        I_{i - r},\ I_{i - r + 1},\ \dots,\ I_{i - r + k - 1}
        $$

    For example, when `stencil_width = 3` and reconstructing $v_{i+1/2}^-$:
        - `template_index = 0` → cells [i,   i+1, i+2]  (rightmost)
        - `template_index = 1` → cells [i-1, i,   i+1]  (middle)
        - `template_index = 2` → cells [i-2, i-1, i  ]  (leftmost)

    The returned coefficients `c[0], c[1], ..., c[k-1]` satisfy:
        $$
        p(x_{\text{point}}) = \sum_{j=0}^{k-1} c[j] \cdot \bar{v}_{i - r + j}
        $$
    where $p(\cdot)$ is the unique polynomial of degree ≤ k−1 that matches the 
    cell averages over the substencil.

    Parameters
    ----------
    template_index : int
        Index of the substencil (0 ≤ template_index < stencil_width).
        Larger values shift the stencil further to the left.
    stencil_width : int
        Number of cells in the substencil (order of accuracy = stencil_width).
    x_point : Fraction
        Relative coordinate where the point value is reconstructed, 
        e.g., Fraction(1, 2) for $i + 1/2$.

    Returns
    -------
    coefficients : np.ndarray of shape (stencil_width,)
        Reconstruction coefficients for the cell averages in the substencil,
        ordered from leftmost to rightmost cell in the stencil.
    """
    
    M = build_moment_matrix(template_index, stencil_width)
    M_inv = inverse_matrix(M)
    monomials = np.array([x_point ** i for i in range(stencil_width)], dtype=object)
    coefficients = monomials @ M_inv
    return coefficients
    
def compute_optimal_reconstruction_stencil(
    stencil_width: int,
    x_point: Fraction
) -> np.ndarray:
    """
    Compute the optimal (high-order) reconstruction stencil centered at cell i,
    using `stencil_width` consecutive cells symmetric around i.

    The stencil covers cells: [i - (k-1)//2, ..., i, ..., i + (k-1)//2]
    and reconstructs the point value at x = i + x_point.

    Example:
        k=5, x_point=1/2 → cells [i-2, i-1, i, i+1, i+2]
        Returns coefficients [c_{-2}, c_{-1}, c_0, c_1, c_2]
    """
    if stencil_width % 2 == 0:
        raise ValueError("Optimal stencil requires odd stencil_width for symmetry.")
    
    r = stencil_width // 2
    
    coefficients = compute_stencil_coefficients_for_point(r, stencil_width, x_point)
    return coefficients

def generate_weno_substencils(stencil_width: int, x_point: Fraction) -> np.ndarray:
    """
    Generate all k = stencil_width substencils for reconstructing a point value at x_point.
    
    The returned matrix has shape (k, k), where:
      - Row r corresponds to the substencil that uses cells:
          [I_{i - r}, I_{i - r + 1}, ..., I_{i - r + k - 1}]
        which is the r-th candidate stencil counting from the RIGHTMOST (r=0) 
        to the LEFTMOST (r=k-1) stencil.
    
    For example, when k=3 and reconstructing v_{i+1/2}^-:
        r=0 → cells [i,   i+1, i+2]  (rightmost)
        r=1 → cells [i-1, i,   i+1]  (middle)
        r=2 → cells [i-2, i-1, i  ]  (leftmost)
    """
    
    stencils = []
    for r in range(stencil_width):
        # r = 0 → rightmost stencil
        # r = stencil_width-1 → leftmost stencil
        
        coef = compute_stencil_coefficients_for_point(r, stencil_width, x_point)
        stencils.append(coef)
    return np.vstack(stencils)
    
def generate_left_stencils(stencil_width: int, offset: Fraction = Fraction(1, 2)):
    """生成左偏模板（用于 vi+1/2）"""
    return generate_weno_substencils(stencil_width, offset)

def generate_right_stencils(stencil_width: int, offset: Fraction = Fraction(1, 2)):
    """生成右偏模板（用于 vi-1/2）"""
    return generate_weno_substencils(stencil_width, -offset)    

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
    
def build_substencil_offset_map(sub_stencils):
    """为每个空间偏移 k，记录 (模板索引 r, 系数)"""
    rows, cols = sub_stencils.shape
    offset_map = defaultdict(list)
    for r in range(rows):
        for j in range(cols):
            k = j - r  # spatial offset: v[i + k]
            coef = sub_stencils[r, j]
            offset_map[k].append((r, coef))
    return offset_map
    
def build_target_offset_map(target_row):
    """
    target_row: 1D array like [1/30, -13/60, 47/60, 9/20, -1/20]
    assumes it corresponds to offsets [-2, -1, 0, 1, 2]
    """
    n = len(target_row)
    base_offset = - (n//2)
    offsets = list(range(base_offset, base_offset + n))  # [-2,-1,0,1,2]
    return {k: target_row[i] for i, k in enumerate(offsets)}  
    
def build_linear_system(sub_stencils, target_offset_map):
    """
    Build A x = b for WENO weights.
    
    Returns:
        A: np.ndarray of shape (num_equations, num_templates)
        b: np.ndarray of shape (num_equations,)
        offsets: list of spatial offsets (for labeling)
    """
    sub_offset_map = build_substencil_offset_map(sub_stencils)
    num_templates = sub_stencils.shape[0]
    
    # Get all spatial offsets that appear in target
    offsets = sorted(target_offset_map.keys())
    
    A = []
    b = []
    
    for k in offsets:
        row = [Fraction(0) for _ in range(num_templates)]
        for r, coef in sub_offset_map.get(k, []):
            row[r] = coef
        A.append(row)
        b.append(target_offset_map[k])
    
    # Convert to float for numpy (or keep as Fraction for exact solve)
    A_float = np.array([[float(x) for x in row] for row in A])
    b_float = np.array([float(x) for x in b])
    
    return A_float, b_float, offsets

def solve_weno_weights(sub_stencils, target_offset_map):
    A, b, offsets = build_linear_system(sub_stencils, target_offset_map)
    # Solve Ax = b in least-squares sense
    x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
    
    print("Solved WENO weights:")
    for i, wi in enumerate(x):
        print(f"d[{i}] = {wi:.6f} ≈ {Fraction(wi).limit_denominator(100)}")
    
    # Verify residual
    if len(residuals) > 0:
        print(f"Residual norm: {np.sqrt(residuals[0]):.2e}")
    else:
        # Exact solution (rank-deficient or square)
        residual = np.linalg.norm(A @ x - b)
        print(f"Residual norm: {residual:.2e}")
    
    return x

def print_weno_equations(sub_stencils, target_dict):
    offset_map = build_substencil_offset_map(sub_stencils)
    all_offsets = sorted(target_dict.keys())
    
    rows, cols = sub_stencils.shape
    
    weights = ", ".join(f"d[{i}]" for i in range(rows))
    print(f"WENO linear system (for weights {weights}):\n")
    
    for k in all_offsets:
        terms = []
        for r, coef in offset_map.get(k, []):
            terms.append(f"d[{r}] * ({coef})")
        
        lhs = " + ".join(terms) if terms else "0"
        rhs = target_dict[k]
        print(f"v[i{k:+}] : {lhs} = {rhs}")
    print()
    
def compute_weno_linear_weights(row_matrix, mymat):
    sub_stencils = mymat
    
    # Build target map
    target_dict = build_target_offset_map(row_matrix)
    print_weno_equations(sub_stencils, target_dict)
    
    # Solve
    weights = solve_weno_weights(mymat, target_dict)
    
def compute_weno_linear_weights_new(order):
    xfrac = Fraction(1,2)

    k = order
    kh = 2*k - 1
    
    mymatL = generate_left_stencils(k)
    row_matL = compute_optimal_reconstruction_stencil(kh, xfrac)
    compute_weno_linear_weights(row_matL, mymatL)
    
    mymatR = generate_right_stencils(k)
    row_matR = compute_optimal_reconstruction_stencil(kh, -xfrac)
    compute_weno_linear_weights(row_matR, mymatR)
    
def solve_weno_linear_weights(optimal_stencil: np.ndarray, sub_stencils: np.ndarray) -> np.ndarray:
    """
    Solve for linear weights d such that:
        optimal_stencil ≈ sum_j d[j] * sub_stencils[j]
    
    Prints the linear system and solved weights.
    """
   
    # Build target map
    target_dict = build_target_offset_map(optimal_stencil)
    print_weno_equations(sub_stencils, target_dict)
    
    # Solve
    weights = solve_weno_weights(sub_stencils, target_dict)    
    return weights    
    
def demo_weno_linear_weights(weno_r: int):
    """
    Demonstrate linear weight computation for WENO-r scheme.
    
    Parameters:
        weno_r (int): Number of substencils (e.g., 3 for WENO5, 2 for WENO3)
    """
    x_half = Fraction(1, 2)
    global_stencil_width = 2 * weno_r - 1  # e.g., 5 for WENO3

    # Left-biased (v_{i+1/2}^-)
    substencils_L = generate_weno_substencils(stencil_width=weno_r, x_point=x_half)
    optimal_L = compute_optimal_reconstruction_stencil(
        stencil_width=global_stencil_width, x_point=x_half
    )
    weights_L = solve_weno_linear_weights(optimal_L, substencils_L)

    # Right-biased (v_{i-1/2}^+)
    substencils_R = generate_weno_substencils(stencil_width=weno_r, x_point=-x_half)
    optimal_R = compute_optimal_reconstruction_stencil(
        stencil_width=global_stencil_width, x_point=-x_half
    )
    weights_R = solve_weno_linear_weights(optimal_R, substencils_R)

    return weights_L, weights_R
    
def demo_weno_linear_weights_maxk():
    maxk = 3
    for k in range(1,maxk+1):
        print(f"\n=== WENO{2*k-1} ===")
        demo_weno_linear_weights(weno_r=k)
        
def derivative_form(n, m):
    """
    返回 x^n 的 m 阶导数形式 (系数, x的指数)
    示例: derivative_form(3, 2) → (6, 1) 表示 6x^1
    """
    if m > n:
        return (0, 0)  # 或返回 None
    
    # 计算系数 n!/(n-m)!
    coeff = math.factorial(n) // math.factorial(n - m)
    power = n - m
    
    return (coeff, power)
    
def print_power_symbol(power_map):
    sorted_keys = sorted(power_map)
    # 遍历所有键，但只在不是最后一项时打印 "+"
    #print(f"len(sorted_keys)={len(sorted_keys)}")
    for idx, key in enumerate(sorted_keys):
        mylist = power_map[key]
        n = len( mylist )
        print(f"(", end = '')
        for i in range(n):
            coef, acoef = mylist[i]
            if i == n-1:
                print(f"{coef}*a{acoef}", end = '')
            else:
                print(f"{coef}*a{acoef} + ", end = '')
        # 判断是否是最后一项（关键修改）
        print(f")*x^{key}",end = '')
        if idx < len(sorted_keys) - 1:
            print(" + ",end = '')  # 不是最后一项，打印 "+"
        else:
            print()     # 是最后一项，只换行不打印 "+"
            
def demo_smoothness_indicator():
    print(f'demo_smoothness_indicator')
   
    n = 5
    m = 2
    coeff, power = derivative_form(5, 2)
    print(f"d^{{{m}}}/dx^{{{m}}}(x^({n}))={coeff}x^{power}")
    
    k = 3
    rows = k-1
    cols = k-1
    matrix = np.empty((rows, cols), dtype=object)
    #print(f'matrix=\n{matrix}')
    
    #x^1 x^2 x^3
    #d^1dx^1 1x^0 2x^1 3x^2
    #d^2dx^2 0x^0 2x^0 6x^1
    #d^3dx^3 0x^0 0x^0 6x^0
    for i in range(rows):
        for j in range(cols):
            coef, power = derivative_form(j+1, i+1)
            acoef = j + 1
            matrix[i][j] = (coef, acoef, power)
            print(f"{coeff}x^{power}",end=' ')
        print()
        
    print(f'matrix=\n{matrix}')
    power_map_list = []
    for i in range(rows):
        power_map = defaultdict(list)
        for j in range(cols):
            coef, acoef, power = matrix[i][j]
            if coef != 0:
                power_map[power].append((coef, acoef))
            print(f"{coef}*a{acoef}*x^{power}",end=' ')
        power_map_list.append(power_map)
        print()
    
    print(f'power_map_list={power_map_list}')
    for i in range(rows):
        print_power_symbol(power_map_list[i])
        
def square_polynomial(power_map):
    """
    多项式平方展开
    参数:
        power_map: {指数: [(数值系数, 符号索引), ...]}
    返回:
        展开后的结果字典
    """
    result = defaultdict(list)
    sorted_exps = sorted(power_map.keys())
    
    # 1. 计算每个项的平方 (a^2)
    for exp in sorted_exps:
        terms = power_map[exp]
        new_exp = exp * 2
        
        for coef, acoef in terms:
            # 存储为 (数值系数, 符号索引, 是否平方标志)
            result[new_exp].append((coef, acoef, True))
    
    # 2. 计算交叉项 (2ab)
    for i in range(len(sorted_exps)):
        for j in range(i + 1, len(sorted_exps)):
            exp_i, exp_j = sorted_exps[i], sorted_exps[j]
            new_exp = exp_i + exp_j
            
            for coef_i, acoef_i in power_map[exp_i]:
                for coef_j, acoef_j in power_map[exp_j]:
                    # 交叉项系数乘以2
                    result[new_exp].append((2 * coef_i * coef_j, 
                                          (acoef_i, acoef_j), False))
    
    return result        
        
def merge_and_simplify(power_map):
    """合并同类项并化简符号乘积"""
    simplified = defaultdict(dict)  # exp: {符号键: 总系数}
    
    for exp, terms in power_map.items():
        for term_info in terms:
            # 解析 term_info: (系数, 符号信息, 是否平方项)
            if len(term_info) == 3:
                coef, symbols_part, is_square = term_info
            else:
                # 兼容旧数据格式
                coef, symbols_part = term_info
                is_square = False
            
            # 关键修复：根据标志计算有效系数
            if is_square:
                # 平方项：系数需要平方 (c*a)^2 → c^2 * a^2
                effective_coef = coef * coef
            else:
                # 交叉项：系数已是最终值 2*c_i*c_j
                effective_coef = coef
            
            # 生成标准化符号键
            if isinstance(symbols_part, tuple):
                # 交叉项：(索引1, 索引2)
                i, j = symbols_part
                symbol_names = [f"a{i}", f"a{j}"]
            else:
                # 平方项：单个符号索引
                symbol_names = [f"a{symbols_part}", f"a{symbols_part}"]
            
            # 排序保证 a1*a2 和 a2*a1 被视为相同
            key = tuple(sorted(symbol_names))
            
            # 累加系数
            simplified[exp][key] = simplified[exp].get(key, 0) + effective_coef
    
    # 转换回列表格式
    result = defaultdict(list)
    for exp, term_dict in simplified.items():
        for symbol_key, total_coef in term_dict.items():
            result[exp].append((total_coef, symbol_key))
    
    return result
    
def format_term(coef, symbols):
    """
    格式化单项式，支持多种输入类型
    symbols 可能是：
      - 整数 (如 1) → 表示 a1
      - 字符串 (如 'a1') → 直接使用
      - 元组 (如 ('a1', 'a2')) → a1*a2 或 a1^2
    """
    # 处理 symbols 为整数的情况（如 1 → 'a1'）
    if isinstance(symbols, int):
        symbol_str = f"a{symbols}"
    # 处理 symbols 为字符串的情况
    elif isinstance(symbols, str):
        symbol_str = symbols
    # 处理 symbols 为元组/列表的情况
    elif isinstance(symbols, (tuple, list)):
        if len(symbols) == 1:
            symbol_str = str(symbols[0])
        elif symbols[0] == symbols[1]:  # 平方项
            symbol_str = f"{symbols[0]}^2"
        else:  # 不同符号相乘
            symbol_str = "*".join(str(s) for s in symbols)
    else:
        raise TypeError(f"Unsupported symbols type: {type(symbols)}")
    
    return f"{coef}*{symbol_str}"
    
def print_expanded_polynomial(power_map):
    """
    打印展开后的多项式，不显示最后一行的"+"
    """
    if not power_map:
        print("0")
        return
    
    sorted_keys = sorted(power_map.keys())
    
    for idx, key in enumerate(sorted_keys):
        terms = power_map[key]
        
        # 构建该项的内部表达式
        inner_terms = []
        for coef, symbols in terms:
            inner_terms.append(format_term(coef, symbols))
        
        # 打印 (内层表达式)*x^key
        if len(inner_terms) == 1:
            print(f"({inner_terms[0]})*x^{key}", end='')
        else:
            print(f"({' + '.join(inner_terms)})*x^{key}", end='')
        
        # 判断是否是最后一项
        if idx < len(sorted_keys) - 1:
            print(" + ", end='')
        else:
            print()  # 最后一项只换行        

def square_polynomial_test():
    # ============= 测试代码 =============

    # 测试1: (1*a1)*x^0 + (2*a2)*x^1
    poly1 = {
        0: [(1, 1)],  # x^0 项: 1*a1
        1: [(2, 2)]   # x^1 项: 2*a2
    }

    print("原始多项式1:")
    print_expanded_polynomial(poly1)
    # 输出: (1*a1)*x^0 + (2*a2)*x^1

    print("\n平方展开后:")
    expanded1 = square_polynomial(poly1)
    print(f"expanded1={expanded1}")
    merged1 = merge_and_simplify(expanded1)
    print(f"merged1={merged1}")
    print_expanded_polynomial(merged1)
    # 期望: (1*a1^2)*x^0 + (4*a1*a2)*x^1 + (4*a2^2)*x^2

    # 测试2: (2*a2)*x^0
    poly2 = {
        0: [(2, 2)]  # x^0 项: 2*a2
    }

    print("\n原始多项式2:")
    print_expanded_polynomial(poly2)
    # 输出: (2*a2)*x^0

    print("\n平方展开后:")
    expanded2 = square_polynomial(poly2)
    merged2 = merge_and_simplify(expanded2)
    print_expanded_polynomial(merged2)
    # 期望: (4*a2^2)*x^0
    
if __name__ == "__main__":
    #demo_smoothness_indicator()
    square_polynomial_test()



