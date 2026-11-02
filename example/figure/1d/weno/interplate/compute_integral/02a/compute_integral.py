import numpy as np
from fractions import Fraction
from collections import defaultdict
from typing import List, Tuple, Dict

# ============ 基础函数（保持不变） ============

def compute_integral(alpha: float, beta: float, power: int) -> float:
    """计算∫_{α}^{β} ξ^power dξ"""
    if power == 0:
        return beta - alpha
    return (beta**(power + 1) - alpha**(power + 1)) / (power + 1)

def compute_alpha_beta(row_index: int, r: int) -> Tuple[float, float]:
    """计算第row_index行的积分区间[α, β]"""
    middle = -r + row_index
    return middle - 0.5, middle + 0.5

def compute_matrix_M(k: int, r: int) -> np.ndarray:
    """计算k×k的矩阵M（数值矩阵）"""
    M = np.zeros((k, k), dtype=float)
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        for j in range(k):
            M[i, j] = compute_integral(alpha, beta, j)
    return M

def solve_for_coefficients(k: int, r: int, i: int = 0) -> List[str]:
    """生成a系数的v表达式（已修正符号连接）"""
    M = compute_matrix_M(k, r)
    try:
        M_inv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        raise ValueError(f"矩阵M (k={k}, r={r})不可逆！")
    
    a_coeffs = []
    for j in range(k):
        terms = []
        for m in range(k):
            coeff = M_inv[j, m]
            if abs(coeff) > 1e-10:
                offset = m - r
                v_str = f"v[i{offset:+d}]" if offset != 0 else "v[i]"
                
                if abs(coeff - 1.0) < 1e-10:
                    term_str = v_str
                elif abs(coeff + 1.0) < 1e-10:
                    term_str = f"-{v_str}"
                else:
                    frac = Fraction(coeff).limit_denominator(1000)
                    term_str = f"{frac}*{v_str}"
                
                terms.append(term_str)
        
        # 智能连接
        if not terms:
            a_coeffs.append("0")
        else:
            expr = terms[0]
            for term in terms[1:]:
                if term.startswith('-'):
                    expr += f" {term}"
                else:
                    expr += f" + {term}"
            a_coeffs.append(expr)
    
    return a_coeffs

# ============ 核心改进：Fraction格式化 ============

def format_coefficient(frac: Fraction, force_display: bool = False) -> str:
    """
    智能格式化Fraction为字符串
    
    参数:
        frac: Fraction对象
        force_display: 是否强制显示系数（即使为1）
    
    返回:
        格式化后的系数字符串
    """
    # 处理零
    if frac == 0:
        return "0"
    
    # 整数系数
    if frac.denominator == 1:
        value = int(frac.numerator)
        if value == 1 and not force_display:
            return ""          # 系数1，省略
        elif value == -1:
            return "-"         # 系数-1，只返回负号
        else:
            return str(value)  # 整数系数，如"3"
    
    # 分数系数
    return f"{frac.numerator}/{frac.denominator}"

def evaluate_polynomial_integral_symbolic(polynomial: Dict[int, List[Tuple[float, List[int]]]], 
                                         a_coeffs: List[str]) -> str:
    """
    对多项式进行积分，并将a_j替换为v表达式
    
    关键修复：全程使用Fraction，避免浮点数
    """
    # ✅ 修复：使用Fraction表示积分限
    a = Fraction(-1, 2)  # -1/2
    b = Fraction(1, 2)   # 1/2
    
    # 用字典累积符号项的系数（使用Fraction）
    result_dict = defaultdict(lambda: Fraction(0, 1))
    
    # 步骤1：积分并累加系数（使用Fraction）
    for exp, expr_list in polynomial.items():
        # 计算积分因子：(b^(exp+1) - a^(exp+1))/(exp+1)
        # 分子和分母都是Fraction
        numerator = b**(exp + 1) - a**(exp + 1)  # Fraction类型
        integral_factor = Fraction(numerator, exp + 1)  # 构造函数接收Fraction和int
        
        for coeff, symbols in expr_list:
            # 生成符号键
            if len(symbols) == 1:
                symbol_key = f"a{symbols[0]}"
            elif symbols[0] == symbols[1]:
                symbol_key = f"a{symbols[0]}^2"
            else:
                symbol_key = "*".join([f"a{s}" for s in symbols])
            
            # 累加分数系数
            # Fraction(coeff)将float转为Fraction（近似）
            # 更精确的做法是：直接传入Fraction(str(coeff))
            contribution = Fraction(str(coeff)) * integral_factor
            result_dict[symbol_key] += contribution
    
    # 步骤2：构建复合表达式
    terms = []
    for symbol_key, total_frac in result_dict.items():
        if total_frac == 0:
            continue
        
        # 获取对应的a_j表达式
        base_symbol = symbol_key.split('*')[0].split('^')[0]
        a_index = int(base_symbol[1:]) - 1
        a_expr = a_coeffs[a_index]
        
        # 格式化系数（接收Fraction）
        coeff_str = format_coefficient(total_frac)
        
        # 构建项
        if '^2' in symbol_key:
            # 平方项：需要括号
            if coeff_str:
                term = f"{coeff_str}*({a_expr})^2"
            else:
                term = f"({a_expr})^2"
        else:
            # 一次项：直接连接
            if coeff_str:
                term = f"{coeff_str}*{a_expr}"
            else:
                term = f"{a_expr}"
        
        terms.append(term)
    
    # 步骤3：智能连接各项
    if not terms:
        return "0"
    
    # 找到第一个非负项
    first_idx = 0
    while first_idx < len(terms) and terms[first_idx].startswith('-'):
        first_idx += 1
    
    if first_idx < len(terms):
        expr = terms[first_idx]
        for t in terms[:first_idx]:
            expr = f"{t} + {expr}"
        for t in terms[first_idx+1:]:
            if t.startswith('-'):
                expr += f" {t}"
            else:
                expr += f" + {t}"
    else:
        expr = terms[0]
        for t in terms[1:]:
            expr += f" {t}"
    
    return expr

# ============ 测试函数 ============

def generate_composite_expressions(k: int, polynomial: Dict[int, List[Tuple[float, List[int]]]], i: int = 0):
    """生成所有r值的复合表达式f_r（Fraction优化版）"""
    f_r_dict = {}
    
    print(f"\n生成k={k}的复合表达式f_r（Fraction优化版）")
    print("="*70)
    
    for r in range(k):
        print(f"\n--- r = {r} ---")
        
        # 1. 生成a系数的v表达式
        a_coeffs = solve_for_coefficients(k, r, i)
        print("a系数的v表达式:")
        for idx, expr in enumerate(a_coeffs):
            print(f"  a_{idx} = {expr}")
        
        # 2. 生成复合表达式f_r
        f_r = evaluate_polynomial_integral_symbolic(polynomial, a_coeffs)
        f_r_dict[r] = f_r
        
        print(f"\nf_{r}（Fraction格式）:")
        print("-" * 50)
        # 美化输出：分行显示
        if '+' in f_r:
            lines = f_r.split('+')
            for line in lines:
                print(f"  {line.strip()}")
        else:
            print(f"  {f_r}")
    
    # 总结
    print("\n" + "="*70)
    print("总结：所有r的复合表达式f_r（Fraction格式）")
    print("="*70)
    for r, expr in f_r_dict.items():
        print(f"\nf_{r} = {expr}")
    
    return f_r_dict

def test_composite():
    """测试你的例子"""
    k = 3
    
    # 多项式：a1^2 + 4*a1*a2*x + (4+4)*a2^2*x^2
    polynomial = {
        0: [(1.0, [1, 1]),   # a1^2
            (4.0, [2, 2])],  # a2^2（系数4）
        1: [(4.0, [1, 2])],  # a1*a2
        2: [(4.0, [2, 2])]   # a2^2（系数4）
    }
    
    print("="*70)
    print("测试：多项式积分后复合a系数表达式（Fraction优化）")
    print("="*70)
    
    print("\n原始多项式（用a_j表示）:")
    for exp, terms in polynomial.items():
        term_strs = []
        for coeff, symbols in terms:
            if len(symbols) == 1:
                term_strs.append(f"{coeff}*a{symbols[0]}")
            elif symbols[0] == symbols[1]:
                term_strs.append(f"{coeff}*a{symbols[0]}^2")
            else:
                term_strs.append(f"{coeff}*a{symbols[0]}*a{symbols[1]}")
        print(f"  x^{exp}: {' + '.join(term_strs)}")
    
    # 生成复合表达式
    f_r_dict = generate_composite_expressions(k, polynomial, i=0)

if __name__ == "__main__":
    test_composite()