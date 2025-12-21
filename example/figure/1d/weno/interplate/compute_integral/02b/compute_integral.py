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
    """生成a系数的v表达式"""
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

# ============ 核心改进：分离积分与代入 ============

def format_coefficient(frac: Fraction, force_display: bool = False) -> str:
    """智能格式化Fraction为字符串"""
    if frac == 0:
        return "0"
    
    if frac.denominator == 1:
        value = int(frac.numerator)
        if value == 1 and not force_display:
            return ""
        elif value == -1:
            return "-"
        else:
            return str(value)
    
    return f"{frac.numerator}/{frac.denominator}"

def integrate_polynomial_x(polynomial: Dict[int, List[Tuple[float, List[int]]]]) -> List[Tuple[Fraction, List[int]]]:
    """
    ✅ 第一步：对x在[-1/2, 1/2]上积分，消去x变量
    
    数学原理：
    - 对x^0项：∫_{-1/2}^{1/2} coeff * a_j * x^0 dx = coeff * a_j * 1
    - 对x^1项：∫_{-1/2}^{1/2} coeff * a_j * x^1 dx = coeff * a_j * 0
    - 对x^2项：∫_{-1/2}^{1/2} coeff * a_j * x^2 dx = coeff * a_j * (1/12)
    
    参数:
        polynomial: {幂次: [(系数, [a索引,...]), ...]}
                     例如: {0: [(1.0, [1,1]), (4.0, [2,2])], 1: [(4.0, [1,2])], 2: [(4.0, [2,2])]}
    
    返回:
        积分后的项列表: [(分数系数, [a索引,...]), ...]
    """
    # ✅ 积分限: [-1/2, 1/2]
    a = Fraction(-1, 2)
    b = Fraction(1, 2)
    
    integrated_terms = []
    
    print("\n积分过程详细计算:")
    print("-" * 50)
    
    # 遍历多项式的每个幂次项
    for exp, expr_list in polynomial.items():
        # ✅ 计算 ∫x^exp dx 在[-1/2, 1/2]上的值
        numerator = b**(exp + 1) - a**(exp + 1)
        integral_factor = Fraction(numerator, exp + 1)
        
        print(f"  x^{exp} 积分因子: ∫ξ^{exp}dξ = {integral_factor}")
        
        # 对当前幂次下的每个a_j表达式项应用积分
        for coeff, symbols in expr_list:
            # 将浮点系数转为精确分数
            coeff_frac = Fraction(str(coeff))
            
            # 积分后的系数 = 原系数 × 积分因子
            new_coeff = coeff_frac * integral_factor
            
            # 生成符号表示
            if len(symbols) == 1:
                symbol_str = f"a{symbols[0]}"
            elif symbols[0] == symbols[1]:
                symbol_str = f"a{symbols[0]}^2"
            else:
                symbol_str = "*".join([f"a{s}" for s in symbols])
            
            print(f"    {coeff_frac} * {symbol_str} × {integral_factor} = {new_coeff} * {symbol_str}")
            
            integrated_terms.append((new_coeff, symbols))
    
    return integrated_terms

def substitute_coefficients(integrated_terms: List[Tuple[Fraction, List[int]]], 
                          a_coeffs: List[str]) -> str:
    """
    ✅ 第二步：将积分后的a_j表达式代入v表达式
    
    参数:
        integrated_terms: 积分后项列表 [(系数, [a索引,...]), ...]
        a_coeffs: a系数的v表达式列表 [a0_expr, a1_expr, a2_expr, ...]
    
    返回:
        最终的复合表达式字符串
    """
    # 累积同类项（如所有a1^2项合并）
    result_dict = defaultdict(lambda: Fraction(0, 1))
    
    print("\n合并同类项:")
    print("-" * 50)
    
    for coeff, symbols in integrated_terms:
        # 生成唯一键用于合并
        if len(symbols) == 1:
            symbol_key = f"a{symbols[0]}"
        elif symbols[0] == symbols[1]:
            symbol_key = f"a{symbols[0]}^2"
        else:
            symbol_key = "*".join([f"a{s}" for s in symbols])
        
        # 累积系数
        result_dict[symbol_key] += coeff
    
    # 构建最终表达式
    terms = []
    for symbol_key, total_coeff in result_dict.items():
        if total_coeff == 0:
            continue
        
        print(f"  {symbol_key}: 系数 = {total_coeff}")
        
        # ✅ 正确映射：多项式中的a1对应a_coeffs[1]，a2对应a_coeffs[2]
        base_symbol = symbol_key.split('*')[0].split('^')[0]  # "a1"或"a2"
        a_index = int(base_symbol[1:])  # a1→1，a2→2
        
        # ✅ 安全检查
        if a_index >= len(a_coeffs):
            raise ValueError(f"多项式索引{a_index}超出系数范围[0, {len(a_coeffs)-1}]")
        
        a_expr = a_coeffs[a_index]
        coeff_str = format_coefficient(total_coeff)
        
        # 构建项字符串
        if '^2' in symbol_key:
            # 平方项需要括号: coeff*(expr)^2
            term = f"{coeff_str}*({a_expr})^2" if coeff_str else f"({a_expr})^2"
        else:
            # 一次项: coeff*expr
            term = f"{coeff_str}*{a_expr}" if coeff_str else f"{a_expr}"
        
        terms.append(term)
    
    return smart_join_terms(terms)

def smart_join_terms(terms: List[str]) -> str:
    """智能连接各项，处理符号"""
    if not terms:
        return "0"
    
    # 找到第一个非负项作为起点
    first_idx = 0
    while first_idx < len(terms) and terms[first_idx].startswith('-'):
        first_idx += 1
    
    if first_idx < len(terms):
        expr = terms[first_idx]
        # 负号项前置
        for t in terms[:first_idx]:
            expr = f"{t} + {expr}"
        # 正号项依次添加
        for t in terms[first_idx+1:]:
            if t.startswith('-'):
                expr += f" {t}"
            else:
                expr += f" + {t}"
    else:
        # 全为负号项
        expr = terms[0]
        for t in terms[1:]:
            expr += f" {t}"
    
    return expr

def evaluate_polynomial_integral_symbolic(polynomial: Dict[int, List[Tuple[float, List[int]]]], 
                                        a_coeffs: List[str]) -> str:
    """
    ✅ 主函数：先积分，后代入
    
    步骤:
        1. 对x在[-1/2, 1/2]上积分，消去x变量
        2. 将a_j替换为v表达式
    """
    integrated_terms = integrate_polynomial_x(polynomial)
    return substitute_coefficients(integrated_terms, a_coeffs)

# ============ 测试函数 ============

def generate_composite_expressions(k: int, polynomial: Dict[int, List[Tuple[float, List[int]]]], i: int = 0):
    """生成所有r值的复合表达式f_r"""
    f_r_dict = {}
    
    print(f"\n生成k={k}的复合表达式f_r")
    print("="*70)
    
    for r in range(k):
        print(f"\n--- r = {r} ---")
        
        # 1. 生成a系数的v表达式
        a_coeffs = solve_for_coefficients(k, r, i)
        print("\na系数的v表达式:")
        for idx, expr in enumerate(a_coeffs):
            print(f"  a_{idx} = {expr}")
        
        # 2. 生成复合表达式f_r
        f_r = evaluate_polynomial_integral_symbolic(polynomial, a_coeffs)
        f_r_dict[r] = f_r
        
        print(f"\nf_{r}（最终复合表达式）:")
        print("-" * 50)
        print_expression_pretty(f_r, indent="  ")
    
    # 总结
    print("\n" + "="*70)
    print("总结：所有r的复合表达式f_r")
    print("="*70)
    for r, expr in f_r_dict.items():
        print(f"\nf_{r} =")
        print_expression_pretty(expr, indent="  ")
    
    return f_r_dict

def print_expression_pretty(expr: str, indent: str = ""):
    """美化打印表达式"""
    if '+' not in expr:
        print(f"{indent}{expr}")
        return
    
    # 按顶层'+'分割
    lines = []
    bracket_depth = 0
    current = ""
    
    for char in expr:
        current += char
        if char == '(':
            bracket_depth += 1
        elif char == ')':
            bracket_depth -= 1
        elif char == '+' and bracket_depth == 0:
            lines.append(current[:-1].strip())
            current = ""
    
    if current:
        lines.append(current.strip())
    
    # 打印
    for i, line in enumerate(lines):
        if line.startswith('-'):
            print(f"{indent}{line}")
        elif i == 0:
            print(f"{indent}{line}")
        else:
            print(f"{indent}+ {line}")

def test_composite():
    """测试例子"""
    k = 3
    
    # 多项式：(a1^2 + 4*a2^2) + 4*a1*a2*x + 4*a2^2*x^2
    polynomial = {
        0: [(1.0, [1, 1]),   # a1^2
            (4.0, [2, 2])],  # a2^2（来自二阶导数部分）
        1: [(4.0, [1, 2])],  # a1*a2
        2: [(4.0, [2, 2])]   # a2^2（来自一阶导数平方的x^2项）
    }
    
    print("="*70)
    print("测试：多项式积分后复合a系数表达式（严格先积分后代入）")
    print("="*70)
    
    print("\n原始多项式（含x幂次）:")
    for exp, terms in polynomial.items():
        term_strs = []
        for coeff, symbols in terms:
            if len(symbols) == 1:
                term_strs.append(f"{coeff}*a{symbols[0]}")
            elif symbols[0] == symbols[1]:
                term_strs.append(f"{coeff}*a{symbols[0]}^2")
            else:
                term_strs.append(f"{coeff}*a{symbols[0]}*a{symbols[1]}")
        print(f"  x^{exp} 项: {' + '.join(term_strs)}")
    
    # 生成复合表达式
    f_r_dict = generate_composite_expressions(k, polynomial, i=0)

if __name__ == "__main__":
    test_composite()