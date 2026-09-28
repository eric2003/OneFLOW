import numpy as np
from fractions import Fraction
from collections import defaultdict
from typing import List, Tuple, Dict
from math import gcd
from functools import reduce
import re

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

# ============ 核心改进：优化整数化与符号处理 ============

def format_coefficient(frac: Fraction) -> str:
    """智能格式化Fraction为字符串"""
    if frac == 0:
        return "0"
    if frac.denominator == 1:
        value = int(frac.numerator)
        return str(value)
    return f"{frac.numerator}/{frac.denominator}"

def compute_lcm(numbers: List[int]) -> int:
    """计算最小公倍数"""
    if not numbers:
        return 1
    return reduce(lambda a, b: a * b // gcd(a, b), numbers)

def parse_and_optimize_expression(expr: str, exponent: int = 1) -> Tuple[str, Fraction]:
    """
    ✅ 修复版：使用正则表达式安全解析表达式
    
    参数:
        expr: 原始表达式，如 "-3/2*v[i] + 2*v[i+1] - 1/2*v[i+2]"
        exponent: 指数（1为线性，2为平方）
    
    返回:
        (优化后表达式, 提取的因子)
    """
    if expr == "0":
        return "0", Fraction(1, 1)
    
    # 1. 如果表达式以正负号开头，添加隐含的+号
    if expr[0] not in '+-':
        expr = '+' + expr
    
    # 2. 使用正则表达式安全解析表达式
    # 模式: ([+-]) 可选符号 + (系数) + (变量名)
    # 系数支持: 整数、分数、带符号
    # 变量名: v[i] 或 v[i+1] 或 v[i-1] 等
    pattern = r'([+-])?(?:(\d+/\d+|\d+|[+-]?\d+/\d+|[+-]?\d+)\*)?(v\[i[+-]?\d*\]|v\[i\])'
    
    matches = re.findall(pattern, expr)
    
    coeff_vars = []  # List of (Fraction系数, 变量字符串)
    denominators = []  # 所有分母
    
    for match in matches:
        sign, coeff_part, var = match
        
        # 处理系数
        if coeff_part:
            # 清理系数字符串
            coeff_str = coeff_part.replace('+', '').replace('-', '')
            if '/' in coeff_str:
                num, den = map(int, coeff_str.split('/'))
                coeff = Fraction(num, den)
            else:
                coeff = Fraction(int(coeff_str), 1)
        else:
            # 没有显式系数，如 "+v[i]"
            coeff = Fraction(1, 1)
        
        # 应用符号
        if sign == '-':
            coeff = -coeff
        
        coeff_vars.append((coeff, var))
        
        # 记录分母
        if coeff.denominator != 1:
            denominators.append(coeff.denominator)
    
    # 检查是否有任何项被解析
    if not coeff_vars:
        print(f"警告：表达式 '{expr}' 未解析出任何项")
        return expr, Fraction(1, 1)
    
    # 3. 计算最小公倍数
    lcm = compute_lcm(denominators) if denominators else 1
    
    # 4. 转换为整数系数
    int_coeffs = []
    for coeff, var in coeff_vars:
        int_coeff = coeff * lcm
        int_coeffs.append((int(int_coeff.numerator), var))
    
    # 5. 符号优化：统计正负号数量
    neg_count = sum(1 for c, _ in int_coeffs if c < 0)
    pos_count = sum(1 for c, _ in int_coeffs if c > 0)
    
    # 如果负数多于正数，提取负号
    factor = Fraction(1, lcm)
    if neg_count > pos_count:
        factor = -factor
        int_coeffs = [(-c, v) for c, v in int_coeffs]
    
    # 6. 根据指数调整因子
    factor_with_exponent = factor ** exponent
    
    # 7. 重建表达式字符串
    expr_parts = []
    for coeff, var in int_coeffs:
        if coeff == 1:
            expr_parts.append(f"+{var}")
        elif coeff == -1:
            expr_parts.append(f"-{var}")
        elif coeff > 0:
            expr_parts.append(f"+{coeff}*{var}")
        else:  # coeff < 0
            expr_parts.append(f"{coeff}*{var}")
    
    # 拼接并清理开头
    result_expr = ''.join(expr_parts)
    if result_expr.startswith('+'):
        result_expr = result_expr[1:]
    
    return result_expr, factor_with_exponent

def solve_for_coefficients_optimized(k: int, r: int, i: int = 0) -> List[Tuple[str, Fraction]]:
    """
    ✅ 生成优化的a系数表达式（整数化+符号优化）
    
    返回: [(优化后表达式, 提取因子), ...]
    """
    M = compute_matrix_M(k, r)
    try:
        M_inv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        raise ValueError(f"矩阵M (k={k}, r={r})不可逆！")
    
    a_coeffs_optimized = []
    for j in range(k):
        nonzero_items = []
        for m in range(k):
            coeff = M_inv[j, m]
            if abs(coeff) > 1e-10:
                offset = m - r
                v_str = f"v[i{offset:+d}]" if offset != 0 else "v[i]"
                nonzero_items.append((coeff, v_str))
        
        if not nonzero_items:
            a_coeffs_optimized.append(("0", Fraction(1, 1)))
            continue
        
        # 构建原始表达式字符串
        expr_parts = []
        for coeff, v_str in nonzero_items:
            frac = Fraction(coeff).limit_denominator(1000)
            if coeff == 1.0:
                expr_parts.append(f"+{v_str}")
            elif coeff == -1.0:
                expr_parts.append(f"-{v_str}")
            else:
                expr_parts.append(f"{frac}*{v_str}")
        
        expr = ''.join(expr_parts)
        if expr.startswith('+'):
            expr = expr[1:]
        
        # 优化表达式（exponent=1，因为是a_j本身）
        optimized_expr, factor = parse_and_optimize_expression(expr, exponent=1)
        
        a_coeffs_optimized.append((optimized_expr, factor))
    
    return a_coeffs_optimized

# ============ 多项式积分与代入 ============

def integrate_polynomial_x(polynomial: Dict[int, List[Tuple[float, List[int]]]]) -> List[Tuple[Fraction, List[int]]]:
    """
    ✅ 第一步：对x在[-1/2, 1/2]上积分
    """
    a = Fraction(-1, 2)
    b = Fraction(1, 2)
    
    integrated_terms = []
    print("\n积分过程详细计算:")
    print("-" * 50)
    
    for exp, expr_list in polynomial.items():
        numerator = b**(exp + 1) - a**(exp + 1)
        integral_factor = Fraction(numerator, exp + 1)
        
        print(f"  x^{exp} 在[-1/2,1/2]上的积分: {integral_factor}")
        
        for coeff, symbols in expr_list:
            coeff_frac = Fraction(str(coeff))
            new_coeff = coeff_frac * integral_factor
            
            # 生成符号表示
            if len(symbols) == 1:
                symbol_str = f"a{symbols[0]}"
            elif symbols[0] == symbols[1]:
                symbol_str = f"a{symbols[0]}^2"
            else:
                symbol_str = "*".join([f"a{s}" for s in symbols])
            
            if new_coeff != 0:
                print(f"    {coeff_frac} * {symbol_str} × {integral_factor} = {new_coeff} * {symbol_str}")
                integrated_terms.append((new_coeff, symbols))
    
    return integrated_terms

def substitute_coefficients_optimized(integrated_terms: List[Tuple[Fraction, List[int]]], 
                                   a_coeffs_opt: List[Tuple[str, Fraction]]) -> str:
    """
    ✅ 第二步：将积分后的a_j表达式代入v表达式（支持优化格式）
    """
    result_dict = defaultdict(lambda: Fraction(0, 1))
    
    print("\n合并同类项:")
    print("-" * 50)
    
    for coeff, symbols in integrated_terms:
        if len(symbols) == 1:
            symbol_key = f"a{symbols[0]}"
        elif symbols[0] == symbols[1]:
            symbol_key = f"a{symbols[0]}^2"
        else:
            symbol_key = "*".join([f"a{s}" for s in symbols])
        
        result_dict[symbol_key] += coeff
    
    # 构建最终表达式
    final_terms = []
    for symbol_key, total_coeff in result_dict.items():
        if total_coeff == 0:
            continue
        
        print(f"  {symbol_key}: 总系数 = {total_coeff}")
        
        base_symbol = symbol_key.split('*')[0].split('^')[0]
        a_index = int(base_symbol[1:])
        
        if a_index >= len(a_coeffs_opt):
            raise ValueError(f"多项式索引{a_index}超出范围[0, {len(a_coeffs_opt)-1}]")
        
        a_expr, a_factor = a_coeffs_opt[a_index]
        
        is_squared = '^2' in symbol_key
        exponent = 2 if is_squared else 1
        
        # ✅ 计算最终系数
        final_coeff = total_coeff * (a_factor ** exponent)
        coeff_str = format_coefficient(final_coeff)
        
        # ✅ 构建最终项
        if is_squared:
            if coeff_str == "1":
                term = f"({a_expr})^2"
            elif coeff_str == "-1":
                term = f"-({a_expr})^2"
            else:
                term = f"{coeff_str}*({a_expr})^2"
        else:
            if coeff_str == "1":
                term = f"{a_expr}"
            elif coeff_str == "-1":
                term = f"-{a_expr}"
            else:
                term = f"{coeff_str}*{a_expr}"
        
        final_terms.append(term)
    
    return smart_join_terms(final_terms)

def smart_join_terms(terms: List[str]) -> str:
    """智能连接各项"""
    if not terms:
        return "0"
    
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

def evaluate_polynomial_integral_symbolic(polynomial: Dict[int, List[Tuple[float, List[int]]]], 
                                        k: int, r: int, i: int = 0) -> str:
    """
    ✅ 完整流程：积分 → 代入（支持优化）
    """
    a_coeffs_opt = solve_for_coefficients_optimized(k, r, i)
    
    print(f"\nr = {r} 的优化a系数:")
    for idx, (expr, factor) in enumerate(a_coeffs_opt):
        if factor == 1:
            print(f"  a_{idx} = ({expr})")
        else:
            print(f"  a_{idx} = {format_coefficient(factor)} * ({expr})")
    
    integrated_terms = integrate_polynomial_x(polynomial)
    return substitute_coefficients_optimized(integrated_terms, a_coeffs_opt)

def generate_composite_expressions(k: int, polynomial: Dict[int, List[Tuple[float, List[int]]]], i: int = 0):
    """生成所有r值的复合表达式f_r（完整优化版）"""
    f_r_dict = {}
    
    print(f"\n生成k={k}的复合表达式f_r（整数化+符号优化）")
    print("="*70)
    
    for r in range(k):
        print(f"\n--- r = {r} ---")
        
        f_r = evaluate_polynomial_integral_symbolic(polynomial, k, r, i)
        f_r_dict[r] = f_r
        
        print(f"\nf_{r}（最终优化表达式）:")
        print("-" * 50)
        print_expression_pretty(f_r, indent="  ")
    
    # 总结
    print("\n" + "="*70)
    print("总结：所有r的优化表达式f_r")
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
    
    polynomial = {
        0: [(1.0, [1, 1]),   # a1^2
            (4.0, [2, 2])],  # a2^2
        1: [(4.0, [1, 2])],  # a1*a2
        2: [(4.0, [2, 2])]   # a2^2
    }
    
    print("="*70)
    print("测试：多项式积分后复合a系数表达式（完整优化版）")
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
    
    f_r_dict = generate_composite_expressions(k, polynomial, i=0)

if __name__ == "__main__":
    test_composite()