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
    M = np.zeros((k, k), dtype=float)
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        for j in range(k):
            M[i, j] = compute_integral(alpha, beta, j)
    return M

def format_coefficient(frac: Fraction) -> str:
    """智能格式化Fraction为字符串"""
    if frac == 0:
        return "0"
    if frac.denominator == 1:
        return str(int(frac.numerator))
    return f"{frac.numerator}/{frac.denominator}"

def compute_lcm(numbers: List[int]) -> int:
    """计算最小公倍数"""
    if not numbers:
        return 1
    return reduce(lambda a, b: a * b // gcd(a, b), numbers)

def parse_and_optimize_expression(expr: str, exponent: int = 1) -> Tuple[str, Fraction]:
    """使用正则表达式安全解析表达式"""
    if expr == "0":
        return "0", Fraction(1, 1)
    
    if expr[0] not in '+-':
        expr = '+' + expr
    
    pattern = r'([+-])?(?:(\d+/\d+|\d+|[+-]?\d+/\d+|[+-]?\d+)\*)?(v\[i[+-]?\d*\]|v\[i\])'
    matches = re.findall(pattern, expr)
    
    coeff_vars = []
    denominators = []
    
    for match in matches:
        sign, coeff_part, var = match
        
        if coeff_part:
            coeff_str = coeff_part.replace('+', '').replace('-', '')
            if '/' in coeff_str:
                num, den = map(int, coeff_str.split('/'))
                coeff = Fraction(num, den)
            else:
                coeff = Fraction(int(coeff_str), 1)
        else:
            coeff = Fraction(1, 1)
        
        if sign == '-':
            coeff = -coeff
        
        coeff_vars.append((coeff, var))
        
        if coeff.denominator != 1:
            denominators.append(coeff.denominator)
    
    if not coeff_vars:
        return expr, Fraction(1, 1)
    
    lcm = compute_lcm(denominators) if denominators else 1
    
    int_coeffs = []
    for coeff, var in coeff_vars:
        int_coeff = coeff * lcm
        int_coeffs.append((int(int_coeff.numerator), var))
    
    neg_count = sum(1 for c, _ in int_coeffs if c < 0)
    pos_count = sum(1 for c, _ in int_coeffs if c > 0)
    
    factor = Fraction(1, lcm)
    if neg_count > pos_count:
        factor = -factor
        int_coeffs = [(-c, v) for c, v in int_coeffs]
    
    factor_with_exponent = factor ** exponent
    
    expr_parts = []
    for coeff, var in int_coeffs:
        if coeff == 1:
            expr_parts.append(f"+{var}")
        elif coeff == -1:
            expr_parts.append(f"-{var}")
        elif coeff > 0:
            expr_parts.append(f"+{coeff}*{var}")
        else:
            expr_parts.append(f"{coeff}*{var}")
    
    result_expr = ''.join(expr_parts)
    if result_expr.startswith('+'):
        result_expr = result_expr[1:]
    
    return result_expr, factor_with_exponent

def solve_for_coefficients_optimized(k: int, r: int, i: int = 0) -> List[Tuple[str, Fraction]]:
    """生成优化的a系数表达式"""
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
        
        optimized_expr, factor = parse_and_optimize_expression(expr, exponent=1)
        a_coeffs_optimized.append((optimized_expr, factor))
    
    return a_coeffs_optimized

# ============ 多项式积分与代入 ============

def integrate_polynomial_x(polynomial: Dict[int, List[Tuple[float, List[int]]]]) -> List[Tuple[Fraction, List[int]]]:
    """对x在[-1/2, 1/2]上积分"""
    a = Fraction(-1, 2)
    b = Fraction(1, 2)
    integrated_terms = []
    
    for exp, expr_list in polynomial.items():
        numerator = b**(exp + 1) - a**(exp + 1)
        integral_factor = Fraction(numerator, exp + 1)
        
        for coeff, symbols in expr_list:
            coeff_frac = Fraction(str(coeff))
            new_coeff = coeff_frac * integral_factor
            
            if new_coeff != 0:
                integrated_terms.append((new_coeff, symbols))
    
    return integrated_terms

def substitute_coefficients_optimized(integrated_terms: List[Tuple[Fraction, List[int]]], 
                                   a_coeffs_opt: List[Tuple[str, Fraction]]) -> Tuple[str, List]:
    """代入并排序"""
    result_dict = defaultdict(lambda: Fraction(0, 1))
    
    for coeff, symbols in integrated_terms:
        if len(symbols) == 1:
            symbol_key = f"a{symbols[0]}"
        elif symbols[0] == symbols[1]:
            symbol_key = f"a{symbols[0]}^2"
        else:
            symbol_key = "*".join([f"a{s}" for s in symbols])
        
        result_dict[symbol_key] += coeff
    
    # 构建带系数的项列表
    terms_with_coeffs = []
    
    for symbol_key, total_coeff in result_dict.items():
        if total_coeff == 0:
            continue
        
        base_symbol = symbol_key.split('*')[0].split('^')[0]
        a_index = int(base_symbol[1:])
        
        if a_index >= len(a_coeffs_opt):
            raise ValueError(f"多项式索引{a_index}超出范围")
        
        a_expr, a_factor = a_coeffs_opt[a_index]
        
        is_squared = '^2' in symbol_key
        exponent = 2 if is_squared else 1
        
        final_coeff = total_coeff * (a_factor ** exponent)
        coeff_abs = abs(final_coeff)
        
        coeff_str = format_coefficient(final_coeff)
        
        # 构建最终项
        if is_squared:
            term_str = f"{coeff_str}*({a_expr})^2" if coeff_str not in ["1", "-1"] else \
                       (f"({a_expr})^2" if coeff_str == "1" else f"-({a_expr})^2")
        else:
            term_str = f"{coeff_str}*{a_expr}" if coeff_str not in ["1", "-1"] else \
                       (f"{a_expr}" if coeff_str == "1" else f"-{a_expr}")
        
        terms_with_coeffs.append((coeff_abs, final_coeff, term_str))
    
    # ✅ 按系数绝对值降序排序
    terms_with_coeffs.sort(key=lambda x: x[0], reverse=True)
    
    # ✅ 提取排序后的项
    final_terms = [term_str for _, _, term_str in terms_with_coeffs]
    
    return smart_join_terms(final_terms), terms_with_coeffs

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
                                        k: int, r: int, i: int = 0) -> Tuple[str, List]:
    """完整流程：积分 → 代入"""
    a_coeffs_opt = solve_for_coefficients_optimized(k, r, i)
    integrated_terms = integrate_polynomial_x(polynomial)
    return substitute_coefficients_optimized(integrated_terms, a_coeffs_opt)

def generate_composite_expressions(k: int, polynomial: Dict[int, List[Tuple[float, List[int]]]], i: int = 0):
    """生成所有r值的复合表达式β_r（完整优化版）"""
    f_r_dict = {}
    all_terms_info = {}
    
    print(f"\n生成k={k}的复合表达式β_r")
    print("="*70)
    
    for r in range(k):
        final_expr, terms_with_coeffs = evaluate_polynomial_integral_symbolic(polynomial, k, r, i)
        f_r_dict[r] = final_expr
        all_terms_info[r] = terms_with_coeffs
        
        print(f"\nβ_{r} = {final_expr}")
    
    # ============ LaTeX格式总结 ============
    print("\n" + "="*70)
    print("LaTeX格式总结（按系数绝对值排序）")
    print("="*70)
    
    latex_dict = {}
    for r, terms_info in all_terms_info.items():
        latex_parts = []
        for coeff_abs, final_coeff, term_str in terms_info:
            # ✅ 修复：正确分离系数和主体，避免重复
            if '*(' in term_str:
                # 分离系数和主体
                coeff_part, main_part = term_str.split('*', 1)
                main_part = main_part.replace('v[i', 'v_{i').replace(']', '}')
                latex_term = main_part
            else:
                # 无显式系数
                latex_term = term_str.replace('v[i', 'v_{i').replace(']', '}')
            
            # 转换系数（确保只出现一次）
            coeff_str = format_coefficient(final_coeff)
            if coeff_str == "1":
                latex_parts.append(latex_term)
            elif coeff_str == "-1":
                latex_parts.append(f"-{latex_term}")
            else:
                latex_parts.append(f"{coeff_str}{latex_term}")
        
        latex_expr = " + ".join(latex_parts)
        latex_expr = latex_expr.replace("+ -", "- ")
        latex_dict[r] = latex_expr
        
        print(f"\n$\\beta_{r} = {latex_expr}$")
    
    # ============ 最终汇总（便于复制） ============
    print("\n" + "="*70)
    print("最终汇总（单行格式）")
    print("="*70)
    
    for r in range(k):
        latex_expr = latex_dict[r]
        print(f"β{r} = {latex_expr}")
    
    # ============ LaTeX代码块（便于复制） ============
    print("\n" + "="*70)
    print("LaTeX代码块")
    print("="*70)
    print("\n```latex")
    for r in range(k):
        latex_expr = latex_dict[r]
        print(f"\\beta_{r} = {latex_expr}")
    print("```")
    
    return f_r_dict

def print_expression_pretty(expr: str, indent: str = "", single_line: bool = False):
    """支持单行输出"""
    if single_line:
        print(f"{indent}{expr}")
        return
    
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
        0: [(1.0, [1, 1]), (4.0, [2, 2])],
        1: [(4.0, [1, 2])],
        2: [(4.0, [2, 2])]
    }
    
    print("="*70)
    print("测试：多项式积分后复合a系数表达式（完整优化版）")
    print("="*70)
    
    f_r_dict = generate_composite_expressions(k, polynomial, i=0)

if __name__ == "__main__":
    test_composite()