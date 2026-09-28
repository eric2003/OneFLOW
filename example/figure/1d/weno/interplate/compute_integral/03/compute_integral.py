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

def format_coefficient(frac: Fraction, latex_mode: bool = False) -> str:
    """
    智能格式化Fraction为字符串
    参数:
        latex_mode: 是否使用LaTeX格式（\\cfrac）
    """
    if frac == 0:
        return "0"
    if frac.denominator == 1:
        return str(int(frac.numerator))
    
    if latex_mode:
        return f"\\cfrac{{{frac.numerator}}}{{{frac.denominator}}}"
    
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

def integrate_polynomial_x(polynomial: Dict[int, List[Tuple[float, List[int]]]]) -> List[Tuple[Fraction, List[int]]]:
    """对x在[-1/2, 1/2]上积分"""
    a = Fraction(-1, 2)
    b = Fraction(1, 2)
    integrated_terms = []
    
    for exp, expr_list in polynomial.items():
        numerator = b**(exp + 1) - a**(exp + 1)
        integral_factor = Fraction(numerator, exp + 1)
        print(f'exp={exp},expr_list={expr_list}')
        
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

# ============ 核心修改：两项表达式优化 ============

def optimize_two_term_expression(expr: str, is_latex: bool = False) -> str:
    """
    ✅ 核心函数：优化两项表达式
    
    规则：
    1. 仅在两项时生效
    2. -v[i-1]+v[i+1] -> v[i-1]-v[i+1]（首项为正，保持顺序）
    3. 下标顺序：i-1在i+1前面
    
    参数:
        expr: 括号内表达式，如 "-v[i-1]+v[i+1]"
        is_latex: 是否为LaTeX格式
    
    返回:
        优化后表达式
    """
    # 检查是否为两项
    var_pattern = r'v_\{i[+-]\d*\}' if is_latex else r'v\[i[+-]?\d*\]'
    if len(re.findall(var_pattern, expr)) != 2:
        return expr
    
    # 匹配项：符号、系数、下标
    if is_latex:
        pattern = r'([+-]?)(\d*)?v_\{i([+-]\d*)\}'
    else:
        pattern = r'([+-]?)(\d*\*)?v\[i([+-]\d*)\]'
    
    matches = re.findall(pattern, expr)
    
    if len(matches) != 2:
        return expr
    
    # 解析两项
    terms = []
    for sign, coeff_part, index in matches:
        if not sign:
            sign = '+'
        
        if is_latex:
            coeff = coeff_part if coeff_part else '1'
        else:
            coeff = coeff_part.rstrip('*') if coeff_part else '1'
        
        terms.append((sign, coeff, index))
    
    # 检查首项是否为负
    if terms[0][0] == '-':
        sign1, coeff1, idx1 = terms[0]  # -v[i-1]
        sign2, coeff2, idx2 = terms[1]  # +v[i+1]
        
        # ✅ 保持顺序：i-1在i+1前面
        # -v[i-1] + v[i+1] -> v[i-1] - v[i+1]
        new_parts = []
        
        # 首项（原第一项，变号）
        if is_latex:
            if coeff1 == '1':
                new_parts.append(f"v_{{i{idx1}}}")
            else:
                new_parts.append(f"{coeff1}v_{{i{idx1}}}")
        else:
            if coeff1 == '1':
                new_parts.append(f"v[i{idx1}]")
            else:
                new_parts.append(f"{coeff1}*v[i{idx1}]")
        
        # 次项（原第二项，符号取反）
        if is_latex:
            if coeff2 == '1':
                new_parts.append(f"-v_{{i{idx2}}}")
            else:
                new_parts.append(f"-{coeff2}v_{{i{idx2}}}")
        else:
            if coeff2 == '1':
                new_parts.append(f"-v[i{idx2}]")
            else:
                new_parts.append(f"-{coeff2}*v[i{idx2}]")
        
        return ''.join(new_parts)
    
    return expr

def apply_two_term_optimization_to_python(term_str: str) -> str:
    """
    ✅ 对Python表达式应用两项优化
    
    示例：
    "1/4*(-v[i-1]+v[i+1])^2" -> "1/4*(v[i-1]-v[i+1])^2"
    """
    if '*(' not in term_str or '^2' not in term_str:
        return term_str
    
    try:
        coeff_part, main_part = term_str.split('*', 1)
        main_part = main_part.strip()
        
        match = re.match(r'^\((.+)\)\^2$', main_part)
        if not match:
            return term_str
        
        inner_expr = match.group(1)
        
        # ✅ 应用两项优化（保持顺序）
        optimized_inner = optimize_two_term_expression(inner_expr, is_latex=False)
        
        if optimized_inner != inner_expr:
            term_str = f"{coeff_part}*({optimized_inner})^2"
    
    except:
        pass
    
    return term_str

def apply_two_term_optimization_to_latex(term_str: str, latex_coeff: str) -> str:
    """
    ✅ 对LaTeX表达式应用两项优化
    
    示例：
    "\\cfrac{1}{4}(-v_{i-1}+v_{i+1})^2" -> "\\cfrac{1}{4}(v_{i-1}-v_{i+1})^2"
    """
    # 提取括号部分
    match = re.search(r'\((.+)\)\^2', term_str)
    if not match:
        return term_str
    
    inner_expr = match.group(1)
    
    # ✅ 应用两项优化（保持顺序）
    optimized_inner = optimize_two_term_expression(inner_expr, is_latex=True)
    
    if optimized_inner != inner_expr:
        term_str = term_str.replace(f"({inner_expr})^2", f"({optimized_inner})^2")
    
    return term_str

def evaluate_polynomial_integral_symbolic(polynomial: Dict[int, List[Tuple[float, List[int]]]], 
                                        k: int, r: int, i: int = 0) -> Tuple[str, List]:
    """完整流程：积分 → 代入"""
    a_coeffs_opt = solve_for_coefficients_optimized(k, r, i)
    print(f'a_coeffs_opt={a_coeffs_opt}')
    integrated_terms = integrate_polynomial_x(polynomial)
    print(f'integrated_terms={integrated_terms}')
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
            latex_coeff = format_coefficient(final_coeff, latex_mode=True)
            
            # ✅ 转换为LaTeX并应用两项优化
            latex_main = apply_two_term_optimization_to_latex(term_str, latex_coeff)
            
            # ✅ 构建完整LaTeX项
            if latex_coeff == "1":
                latex_term = latex_main
            elif latex_coeff == "-1":
                latex_term = f"-{latex_main}"
            else:
                latex_term = f"{latex_coeff}{latex_main}"
            
            latex_parts.append(latex_term)
        
        latex_expr = " + ".join(latex_parts)
        latex_expr = latex_expr.replace("+ -", "- ")
        latex_dict[r] = latex_expr
        
        print(f"\n$\\beta_{r} = {latex_expr}$")
    
    # ============ 最终汇总（单行格式 - 应用两项优化） ============
    print("\n" + "="*70)
    print("最终汇总（单行格式）")
    print("="*70)
    
    python_dict = {}
    for r, terms_info in all_terms_info.items():
        python_parts = []
        for coeff_abs, final_coeff, term_str in terms_info:
            # ✅ 应用两项优化到Python表达式
            optimized_term = apply_two_term_optimization_to_python(term_str)
            python_parts.append(optimized_term)
        
        python_expr = " + ".join(python_parts).replace("+ -", "- ")
        python_dict[r] = python_expr
        
        print(f"β{r} = {python_expr}")
    
    # ============ LaTeX代码块 ============
    print("\n" + "="*70)
    print("LaTeX代码块（统一在array环境中）")
    print("="*70)
    print("\n```latex")
    print("\\begin{array}{l}")
    for r in range(k):
        if r == k - 1:
            print(f"    \\beta_{r} = {latex_dict[r]}")
        else:
            print(f"    \\beta_{r} = {latex_dict[r]}\\\\")
    print("\\end{array}")
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