from fractions import Fraction
from collections import Counter
from collections import defaultdict
from typing import List, Tuple, Dict
import numpy as np
import math
from math import gcd
from functools import reduce

# 类型别名
Term = Tuple[int, List[int]]  # (系数, 符号下标列表)
Expression = List[Term]        # Term 的列表
Polynomial = Dict[int, Expression]  # {指数: 表达式}

def extract_max_common_factorBAK(numbers):
    """
    从数值列表中提取最大的公共因子，使得剩余部分为互质整数列表
    
    参数:
        numbers: 数值列表（可包含整数、浮点数、字符串分数等）
        
    返回:
        tuple: (factor, simplified_list)
               factor: Fraction类型，提取的公共因子
               simplified_list: 整数列表，最简形式（gcd=1）
    """
    # 1. 将所有输入转换为Fraction，确保精确的有理数运算
    #    Fraction(0.5) = 1/2, Fraction(-1) = -1/1, Fraction("1/3") = 1/3
    print(f'numbers={numbers}')
    fractions = [Fraction(x) for x in numbers]
    
    # 2. 特殊情况处理
    if not fractions:  # 空列表
        return Fraction(1, 1), []
    
    if all(f == 0 for f in fractions):  # 全零列表
        return Fraction(1, 1), [0] * len(fractions)
    
    # 3. 提取所有分子和分母
    numerators = [f.numerator for f in fractions]
    denominators = [f.denominator for f in fractions]
    
    # 4. 计算分子的最大公约数(GCD)
    #    reduce函数连续应用gcd: gcd(gcd(p1,p2), p3), ...
    numerator_gcd = reduce(gcd, numerators)
    
    # 5. 计算分母的最小公倍数(LCM)
    def lcm(a, b):
        """计算两个数的最小公倍数：lcm(a,b) = |a×b| / gcd(a,b)"""
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // gcd(a, b)
    
    #    连续应用lcm: lcm(lcm(q1,q2), q3), ...
    denominator_lcm = reduce(lcm, denominators)
    
    # 6. 最大公共因子 = 分子GCD / 分母LCM
    common_factor = Fraction(numerator_gcd, denominator_lcm)
    
    # 7. 计算简化后的列表
    simplified_fractions = [f / common_factor for f in fractions]
    
    # 8. 验证并转换为整数列表
    #    理论上所有分母都应为1，因为除以了最大公共因子
    simplified_integers = [sf.numerator for sf in simplified_fractions]
    
    # 9. 额外验证：确保结果整数列表互质（gcd=1）
    verification_gcd = reduce(gcd, simplified_integers)
    
    #    如果verification_gcd≠1，说明还能继续提取，调整结果
    if verification_gcd != 1:
        true_factor = common_factor * verification_gcd
        simplified_integers = [x // verification_gcd for x in simplified_integers]
        return true_factor, simplified_integers
    
    return common_factor, simplified_integers
    
def extract_max_common_factorBAK1(numbers):
    """
    从数值列表中提取最大的公共因子，使得剩余部分为互质整数列表
    
    参数:
        numbers: 数值列表（支持int、float、numpy标量、字符串等）
        
    返回:
        tuple: (factor, simplified_list)
               factor: Fraction类型，提取的公共因子
               simplified_list: 整数列表，最简形式（gcd=1）
    """
    
    def _to_python_number(x):
        """
        关键修复：将numpy标量转换为Python原生类型
        """
        if isinstance(x, (np.integer, np.floating, np.ndarray)):
            # numpy标量转Python标量
            return x.item()
        return x
        
        
    print(f'numbers={numbers}')
    
    # 1. 转换前先处理numpy类型
    processed_numbers = [_to_python_number(x) for x in numbers]
    
    # 2. 转换为Fraction（现在安全了）
    fractions = [Fraction(x) for x in processed_numbers]
    print(f'fractions={fractions}')
    
    # 3. 后续逻辑不变
    if not fractions:
        return Fraction(1, 1), []
    
    if all(f == 0 for f in fractions):
        return Fraction(1, 1), [0] * len(fractions)
    
    # 提取分子和分母
    numerators = [f.numerator for f in fractions]
    denominators = [f.denominator for f in fractions]
    
    # 计算分子GCD
    numerator_gcd = reduce(gcd, numerators)
    
    # 计算分母LCM
    def lcm(a, b):
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // gcd(a, b)
    
    denominator_lcm = reduce(lcm, denominators)
    
    # 最大公共因子 = 分子GCD / 分母LCM
    common_factor = Fraction(numerator_gcd, denominator_lcm)
    
    # 简化列表
    simplified_fractions = [f / common_factor for f in fractions]
    simplified_integers = [sf.numerator for sf in simplified_fractions]
    
    # 验证互质
    verification_gcd = reduce(gcd, simplified_integers)
    if verification_gcd != 1:
        true_factor = common_factor * verification_gcd
        simplified_integers = [x // verification_gcd for x in simplified_integers]
        return true_factor, simplified_integers
    
    return common_factor, simplified_integers    
    
def extract_max_common_factor(numbers, max_denominator=1000000, tolerance=1e-12):
    """
    提取最大公共因子，支持numpy浮点数精度修复
    
    参数:
        numbers: 数值列表（支持numpy标量、float、int等）
        max_denominator: 限制分母的最大值，用于修复浮点误差
        tolerance: 接近零值的容差阈值
    """
    
    def _to_python_number(x):
        """转换numpy标量为Python原生类型"""
        if isinstance(x, (np.integer, np.floating)):
            return x.item()
        if isinstance(x, np.ndarray) and x.size == 1:
            return x.item()
        return x
    
    def _smart_fraction(x):
        """智能转换为Fraction，自动修复精度问题"""
        # 1. 转换为Python原生数值
        val = _to_python_number(x)
        
        # 2. 处理接近零的值
        if abs(val) < tolerance:
            return Fraction(0, 1)
        
        # 3. 对浮点数先限制分母再转换
        if isinstance(val, float):
            # 先创建Fraction，再限制分母复杂度
            return Fraction(val).limit_denominator(max_denominator)
        
        # 4. 其他类型直接转换
        return Fraction(val)
    
    # 主逻辑
    fractions = [_smart_fraction(x) for x in numbers]
    
    if not fractions:
        return Fraction(1, 1), []
    
    if all(f == 0 for f in fractions):
        return Fraction(1, 1), [0] * len(fractions)
    
    # 计算分子GCD和分母LCM
    numerators = [f.numerator for f in fractions]
    denominators = [f.denominator for f in fractions]
    
    numerator_gcd = reduce(gcd, numerators)
    
    def lcm(a, b):
        if a == 0 or b == 0:
            return 0
        return abs(a * b) // gcd(a, b)
    
    denominator_lcm = reduce(lcm, denominators)
    
    # 最大公共因子 = 分子GCD / 分母LCM
    common_factor = Fraction(numerator_gcd, denominator_lcm)
    
    # 简化并验证互质
    simplified = [f / common_factor for f in fractions]
    simplified_integers = [sf.numerator for sf in simplified]
    
    # 确保gcd=1
    final_gcd = reduce(gcd, simplified_integers)
    if final_gcd != 1:
        common_factor *= final_gcd
        simplified_integers = [x // final_gcd for x in simplified_integers]
    
    return common_factor, simplified_integers    

def term_multiply(t1: Term, t2: Term) -> Term:
    """
    两个 Term 相乘
    示例: (2, [1]) * (3, [2]) = (6, [1, 2])
    """
    coeff1, symbols1 = t1
    coeff2, symbols2 = t2
    # 合并符号列表并排序
    new_symbols = sorted(symbols1 + symbols2)
    return (coeff1 * coeff2, new_symbols)

def expression_add(expr1: Expression, expr2: Expression) -> Expression:
    """
    两个 Expression 相加，合并同类项
    示例: [(2,[1])] + [(3,[2]), (2,[1])] = [(4,[1]), (3,[2])]
    """
    # 用字典合并：键是符号元组，值是系数和
    term_dict = defaultdict(int)
    
    for coeff, symbols in expr1 + expr2:
        key = tuple(symbols)
        term_dict[key] += coeff
    
    # 过滤系数为0的项，转换回列表
    result = [(coeff, list(symbols)) for symbols, coeff in term_dict.items() if coeff != 0]
    return result

def polynomial_square(polynomial: Polynomial) -> Polynomial:
    """
    多项式平方展开
    算法: 遍历所有指数对 (i, j)，对应项相乘后指数相加
    """
    result: Polynomial = defaultdict(list)
    exps = sorted(polynomial.keys())
    
    # 1. 平方项 (i == j)
    for exp in exps:
        expr = polynomial[exp]
        new_exp = exp * 2
        
        # 表达式与自身相乘
        for i, term1 in enumerate(expr):
            # 平方项
            squared_term = term_multiply(term1, term1)
            result[new_exp].append(squared_term)
            
            # 交叉项 (2ab)
            for term2 in expr[i+1:]:
                cross_term = term_multiply(term1, term2)
                # 系数乘以2
                double_cross = (2 * cross_term[0], cross_term[1])
                result[new_exp].append(double_cross)
    
    # 2. 交叉项 (i != j)
    for i in range(len(exps)):
        for j in range(i + 1, len(exps)):
            exp_i, exp_j = exps[i], exps[j]
            new_exp = exp_i + exp_j
            
            for term1 in polynomial[exp_i]:
                for term2 in polynomial[exp_j]:
                    product = term_multiply(term1, term2)
                    # 乘以2
                    result[new_exp].append((2 * product[0], product[1]))
    
    # 合并同类项
    final_result = {}
    for exp, expr in result.items():
        final_result[exp] = expression_add(expr, [])
    
    return final_result

def integrate_polynomial(polynomial: Polynomial) -> Polynomial:
    """
    对多项式进行积分: ∫ x^k dx = x^(k+1)/(k+1)
    符号部分保持不变，数值系数除以 (k+1)
    """
    integrated: Polynomial = {}
    
    for exp, expr in polynomial.items():
        new_exp = exp + 1
        integrated[new_exp] = []
        
        for coeff, symbols in expr:
            # ∫ c*x^exp dx = (c/(exp+1))*x^(exp+1)
            # 系数除以 (exp+1)，可能是分数
            new_coeff = coeff / (exp + 1)
            integrated[new_exp].append((new_coeff, symbols))
    
    return integrated

def format_term(term: Term) -> str:
    """格式化单个 Term"""
    coeff, symbols = term
    
    if not symbols:  # 常数项
        return str(coeff)
    
    # 统计符号出现次数，处理 a1*a1 → a1^2
    symbol_counts = defaultdict(int)
    for idx in symbols:
        symbol_counts[idx] += 1
    
    parts = []
    for idx in sorted(symbol_counts.keys()):
        count = symbol_counts[idx]
        if count == 1:
            parts.append(f"a{idx}")
        else:
            parts.append(f"a{idx}^{count}")
    
    symbol_str = "*".join(parts)
    
    # 处理系数为1的情况（省略系数）
    if coeff == 1:
        return symbol_str
    # 处理分数系数
    elif isinstance(coeff, float) and not coeff.is_integer():
        return f"({coeff})*{symbol_str}"
    else:
        return f"{int(coeff)}*{symbol_str}"

def sum_integrals_same_bounds(polynomials: List[Polynomial], 
                              a: float = -0.5, 
                              b: float = 0.5) -> Expression:
    """
    多个多项式在相同区间[a, b]上积分后求和
    
    参数:
        polynomials: 多项式列表 [poly1, poly2, ...]
        a, b: 积分下限和上限（所有多项式相同）
    
    返回:
        合并后的符号表达式（不含x）
    """
    # 初始化结果为空表达式
    total_expression = []  # 相当于0
    
    # 遍历每个多项式
    for idx, poly in enumerate(polynomials):
        # 对当前多项式在[a, b]上积分
        integral_result = evaluate_polynomial_integral(poly, a, b)
        
        # 打印每个多项式的积分结果（调试用）
        print(f"  多项式{idx+1}积分: {format_expression(integral_result)}")
        
        # 累加到总表达式中
        total_expression = expression_add(total_expression, integral_result)
    
    return total_expression
    
def format_expression(expr: Expression) -> str:
    """格式化纯符号表达式"""
    if not expr:
        return "0"
    
    term_strs = []
    for coeff, symbols in expr:
        if len(symbols) == 1:
            term_strs.append(f"{coeff}*a{symbols[0]}")
        elif symbols[0] == symbols[1]:
            term_strs.append(f"{coeff}*a{symbols[0]}^2")
        else:
            symbol_str = "*".join([f"a{s}" for s in symbols])
            term_strs.append(f"{coeff}*{symbol_str}")
    
    return " + ".join(term_strs)    

def print_polynomial(polynomial: Polynomial, title: str = ""):
    """打印多项式，带标题"""
    if title:
        print(f"\n{title}:")
    
    if not polynomial:
        print("0")
        return
    
    sorted_exps = sorted(polynomial.keys())
    
    for idx, exp in enumerate(sorted_exps):
        expr_str = format_expression(polynomial[exp])
        
        # 处理常数项 (x^0)
        if exp == 0:
            print(f"({expr_str})", end='')
        else:
            print(f"({expr_str})*x^{exp}", end='')
        
        # 打印连接符
        if idx < len(sorted_exps) - 1:
            print(" + ", end='')
        else:
            print()
            
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
    
def evaluate_polynomial_integral(polynomial: Polynomial, 
                                 a: float = -0.5, 
                                 b: float = 0.5) -> Expression:
    """
    在区间 [a, b] 上对多项式进行定积分
    返回：纯符号表达式（不再含x）
    
    示例: 
        ∫[-0.5,0.5] [a1^2 + 4*a1*a2*x + 4*a2^2*x^2] dx
        = 1*a1^2 + 0*a1*a2 + (1/3)*a2^2
    """
    # 用字典存储符号表达式：键=符号元组，值=总系数
    # 例如: {(1,1): 1.0, (2,2): 0.3333}
    result_dict = defaultdict(float)
    
    # 遍历多项式的每一项：指数 -> 表达式
    # polynomial = {0: [(1, [1,1])], 1: [(4, [1,2])], 2: [(4, [2,2])]}
    for exp, expr in polynomial.items():
        # exp: x的指数（0, 1, 2...）
        # expr: 对应指数的表达式列表
        
        # 计算 ∫ x^exp dx = [x^(exp+1)/(exp+1)] 在 a 到 b 的值
        # integral_factor = (b^(exp+1) - a^(exp+1)) / (exp+1)
        integral_factor = (b ** (exp + 1) - a ** (exp + 1)) / (exp + 1)
        # 示例: exp=0 → (0.5^1 - (-0.5)^1)/1 = (0.5 + 0.5) = 1
        #       exp=1 → (0.5^2 - (-0.5)^2)/2 = (0.25 - 0.25)/2 = 0
        #       exp=2 → (0.5^3 - (-0.5)^3)/3 = (0.125 + 0.125)/3 = 0.25/3 ≈ 0.08333
        
        # 遍历该指数下的所有项
        for coeff, symbols in expr:
            # coeff: 数值系数（如 4）
            # symbols: 符号下标列表（如 [1,2] 表示 a1*a2）
            
            # 生成符号键：排序后的符号元组（确保 a1*a2 和 a2*a1 相同）
            # tuple(sorted([1,2])) = (1,2)
            symbol_key = tuple(sorted(symbols))
            
            # 累加积分后的系数
            # 总系数 = 原系数 * 积分因子
            # 例: exp=2, coeff=4, integral_factor≈0.08333
            #     contribution = 4 * 0.08333 ≈ 0.3333
            contribution = coeff * integral_factor
            
            result_dict[symbol_key] += contribution
            # 如果是相同符号项，系数会累加（如多处出现 a1^2）
    
    # 转换回 Expression 格式：[(系数, [符号列表]), ...]
    # 过滤掉系数为0的项（如奇函数项）
    result_expression = [
        (coeff, list(symbols)) 
        for symbols, coeff in result_dict.items() 
        if abs(coeff) > 1e-10  # 避免浮点误差
    ]
    
    return result_expression
    
def evaluate_and_print(polynomial: Polynomial, title: str = ""):
    """求值并打印结果"""
    if title:
        print(f"\n{title}")
    
    # 在 [-0.5, 0.5] 上积分
    result = evaluate_polynomial_integral(polynomial, -0.5, 0.5)
    
    # 格式化输出
    result_str = format_expression(result)
    print(f"∫[-1/2,1/2] P(x) dx = {result_str}")
    
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

def print_polynomial_old_style(polynomial, title=""):
    """
    改造重点1：创建兼容旧格式的打印函数
    总是显示 *x^e，包括 *x^0
    """
    if title:
        print(f"\n{title}")
    
    if not polynomial:
        print("0")
        return
    
    # 按指数排序
    sorted_exps = sorted(polynomial.keys())
    
    for idx, exp in enumerate(sorted_exps):
        # 格式化x^e项的系数部分
        expr = polynomial[exp]
        term_strs = []
        for coef, symbols in expr:
            # 如果符号列表只有一个元素
            if len(symbols) == 1:
                term_strs.append(f"{coef}*a{symbols[0]}")
            else:
                # 多个符号相乘（虽然这里用不到，但为完整性保留）
                symbol_str = "*".join([f"a{s}" for s in symbols])
                term_strs.append(f"{coef}*{symbol_str}")
        
        # 总是显示 *x^exp，包括 *x^0
        print(f"({' + '.join(term_strs)})*x^{exp}", end="")
        
        # 不是最后一项就打印" + "
        if idx < len(sorted_exps) - 1:
            print(" + ", end="")
        else:
            print()  # 最后一项换行

def verify_format(poly: Polynomial):
    """验证格式是否正确"""
    for exp, expr in poly.items():
        print(f"指数 {exp}: {expr}")
        for term in expr:
            assert isinstance(term, tuple), "Term必须是元组"
            assert len(term) == 2, "Term必须有2个元素"
            coeff, symbols = term
            assert isinstance(coeff, (int, float)), "系数必须是数字"
            assert isinstance(symbols, list), "符号必须是列表"
            print(f"  - 系数: {coeff}, 符号: {symbols}")

def test_polynomial_operations():
    print("="*60)
    print("多项式操作测试")
    print("="*60)
    
    # 测试1: (1*a1)*x^0 + (2*a2)*x^1
    poly1 = {
        0: [(1, [1])],      # x^0 项: 1*a1
        1: [(2, [2])]       # x^1 项: 2*a2
    }
    
    print_polynomial(poly1, "原始多项式1")
    # 输出: (1*a1) + (2*a2)*x^1
    
    # 平方展开
    squared1 = polynomial_square(poly1)
    print_polynomial(squared1, "平方展开后")
    # 输出: (1*a1^2) + (4*a1*a2)*x^1 + (4*a2^2)*x^2
    
    # 积分
    integrated1 = integrate_polynomial(squared1)
    print_polynomial(integrated1, "积分后")
    # 输出: (1*a1^2)*x^1 + (2*a1*a2)*x^2 + (4/3*a2^2)*x^3
    
    # 测试2: (2*a2)*x^0
    poly2 = {
        0: [(2, [2])]  # x^0 项: 2*a2
    }
    
    print_polynomial(poly2, "\n原始多项式2")
    # 输出: (2*a2)
    
    squared2 = polynomial_square(poly2)
    print_polynomial(squared2, "平方展开后")
    # 输出: (4*a2^2)
    
    integrated2 = integrate_polynomial(squared2)
    print_polynomial(integrated2, "积分后")
    # 输出: (4*a2^2)*x^1
    
    # 测试3: 包含多个项的表达式
    poly3 = {
        0: [(1, [1]), (1, [2])],  # x^0 项: a1 + a2
        2: [(3, [3])]             # x^2 项: 3*a3
    }
    
    print_polynomial(poly3, "\n原始多项式3")
    # 输出: (1*a1 + 1*a2) + (3*a3)*x^2
    
    squared3 = polynomial_square(poly3)
    print_polynomial(squared3, "平方展开后")
    # 输出: (1*a1^2 + 2*a1*a2 + 1*a2^2) + (6*a1*a3 + 6*a2*a3)*x^2 + (9*a3^2)*x^4
    
    integrated3 = integrate_polynomial(squared3)
    print_polynomial(integrated3, "积分后")
    
    
def test_integral():
    print("="*60)
    print("测试：平方后的积分")
    print("="*60)
    
    # 原始多项式: (1*a1^2) + (4*a1*a2)*x + (4*a2^2)*x^2
    poly_squared = {
        0: [(1.0, [1, 1])],  # a1^2 * x^0
        1: [(4.0, [1, 2])],  # 4*a1*a2 * x^1
        2: [(4.0, [2, 2])]   # 4*a2^2 * x^2
    }
    
    # 显示原始多项式
    print("原始多项式:")
    #print_polynomial_old_style(poly_squared)
    print_polynomial(poly_squared)
    
    # 计算定积分
    evaluate_and_print(poly_squared, "定积分计算")
    # 期望结果: 1*a1^2 + 0*a1*a2 + (1/3)*a2^2
    
def test_general_integral():
    """测试更复杂的情况"""
    print("\n" + "="*60)
    print("测试：一般多项式积分")
    print("="*60)
    
    # 多项式: (2*a1 + 3*a2)*x^0 + (4*a1*a3)*x^2
    poly = {
        0: [(2.0, [1]), (3.0, [2])],  # (2*a1 + 3*a2)
        2: [(4.0, [1, 3])]             # 4*a1*a3 * x^2
    }
    
    print("原始多项式:")
    #print_polynomial_old_style(poly)
    print_polynomial(poly)
    
    # 在 [-0.5, 0.5] 上积分
    evaluate_and_print(poly, "定积分计算")
    # 期望:
    # x^0 项: (2*a1 + 3*a2) * 1 = 2*a1 + 3*a2
    # x^2 项: 4*a1*a3 * (0.5^3 - (-0.5)^3)/3 = 4*a1*a3 * 0.08333 = 0.3333*a1*a3
    
def test_same_bounds():
    print("="*60)
    print("情况一：多个多项式在同一区间积分后求和")
    print("="*60)
    
    # 多项式1: (1*a1^2) + (4*a1*a2)*x + (4*a2^2)*x^2
    poly1 = {
        0: [(1.0, [1, 1])],
        1: [(4.0, [1, 2])],
        2: [(4.0, [2, 2])]
    }
    
    # 多项式2: (2*a3) + (3*a1*a3)*x
    poly2 = {
        0: [(2.0, [3])],
        1: [(3.0, [1, 3])]
    }
    
    # 多项式3: (5*a2*a3)
    poly3 = {
        0: [(5.0, [2, 3])]
    }
    
    polynomials = [poly1, poly2, poly3]
    
    # 打印原始多项式
    print("\n原始多项式列表:")
    for i, p in enumerate(polynomials, 1):
        print(f"  P{i}(x) = ", end="")
        print_polynomial_old_style(p, "")
    
    # 积分并求和
    print("\n积分求和过程（区间[-1/2, 1/2]）:")
    total_result = sum_integrals_same_bounds(polynomials, -0.5, 0.5)
    
    # 打印最终结果
    print(f"\n最终结果:")
    print(f"Σ ∫ P_i(x) dx = {format_expression(total_result)}")
    
def compute_alpha_beta(row_index: int, r: int) -> Tuple[float, float]:
    """计算第row_index行的积分区间[α, β]"""
    middle = -r + row_index
    return middle - 0.5, middle + 0.5
    
def compute_integral(alpha: float, beta: float, power: int) -> float:
    """计算∫_{α}^{β} ξ^power dξ"""
    if power == 0:
        return beta - alpha
    return (beta**(power + 1) - alpha**(power + 1)) / (power + 1) 
    
def compute_mass_matrix(k: int, r: int) -> np.ndarray:
    """计算k×k的矩阵M（数值矩阵）"""
    M = np.zeros((k, k), dtype=float)
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        for j in range(k):
            M[i, j] = compute_integral(alpha, beta, j)
    return M
    
def create_differential_matrix(k):
    rows = k - 1
    cols = k - 1
    # matrix每个元素存Term + power
    matrix = np.empty((rows, cols), dtype=object)
    
    # 构建matrix并打印导数系数表
    #              x^1, x^2, ..., x^(k-1)
    # d^1/dx^1     1  , 2x^1,...,(k-1)x^(k-2)
    # d^2/dx^2     0  , 2   ,...,(k-2)(k-1)x^(k-3)
    # ....
    # d^k-1/dx^k-1 0   ,0   ,..., k!
    # coef, power = n!/(n-m)!, n-m
    for i in range(rows):
        for j in range(cols):
            #返回 x^n 的 m 阶导数形式 (系数, x的指数)
            coef, power = derivative_form(j + 1, i + 1)
            acoef = j + 1
            
            if coef != 0:
                # 新结构：Term = (系数, [符号下标列表])
                term = (coef, [acoef])
            else:
                term = (0, [])
            
            # 存储 (Term, power) 元组
            matrix[i][j] = (term, power)
    
    #print(f'matrix=\n{matrix}')
    return matrix
    
def build_polynomial_list(matrix, num_rows, num_cols):
    """
    从差分矩阵构建多项式列表。
    每个多项式是一个字典：{power: [terms]}，其中term = (coef, symbols)。
    """
    polynomial_list = []
    for i in range(num_rows):
        polynomial = defaultdict(list)
        for j in range(num_cols):
            term, power = matrix[i][j]
            coef, symbols = term
            if coef != 0:
                polynomial[power].append(term)
        polynomial_list.append(dict(polynomial))
    return polynomial_list
    
    
def compute_squared_polynomials(polynomial_list):
    """
    计算多项式列表中每个多项式的平方。
    返回平方后的多项式列表。
    """
    squared_list = []
    for poly in polynomial_list:
        squared = polynomial_square(poly)
        squared_list.append(squared)
    return squared_list

def print_original_polynomials(squared_polynomials):
    """
    以旧风格打印平方后的多项式列表。
    """
    print("\n原始多项式列表:")
    for i, poly in enumerate(squared_polynomials, 1):
        print(f" P{i}(x) = ", end="")
        print_polynomial_old_style(poly, "")
        
def solve_for_coefficients(M):
    rows, cols = M.shape
    #print(f'rows,cols={rows},{cols}')
    a_coeffs = np.empty((rows, cols), dtype=object)
    for i in range(rows):
        for j in range(cols):
            coeff = M[i, j]
            a_coeffs[i,j] = (coeff, j)
    #print(f'a_coeffs={a_coeffs}')
    return a_coeffs
    
def float_to_fraction_str(num, max_denominator=1000):
    """将浮点数转换为最简分数字符串"""
    frac = Fraction(num).limit_denominator(max_denominator)
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"    
    
def coef_to_str(coeff, id, isfirst):
    csign = '-'
    #print(f'isfirst={isfirst}')
    if coeff >= 0:
        if not isfirst:
            csign = '+'
        else:
            csign = ' '
    else:
        csign = '-'
        
    v_str = f'{csign}{abs(coeff)}*v[{id}]'
    return v_str
    
def id_with_sign(id):
    id_sign = '-'
    if id >= 0:
        id_sign = '+'
    return id_sign, f"{abs(id)}"
    
def coef_to_fraction_str(coeff, id, isfirst):
    csign = '-'
    if coeff >= 0:
        if not isfirst:
            csign = '+'
        else:
            csign = ' '
    else:
        csign = '-'
        
    max_denominator = 1000
    frac = Fraction(abs(coeff)).limit_denominator(max_denominator)
    frac_str = f"{frac.numerator}/{frac.denominator}"
    if frac.denominator == 1:
        frac_str = f"{frac.numerator}"
        
    id_sign, abs_id = id_with_sign(id)
    
    v_str = f"{csign}{frac_str}*v[i{id_sign}{abs_id}]"
    return v_str

def print_coeffs_expression(a_coeffs,k,r):
    rows, cols = a_coeffs.shape
    print(f'\nr={r},k={k}')
    for i in range(rows):
        expr_parts = [f'a[{i}] = ']
        for j in range(cols):
            coeff, id = a_coeffs[i, j]
            #v_str = coef_to_str(coeff, id, j==0)
            v_str = coef_to_fraction_str(coeff, -r+id, j==0)
            expr_parts.append(f'{v_str}')
        expr = ''.join(expr_parts)
        print(f'{expr}')

    #print(f'a_coeffs={a_coeffs}')
    return a_coeffs
    
def print_separator(length=70, char='='):
    """打印指定长度和字符的分隔线"""
    print(char * length)
    
def get_index_range(k, r):
    # Generate index range string    
    if r == 0:
        return f"[i,i+{k-1}]"
    elif r == k-1:
        return f"[i-{k-1}, i]"
    else:
        return f"[i-{r},i+{k-1-r}]"
    
def print_polynomial_coefficients(k, coeffs_list, v_name='v'):
    """
    Print reconstruction coefficients in a professional academic format (English)
    """
    # Print header
    print_separator()
    print(f"Polynomial Reconstruction: p(ξ) = a₀ + a₁ξ + a₂ξ² + ... + aₖ₋₁ξ^{k-1}")
    print(f"Configuration Parameters: k = {k} (Polynomial Degree = {k-1})")
    print_separator()
            
    # Process each r value
    r_values = list(range(k))
    for idx, (r, coeffs) in enumerate(zip(r_values, coeffs_list)):
        print(f"\n[r = {r}] (Stencil Center Offset, Covering Cells {get_index_range(k,r)})")
        print_separator(60, "-")

        M_inv = coeffs  # Assuming input is already M^{-1}
        
        for a_idx in range(k):
            terms = []
            
            for col in range(k):
                coeff, id = M_inv[a_idx, col]
                
                # Skip near-zero coefficients
                if abs(coeff) < 1e-12:
                    continue
                
                # Convert to fraction
                frac = Fraction(coeff).limit_denominator(1000)
                if frac.denominator == 1:
                    coeff_str = str(frac.numerator)
                else:
                    coeff_str = f"{frac.numerator}/{frac.denominator}"
                
                # Handle sign
                if float(coeff) >= 0:
                    sign = " + " if terms else "  "
                else:
                    sign = " - "
                    if coeff_str.startswith('-'):
                        coeff_str = coeff_str[1:]
                
                # Handle coefficient of ±1
                if coeff_str == '1':
                    term = f"{sign}{v_name}[i"
                else:
                    term = f"{sign}{coeff_str}·{v_name}[i"
                
                # Determine index offset
                offset = col - r
                if offset > 0:
                    term += f"+{offset}"
                elif offset < 0:
                    term += f"{offset}"
                term += "]"
                
                terms.append(term)
            
            expression = "".join(terms) if terms else "  0"
            print(f"a{a_idx} = {expression}")
    
    print()
    print_separator()
    
def solve_polynomial_coefficients(k, r):
    M = compute_mass_matrix(k,r)
    #print(f'mass_matrix=\n{M}')
    
    try:
        M_inv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        raise ValueError(f"矩阵M (k={k}, r={r})不可逆！")
    #print(f'M_inv=\n{M_inv}')
    
    a_coeffs = solve_for_coefficients(M_inv)
    return a_coeffs
    
    
def solve_smoothness_indicator(k):
    # 创建差分矩阵
    matrix = create_differential_matrix(k)
    num_rows = k - 1
    num_cols = k - 1
    
    #print(f'差分矩阵:\n{matrix}')
    
    # 从矩阵构建多项式列表
    polynomial_list = build_polynomial_list(matrix, num_rows, num_cols)
    #print(f'\n多项式列表: {polynomial_list}')
    
    # 计算每个多项式的平方
    squared_polynomials = compute_squared_polynomials(polynomial_list)
    #print(f"squared_polynomials={squared_polynomials}")
    
    # 打印平方后的多项式（原始多项式列表）
    print_original_polynomials(squared_polynomials)
    
    # 在指定区间上积分求和
    print("\n积分求和过程（区间[-1/2, 1/2]）:")
    lower_bound = -0.5
    upper_bound = 0.5
    total_result = sum_integrals_same_bounds(squared_polynomials, lower_bound, upper_bound)
    
    # 打印最终结果
    print(f"\n最终结果:")
    formatted_result = format_expression(total_result)
    print(f"Σ ∫ P_i(x) dx = {formatted_result}")
    #print(f'total_result={total_result}')
    return total_result
    
def sort_indices_with_counts(index_list):
    """
    统计下标频次并排序
    
    返回: (排序后的下标列表, 对应的次数列表)
    """
    freq_dict = Counter(index_list)
    sorted_items = sorted(freq_dict.items())
    indices, counts = zip(*sorted_items)  # 解压元组
    return list(indices), list(counts)
    
def polynomial_coefficients_str(coeffs,k,r):
    expr_parts = []
    for j in range(k):
        coeff, id = coeffs[j]
        v_str = coef_to_fraction_str(coeff, -r+id, j==0)
        expr_parts.append(f"{v_str}")
    expr = ''.join(expr_parts)
    print(f'{expr}')
    return expr
    
def get_numeric_list(numbers):
    """
    从元组列表中提取第一个数值元素
    
    参数：
        numbers: list[tuple] - 元组列表，每个元组第一个元素为np.float64
    返回：
        list[np.float64] - 数值列表
    """
    # 列表推导式 + 简单校验，避免索引越界
    return [item[0] for item in numbers if isinstance(item, tuple) and len(item) >= 1]
    
def unpack_tuple_list(numbers):
    print("输入类型:", type(numbers))  # 调试：查看是list还是np.ndarray
    print("输入内容:", numbers)    
    float_list = []
    index_list = []
    # 遍历+类型校验，避免非法数据报错
    for item in numbers:
        if isinstance(item, tuple) and len(item) >= 2:
            float_val, index = item[0], item[1]
            float_list.append(float_val)
            index_list.append(index)
    return float_list, index_list
    
def zip_lists_to_tuples(value_list, index_list):
    """
    极简版合并列表，兼容任意类型（无校验，适合内部可信数据）
    
    参数：
        value_list: list - 任意类型数值列表
        index_list: list - 任意类型索引列表
    返回：
        list[tuple] - 合并后的元组列表
    """
    return list(zip(value_list, index_list))

def format_expression_coefficients(expr, a_coeffs, k, r):
    """格式化纯符号表达式"""
    if not expr:
        return ""
        
    print(f"expr={expr}")
    term_strs = []
    for coeff, symbols in expr:
        indices, counts = sort_indices_with_counts(symbols)
        #print(f"排序下标: {indices}")
        #print(f"出现次数: {counts}")
        
        nSize = len(indices)
        symbol_str = []
        totalfactor = 1
        print(f'counts={counts}')
        for i in range(nSize):
            id = indices[i]
            co = counts[i]
            #a_coeff = get_numeric_list(a_coeffs[id])
            print(f'a_coeffs[id]={a_coeffs[id]}')
            floatlist, idlist = unpack_tuple_list(a_coeffs[id])
            factor, simplified = extract_max_common_factor(floatlist)
            a_coeff_new = zip_lists_to_tuples(simplified, idlist)
            factors = pow(factor, co)
            totalfactor *= factors
            coefficients_str = polynomial_coefficients_str(a_coeff_new,k,r)
            symbol_str.append(f"({coefficients_str})^{co}")
        symbol_str_final = "*".join(symbol_str)
        
        print(f"symbol_str_final: {symbol_str_final}")
        term_strs.append(f"{coeff*factors}*{symbol_str_final}")
    
    return " + ".join(term_strs)
    
def print_coeffs_expression(a_coeffs,k,r):
    rows, cols = a_coeffs.shape
    print(f'\nr={r},k={k}')
    for i in range(rows):
        expr_parts = [f'a[{i}] = ']
        for j in range(cols):
            coeff, id = a_coeffs[i, j]
            #v_str = coef_to_str(coeff, id, j==0)
            v_str = coef_to_fraction_str(coeff, -r+id, j==0)
            expr_parts.append(f'{v_str}')
        expr = ''.join(expr_parts)
        print(f'{expr}')

    #print(f'a_coeffs={a_coeffs}')
    return a_coeffs    
    
def print_smoothness_indicator(expression,a_coeffs,k,r):
    print(f"expression={expression}")
    print(f"Configuration Parameters: k = {k} (Polynomial Degree = {k-1})")
    print(f"\n[r = {r}] (Stencil Center Offset, Covering Cells {get_index_range(k,r)})")
    print(f"expression = {format_expression(expression)}")
    print(f"β{r} = {format_expression(expression)}")
    expr_str = format_expression_coefficients(expression, a_coeffs, k, r)
    print(f"expr_str = {expr_str}")
    
def demo_smoothness_indicator(k):
    total_result = solve_smoothness_indicator(k)
    print(f'total_result={total_result}')
    
    coeffs_list = []
    for r in range(k):
        a_coeffs = solve_polynomial_coefficients(k, r)
        coeffs_list.append( a_coeffs )
        #print(f'a_coeffs={a_coeffs}')
        #print_coeffs_expression(a_coeffs,k,r)
        print_smoothness_indicator(total_result,a_coeffs,k,r)
    print_polynomial_coefficients(k, coeffs_list)

if __name__ == "__main__":
    demo_smoothness_indicator(3)