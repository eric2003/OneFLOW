from collections import defaultdict
from typing import List, Tuple, Dict
import numpy as np
import math

# 类型别名
Term = Tuple[int, List[int]]  # (系数, 符号下标列表)
Expression = List[Term]        # Term 的列表
Polynomial = Dict[int, Expression]  # {指数: 表达式}

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
        #        exp=1 → (0.5^2 - (-0.5)^2)/2 = (0.25 - 0.25)/2 = 0
        #        exp=2 → (0.5^3 - (-0.5)^3)/3 = (0.125 + 0.125)/3 = 0.25/3 ≈ 0.08333
        
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

def demo_smoothness_indicator_old():
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

def demo_smoothness_indicator():
    k = 3
    rows = k - 1
    cols = k - 1
    
    # matrix每个元素存Term + power
    matrix = np.empty((rows, cols), dtype=object)
    
    # 构建matrix并打印导数系数表
    print("\n=== 导数系数表 (coef*x^power) ===")
    for i in range(rows):
        for j in range(cols):
            coef, power = derivative_form(j + 1, i + 1)
            acoef = j + 1
            
            if coef != 0:
                # 新结构：Term = (系数, [符号下标列表])
                term = (coef, [acoef])
            else:
                term = (0, [])
            
            # 存储 (Term, power) 元组
            matrix[i][j] = (term, power)
            print(f"{coef}x^{power}", end=' ')
        print()
    
    print(f'\nmatrix=\n{matrix}')
    
    # 从matrix构建Polynomial List
    print("\n=== 符号表达式表 (coef*a{acoef}*x^power) ===")
    polynomial_list = []
    for i in range(rows):
        # Polynomial = {指数: Expression}
        polynomial = defaultdict(list)
        
        for j in range(cols):
            term, power = matrix[i][j]
            coef, symbols = term
            
            if coef != 0:
                # 添加到对应指数的项列表
                polynomial[power].append(term)
                print(f"{coef}*a{symbols[0]}*x^{power}", end=' ')
            else:
                print("0 ", end=' ')
        
        polynomial_list.append(dict(polynomial))
        print()
    
    print(f'\npolynomial_list={polynomial_list}')
    
    for i, poly in enumerate(polynomial_list):
        print(f"-- 第{i+1}行导数 --")
        print_polynomial(poly)
        squared = polynomial_square(poly)
        print_polynomial(squared, "平方展开后")
        integrated = integrate_polynomial(squared)
        print_polynomial(integrated, "积分后")
            
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
    
if __name__ == "__main__":
    #demo_smoothness_indicator()
    test_same_bounds()
