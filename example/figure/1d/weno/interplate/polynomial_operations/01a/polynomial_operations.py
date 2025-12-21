from collections import defaultdict
from typing import List, Tuple, Dict
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

def format_expression(expr: Expression) -> str:
    """格式化 Expression"""
    if not expr:
        return "0"
    
    term_strs = [format_term(term) for term in expr]
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

# ============= 测试 =============

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
    
    
def test_verify_format():
    poly1 = {0: [(1, [1])], 1: [(2, [2])]}
    
    squared = polynomial_square(poly1)
    # squared 仍然是 Polynomial 格式：
    # {
    #     0: [(1, [1, 1])],
    #     1: [(4, [1, 2])],
    #     2: [(4, [2, 2])]
    # }    
    
    integrated = integrate_polynomial(squared)
    # integrated 仍然是 Polynomial 格式：
    # {
    #     1: [(1.0, [1, 1])],
    #     2: [(2.0, [1, 2])],
    #     3: [(1.333..., [2, 2])]
    # }
    
    # 测试
    verify_format(poly1)      # ✓ 通过
    verify_format(squared)    # ✓ 通过
    verify_format(integrated) # ✓ 通过    
    
# 运行测试
#test_polynomial_operations()
test_verify_format()