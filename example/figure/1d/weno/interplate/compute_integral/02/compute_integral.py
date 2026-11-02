import numpy as np
from fractions import Fraction
from collections import defaultdict
from typing import List, Tuple, Dict


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
    """
    求解系数a的符号表达式（输出格式：v[i±offset]）
    
    关键改进：负系数前直接显示减号，不显示+ -
    """
    M = compute_matrix_M(k, r)
    try:
        M_inv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        raise ValueError(f"矩阵M (k={k}, r={r})不可逆！")
    
    a_coeffs = []
    
    # 对每个系数a_j
    for j in range(k):
        terms = []
        
        # 对求和项m
        for m in range(k):
            coeff = M_inv[j, m]
            
            if abs(coeff) > 1e-10:  # 非零系数
                # 计算相对偏移量
                offset = m - r
                
                # 格式化符号下标
                if offset == 0:
                    v_str = "v[i]"
                elif offset > 0:
                    v_str = f"v[i+{offset}]"
                else:  # offset < 0
                    v_str = f"v[i{offset}]"
                
                # 格式化系数部分
                if abs(coeff - 1.0) < 1e-10:
                    term_str = v_str  # 系数为1
                elif abs(coeff + 1.0) < 1e-10:
                    term_str = f"-{v_str}"  # 系数为-1
                else:
                    frac = Fraction(coeff).limit_denominator(1000)
                    term_str = f"{frac}*{v_str}"
                
                terms.append(term_str)
        
        # 关键改进：智能连接各项，避免+ -问题
        if not terms:
            a_coeffs.append("0")
        else:
            # 第一个项直接添加（保留其原始符号）
            expr = terms[0]
            
            # 后续项根据符号智能连接
            for term in terms[1:]:
                if term.startswith('-'):
                    # 负号项：直接加空格和项（负号自带）
                    expr += f" {term}"
                else:
                    # 正号项：加 + 和项
                    expr += f" + {term}"
            
            a_coeffs.append(expr)
    
    return a_coeffs
    
# ============ 新增：符号复合功能 ============

def evaluate_polynomial_integral_symbolic(polynomial: Dict[int, List[Tuple[float, List[int]]]], 
                                         a_coeffs: List[str]) -> str:
    """
    对多项式进行积分，并将a_j替换为v表达式
    
    参数:
        polynomial: {指数: [(系数, [符号下标列表])]}
        a_coeffs: a_j的v表达式列表
    
    返回:
        复合表达式字符串（如"1.0*(-1/2*v[i-1] - 1/2*v[i+1])^2"）
    """
    # 在[-0.5, 1/2]上积分
    a, b = -0.5, 0.5
    
    # 用字典累积符号项的系数
    result_dict = defaultdict(float)
    
    for exp, expr_list in polynomial.items():
        integral_factor = (b ** (exp + 1) - a ** (exp + 1)) / (exp + 1)
        
        for coeff, symbols in expr_list:
            # 生成符号键：如"a1*a1"或"a1*a2"
            if len(symbols) == 1:
                symbol_key = f"a{symbols[0]}"  # 如"a1"
            elif symbols[0] == symbols[1]:
                symbol_key = f"a{symbols[0]}^2"  # 如"a1^2"
            else:
                symbol_key = "*".join([f"a{s}" for s in symbols])  # 如"a1*a2"
            
            # 累加积分后的系数
            contribution = coeff * integral_factor
            result_dict[symbol_key] += contribution
    
    # 构建复合表达式
    terms = []
    for symbol_key, total_coeff in result_dict.items():
        if abs(total_coeff) < 1e-10:
            continue
        
        # 获取符号对应的a_j表达式
        # symbol_key如"a1^2" -> 需要找到a1的表达式
        base_symbol = symbol_key.split('*')[0].split('^')[0]  # 提取"a1"
        a_index = int(base_symbol[1:]) - 1  # a1 -> 索引0
        
        if 0 <= a_index < len(a_coeffs):
            a_expr = a_coeffs[a_index]
            
            # 构建项
            if '^2' in symbol_key:
                # 平方项：用括号
                term = f"{total_coeff}*({a_expr})^2"
            else:
                # 一次项：直接用
                term = f"{total_coeff}*{a_expr}"
            
            terms.append(term)
    
    # 智能连接各项
    if not terms:
        return "0"
    
    expr = terms[0]
    for term in terms[1:]:
        if term.startswith('-'):
            expr += f" {term}"
        else:
            expr += f" + {term}"
    
    return expr

def generate_composite_expressions(k: int, polynomial: Dict[int, List[Tuple[float, List[int]]]], i: int = 0):
    """
    生成所有r值的复合表达式f_r
    
    参数:
        k: 矩阵维度
        polynomial: 多项式（用a_j表示）
        i: 基础索引
    
    返回:
        f_r_dict: {r: 复合表达式字符串}
    """
    f_r_dict = {}
    
    print(f"\n生成k={k}的复合表达式f_r")
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
        
        print(f"\nf_{r} = ∫ P(x) dx （代入a系数后）")
        print("-" * 50)
        print(f"f_{r} = {f_r}")
    
    return f_r_dict

# ============= 测试：用你的例子 =============

def test_your_example():
    """
    测试你的例子：P1(x) = (1*a1^2)*x^0 + (4*a1*a2)*x^1 + (4*a2^2)*x^2
              P2(x) = (4*a2^2)*x^0
    """
    k = 3
    
    # 构建多项式（用a1, a2表示）
    # 注意：a1对应索引1，a2对应索引2
    polynomial = {
        0: [(1.0, [1, 1]),   # a1^2
            (4.0, [2, 2])],  # a2^2（来自P2）
        1: [(4.0, [1, 2])],  # a1*a2
        2: [(4.0, [2, 2])]   # a2^2（来自P1的平方项）
    }
    
    print("="*70)
    print("测试：多项式积分后复合a系数表达式")
    print("="*70)
    
    print("\n原始多项式（用a_j表示）:")
    # 简单打印多项式结构
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
    
    # 总结输出
    print("\n" + "="*70)
    print("总结：所有r的复合表达式f_r")
    print("="*70)
    for r, expr in f_r_dict.items():
        print(f"\nf_{r} = {expr}")

def format_coefficients(a_coeffs: List[str], k: int, r: int, i: int = 0):
    """美化打印系数表达式"""
    print(f"\n{'='*70}")
    print(f"k={k}, r={r}, i={i} 的求解结果")
    print(f"{'='*70}\n")
    
    print("矩阵M:")
    M = compute_matrix_M(k, r)
    print(M)
    
    print(f"\n逆矩阵 M⁻¹:")
    M_inv = np.linalg.inv(M)
    print(M_inv)
    
    print(f"\n符号向量 v:")
    print(f"v = [v[i-r+0], v[i-r+1], ..., v[i-r+{k-1}]]")
    
    print(f"\n求解得到的系数 a:")
    print("-" * 50)
    for j, expr in enumerate(a_coeffs):
        print(f"  a_{j} = {expr}")
    
    print("\n向量形式:")
    print("  a = M⁻¹ * v")

def example_k_3():
    """k=3的完整示例"""
    k = 3
    r = 1
    i = 0
    
    print("="*70)
    print("示例：k=3, r=1, i=0（修正符号连接）")
    print("="*70)
    
    a_coeffs = solve_for_coefficients(k, r, i)
    format_coefficients(a_coeffs, k, r, i)
    
    # 对比不同r值
    print("\n" + "="*70)
    print("对比：r=0, 1, 2 的系数表达式（已修正）")
    print("="*70)
    
    for test_r in [0, 1, 2]:
        if test_r < k:
            coeffs = solve_for_coefficients(k, test_r, i)
            print(f"\nr={test_r}:")
            for j, expr in enumerate(coeffs):
                print(f"  a_{j} = {expr}")

if __name__ == "__main__":
    #example_k_3()
    test_your_example()
    