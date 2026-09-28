import numpy as np
from fractions import Fraction
from typing import List, Tuple

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
    example_k_3()