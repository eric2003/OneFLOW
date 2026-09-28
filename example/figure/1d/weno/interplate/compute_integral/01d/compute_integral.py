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
    
    核心修改：将v[-1], v[0]改为v[i-1], v[i+0]的形式
    """
    # 计算矩阵M和逆矩阵
    M = compute_matrix_M(k, r)
    try:
        M_inv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        raise ValueError(f"矩阵M (k={k}, r={r})不可逆！")
    
    a_coeffs = []
    
    # 对每个系数a_j（j是a的下标）
    for j in range(k):
        terms = []
        
        # 对求和项m（m是v的下标偏移）
        for m in range(k):
            coeff = M_inv[j, m]  # 数值系数
            
            if abs(coeff) > 1e-10:  # 只处理非零系数
                # 计算v的下标相对于i的偏移量
                # v_index = i - r + m
                # offset = v_index - i = m - r
                offset = m - r  # 这是核心！不再是具体的数值下标
                
                # 根据偏移量格式化符号下标
                if offset == 0:
                    v_str = f"v[i]"  # 偏移为0
                elif offset > 0:
                    v_str = f"v[i+{offset}]"  # 正偏移
                else:  # offset < 0
                    v_str = f"v[i{offset}]"  # 负偏移（offset自带负号）
                
                # 格式化系数
                if abs(coeff - 1.0) < 1e-10:
                    term_str = v_str  # 系数为1，不显示
                elif abs(coeff + 1.0) < 1e-10:
                    term_str = f"-{v_str}"  # 系数为-1，只显示负号
                else:
                    # 分数格式化
                    frac = Fraction(coeff).limit_denominator(1000)
                    term_str = f"{frac}*{v_str}"
                
                terms.append(term_str)
        
        # 合并项
        if not terms:
            a_coeffs.append("0")
        else:
            a_coeffs.append(" + ".join(terms))
    
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
    if i == 0:  # 简化显示
        print(f"  = [v[-{r}], v[-{r}+1], ..., v[{k-1-r}]]")
    
    print(f"\n求解得到的系数 a:")
    print("-" * 50)
    for j, expr in enumerate(a_coeffs):
        print(f"  a_{j} = {expr}")
    
    print("\n向量形式:")
    print("  [a_0, a_1, ..., a_{k-1}]^T = M⁻¹ * [v[i-r], v[i-r+1], ..., v[i-r+k-1]]^T")

def example_k_3():
    """k=3的完整示例"""
    k = 3
    r = 1
    i = 0
    
    print("="*70)
    print("示例：k=3, r=1, i=0（输出格式为v[i±offset]）")
    print("="*70)
    
    # 求解系数
    a_coeffs = solve_for_coefficients(k, r, i)
    format_coefficients(a_coeffs, k, r, i)
    
    # 对比不同r值
    print("\n" + "="*70)
    print("对比：r=0, 1, 2 的系数表达式")
    print("="*70)
    
    for test_r in [0, 1, 2]:
        if test_r < k:
            coeffs = solve_for_coefficients(k, test_r, i)
            print(f"\nr={test_r}:")
            for j, expr in enumerate(coeffs):
                print(f"  a_{j} = {expr}")

if __name__ == "__main__":
    example_k_3()