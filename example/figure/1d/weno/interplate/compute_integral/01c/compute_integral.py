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
    求解系数a的符号表达式
    
    参数:
        k: 矩阵维度
        r: 参数值 (0 ≤ r < k)
        i: 基础索引（用于构造v的下标）
    
    返回:
        a_coeffs: a0, a1, ..., a_{k-1}的符号表达式列表
    """
    # 1. 计算矩阵M
    M = compute_matrix_M(k, r)
    
    # 2. 计算逆矩阵
    try:
        M_inv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        raise ValueError(f"矩阵M (k={k}, r={r})不可逆！")
    
    # 3. 构造符号向量v的下标
    # v_indices = [i-r+0, i-r+1, ..., i-r+k-1]
    v_indices = [i - r + m for m in range(k)]
    
    # 4. 求解每个a_j
    a_coeffs = []
    
    for j in range(k):  # j是a的下标
        # a_j = Σ(M_inv[j][m] * v_{i-r+m})
        terms = []
        
        for m in range(k):  # m是求和变量
            coeff = M_inv[j, m]  # 数值系数
            
            # 只处理非零系数
            if abs(coeff) > 1e-10:
                # 格式化系数
                if abs(coeff - 1.0) < 1e-10:
                    term_str = f"v[{v_indices[m]}]"
                elif abs(coeff + 1.0) < 1e-10:
                    term_str = f"-v[{v_indices[m]}]"
                else:
                    # 尝试转换为分数显示
                    frac = Fraction(coeff).limit_denominator(1000)
                    if frac.denominator == 1:
                        term_str = f"{frac.numerator}*v[{v_indices[m]}]"
                    else:
                        term_str = f"{frac}*v[{v_indices[m]}]"
                
                terms.append(term_str)
        
        # 合并项
        if not terms:
            a_coeffs.append("0")
        elif len(terms) == 1:
            a_coeffs.append(terms[0])
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
    print(f"v = [v[{i-r+0}], v[{i-r+1}], ..., v[{i-r+k-1}]]")
    
    print(f"\n求解得到的系数 a:")
    print("-" * 50)
    for j, expr in enumerate(a_coeffs):
        print(f"  a_{j} = {expr}")
    
    print("\n向量形式:")
    print("  [a_0, a_1, ..., a_{k-1}]^T = M⁻¹ * [v[i-r+0], v[i-r+1], ..., v[i-r+k-1]]^T")

# ============= 示例：k=3 =============

def example_k_3():
    """k=3的完整示例"""
    k = 3
    r = 1
    i = 0  # 基础索引
    
    print("="*70)
    print("示例：k=3, r=1, i=0")
    print("="*70)
    
    # 计算矩阵M
    M = compute_matrix_M(k, r)
    print("\n步骤1: 计算矩阵M")
    print("-" * 40)
    print("M[i][j] = ∫_{α_i}^{β_i} ξ^j dξ")
    print("\n其中区间：")
    for row in range(k):
        alpha, beta = compute_alpha_beta(row, r)
        print(f"  第{row}行: α_{row}={alpha:.2f}, β_{row}={beta:.2f}")
    
    print("\n得到的矩阵M:")
    print(M)
    
    # 计算逆矩阵
    M_inv = np.linalg.inv(M)
    print("\n步骤2: 计算逆矩阵 M⁻¹")
    print("-" * 40)
    print(M_inv)
    
    # 验证 M * M⁻¹ = I
    identity = M @ M_inv
    print("\n验证 M × M⁻¹ = I:")
    print(identity)
    
    # 求解系数
    print("\n步骤3: 求解 a = M⁻¹ v")
    print("-" * 40)
    print("符号向量 v = [v[i-r+0], v[i-r+1], v[i-r+2]]")
    print(f"         = [v[{0-r+0}], v[{0-r+1}], v[{0-r+2}]]")
    print(f"         = [v[{-r}], v[{-r+1}], v[{-r+2}]]")
    
    a_coeffs = solve_for_coefficients(k, r, i)
    format_coefficients(a_coeffs, k, r, i)
    
    # 额外：展示r=0和r=2的对比
    print("\n" + "="*70)
    print("对比：r=0, 1, 2 的系数表达式")
    print("="*70)
    
    for test_r in [0, 1, 2]:
        if test_r < k:  # r必须小于k
            coeffs = solve_for_coefficients(k, test_r, i)
            print(f"\nr={test_r}:")
            for j, expr in enumerate(coeffs):
                print(f"  a_{j} = {expr}")

if __name__ == "__main__":
    example_k_3()