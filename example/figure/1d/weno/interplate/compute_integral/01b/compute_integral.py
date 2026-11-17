import numpy as np
from fractions import Fraction
from typing import Tuple

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
    """计算k×k的矩阵M"""
    M = np.zeros((k, k), dtype=float)
    
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        for j in range(k):
            M[i, j] = compute_integral(alpha, beta, j)
    
    return M

def matrix_to_fractions(M: np.ndarray) -> list:
    """
    将矩阵转换为Fraction格式，便于美观显示
    
    处理技巧：
    - 分母不超过1000
    - 避免过大分数
    """
    fraction_matrix = []
    for row in M:
        fraction_row = []
        for val in row:
            # 将小数转换为分数，限制分母大小
            frac = Fraction(val).limit_denominator(1000)
            fraction_row.append(frac)
        fraction_matrix.append(fraction_row)
    
    return fraction_matrix

def format_fraction_matrix(frac_matrix: list) -> str:
    """格式化Fraction矩阵为字符串"""
    lines = []
    for row in frac_matrix:
        # 将每个Fraction转换为字符串，并统一宽度
        row_str = "  ".join([f"{str(frac):>10}" for frac in row])
        lines.append(f"[ {row_str} ]")
    
    return "\n".join(lines)

def format_float_matrix(M: np.ndarray, precision: int = 4) -> str:
    """格式化浮点矩阵为字符串"""
    # 设置numpy打印选项
    with np.printoptions(precision=precision, suppress=True, linewidth=100):
        return str(M)

def compute_and_display(k: int, r: int):
    """
    主函数：计算并显示矩阵M及其逆矩阵
    
    参数:
        k: 矩阵维度
        r: 参数值 (0 ≤ r < k)
    """
    print(f"\n{'='*70}")
    print(f"k = {k}, r = {r} 的矩阵M 及其逆矩阵")
    print(f"{'='*70}\n")
    
    # 1. 计算矩阵M
    M = compute_matrix_M(k, r)
    
    # 2. 计算逆矩阵
    try:
        M_inv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        print("⚠️  矩阵不可逆！")
        return
    
    # 3. 转换为Fraction格式
    M_frac = matrix_to_fractions(M)
    M_inv_frac = matrix_to_fractions(M_inv)
    
    # 4. 显示实数格式
    print("📊 实数格式:")
    print("-" * 30)
    print("矩阵 M:")
    print(format_float_matrix(M, precision=6))
    print("\n逆矩阵 M⁻¹:")
    print(format_float_matrix(M_inv, precision=6))
    
    # 5. 显示Fraction格式（美观）
    print(f"\n{'='*70}")
    print("🎨 Fraction分数格式（美观）:")
    print("-" * 30)
    print("矩阵 M:")
    print(format_fraction_matrix(M_frac))
    print("\n逆矩阵 M⁻¹:")
    print(format_fraction_matrix(M_inv_frac))
    
    # 6. 验证 M × M⁻¹ = I
    identity = M @ M_inv
    print(f"\n{'='*70}")
    print("✅ 验证 M × M⁻¹ = I (单位矩阵):")
    print("-" * 30)
    print(format_float_matrix(identity, precision=10))

# ==================== 运行示例 ====================
if __name__ == "__main__":
    # 示例1：3×3矩阵，r=0
    compute_and_display(k=3, r=0)
    
    # 示例2：3×3矩阵，r=1
    compute_and_display(k=3, r=1)
    
    # 示例3：4×4矩阵，r=2
    compute_and_display(k=4, r=2)