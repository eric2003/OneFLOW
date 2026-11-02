import numpy as np
from typing import Tuple

def compute_integral(alpha: float, beta: float, power: int) -> float:
    """
    计算∫_{α}^{β} ξ^power dξ的精确值
    """
    if power < 0:
        raise ValueError(f"幂次必须为非负整数，但得到 {power}")
    
    if power == 0:
        return beta - alpha
    
    return (beta**(power + 1) - alpha**(power + 1)) / (power + 1)


def compute_alpha_beta(row_index: int, r: int) -> Tuple[float, float]:
    """
    计算第row_index行的积分区间[α, β]
    """
    middle = -r + row_index
    alpha = middle - 0.5
    beta = middle + 0.5
    
    return alpha, beta


def compute_matrix_M(k: int, r: int) -> np.ndarray:
    """
    计算k×k的矩阵M
    
    矩阵定义: M[i][j] = ∫_{α_i}^{β_i} ξ^j dξ
    """
    print(f"\n开始计算矩阵M (k={k}, r={r})...")
    print("=" * 50)
    
    # 创建全零矩阵
    M = np.zeros((k, k), dtype=float)
    
    # 逐行计算
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        
        print(f"\n第 {i} 行: 区间 [{alpha:.2f}, {beta:.2f}]")
        print("-" * 40)
        
        # 逐列计算
        for j in range(k):
            value = compute_integral(alpha, beta, j)
            M[i, j] = value
            
            print(f"  M[{i}][{j}] = ∫_{alpha:.2f}^{beta:.2f} ξ^{j} dξ = {value:.6f}")
        
        print(f"  第{i}行完成: {M[i, :]}")
    
    print("\n" + "=" * 50)
    print("矩阵计算完成！")
    
    return M


def print_matrix_nicely(M: np.ndarray, r: int, precision: int = 4):
    """
    美化打印矩阵
    """
    k = M.shape[0]
    
    print(f"\n{'='*60}")
    print(f"最终矩阵M (k={k}, r={r})")
    print(f"{'='*60}\n")
    
    # ✅ 修复：suppress_small 改为 suppress
    np.set_printoptions(precision=precision, suppress=True, linewidth=120)
    
    print("矩阵数值:")
    print(M)
    print()
    
    # 打印区间信息
    print("积分区间详情:")
    print(f"{'行号 i':<8} {'α_i':<10} {'β_i':<10} {'宽度 β_i-α_i':<15}")
    print("-" * 45)
    
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        width = beta - alpha
        print(f"{i:<8} {alpha:<10.2f} {beta:<10.2f} {width:<15.2f}")
    
    # 验证第一列
    print(f"\n验证：每行第一列M[i][0] = β_i - α_i = 1.0")
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        expected = beta - alpha
        actual = M[i, 0]
        status = "✓ 正确" if abs(actual - expected) < 1e-6 else "✗ 错误"
        print(f"  第{i}行: M[{i}][0] = {actual:.4f}, 期望值 = {expected:.4f}  {status}")


def demonstrate_3x3_example():
    """
    演示3×3矩阵的例子，r=1
    """
    print("="*70)
    print("演示：计算3×3矩阵M (r=1)")
    print("="*70)
    
    # 计算3×3矩阵
    M = compute_matrix_M(k=3, r=1)
    
    # 美化打印
    print_matrix_nicely(M, r=1, precision=4)


if __name__ == "__main__":
    demonstrate_3x3_example()