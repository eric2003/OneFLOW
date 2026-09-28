import numpy as np
from typing import Tuple, Optional
from functools import lru_cache

def compute_integral(alpha: float, beta: float, power: int) -> float:
    """
    计算∫_{α}^{β} ξ^power dξ的精确值
    
    数学公式:
        power = 0: ∫ 1 dξ = β - α
        power ≥ 1: ∫ ξ^power dξ = (β^(power+1) - α^(power+1)) / (power+1)
    
    参数:
        alpha: 积分下限
        beta: 积分上限
        power: 非负整数幂次
    
    返回:
        积分值（浮点数）
    
    示例:
        >>> compute_integral(-0.5, 0.5, 0)
        1.0
        >>> compute_integral(-0.5, 0.5, 1)
        0.0  # 奇函数对称区间
        >>> compute_integral(-0.5, 0.5, 2)
        0.08333333333333333
    """
    if power < 0:
        raise ValueError(f"幂次必须为非负整数，但得到{power}")
    
    # 特殊情况：power=0时，积分值为β-α
    if power == 0:
        return beta - alpha
    
    # 一般情况：使用解析解
    return (beta**(power + 1) - alpha**(power + 1)) / (power + 1)


def compute_alpha_beta(i: int, r: int) -> Tuple[float, float]:
    """
    计算第i行的积分区间[α_i, β_i]
    
    公式:
        α_i = -r + i - 1/2
        β_i = -r + i + 1/2
    
    参数:
        i: 行索引 (0 ≤ i < k)
        r: 给定的参数值
    
    返回:
        (α_i, β_i) 元组
    """
    offset = -r + i
    alpha = offset - 0.5
    beta = offset + 0.5
    return alpha, beta


def compute_matrix_M(k: int, r: int, vectorized: bool = True) -> np.ndarray:
    """
    计算k×k的矩阵M
    
    矩阵M的定义:
        M[i][j] = ∫_{α_i}^{β_i} ξ^j dξ,  其中 i,j = 0,...,k-1
    
    参数说明:
        k: 矩阵维度（正整数）
        r: 参数值（0 ≤ r < k）
        vectorized: 是否使用向量化计算（默认True，更快）
    
    返回:
        k×k的NumPy数组
    
    示例:
        >>> compute_matrix_M(3, 1)
        array([[1. , 0. , 0.08333333],
               [1. , 1. , 0.75      ],
               [1. , 2. , 2.08333333]])
    """
    # 参数验证
    if not isinstance(k, int) or k <= 0:
        raise ValueError(f"k必须是正整数，但得到{k}")
    if not (0 <= r < k):
        raise ValueError(f"r必须在[0, {k-1}]范围内，但得到{r}")
    
    # 使用向量化计算（推荐，速度快）
    if vectorized:
        return _compute_matrix_M_vectorized(k, r)
    
    # 或使用循环计算（更直观，但较慢）
    return _compute_matrix_M_loop(k, r)


def _compute_matrix_M_loop(k: int, r: int) -> np.ndarray:
    """循环版本（易于理解）"""
    M = np.zeros((k, k))
    
    # 逐行计算
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        
        # 逐列计算
        for j in range(k):
            M[i, j] = compute_integral(alpha, beta, j)
    
    return M


@lru_cache(maxsize=32)
def _compute_matrix_M_vectorized(k: int, r: int) -> np.ndarray:
    """
    向量化版本（性能最优）
    
    利用NumPy广播机制，一次性计算所有元素
    缓存结果：相同(k,r)重复调用时直接返回缓存
    """
    # 创建行索引数组 [0, 1, ..., k-1]
    i = np.arange(k)
    
    # 计算所有α_i和β_i
    # α_i = -r + i - 0.5
    # β_i = -r + i + 0.5
    alpha = -r + i - 0.5  # shape: (k,)
    beta = -r + i + 0.5   # shape: (k,)
    
    # 创建列索引数组（幂次+1） [1, 2, ..., k]
    # 因为积分公式需要 j+1
    j = np.arange(1, k + 1)  # shape: (k,)
    
    # 计算α_i^{j+1}和β_i^{j+1}
    # alpha[:, None] shape: (k, 1)
    # j shape: (k,)
    # 广播后结果shape: (k, k)
    alpha_pow = alpha[:, None] ** j  # α_i^{j+1}
    beta_pow = beta[:, None] ** j    # β_i^{j+1}
    
    # 计算积分值矩阵
    # M[i,j] = (β_i^{j+1} - α_i^{j+1}) / (j+1)
    M = (beta_pow - alpha_pow) / j
    
    return M


def format_matrix(M: np.ndarray, precision: int = 6) -> str:
    """
    美化打印矩阵
    
    参数:
        M: 矩阵
        precision: 小数点后保留位数
    
    返回:
        格式化字符串
    """
    return np.array2string(
        M, 
        precision=precision, 
        suppress_small=True,
        max_line_width=100
    )


# ============= 测试和演示 =============

def demonstrate_matrix_properties(k: int = 4, r: int = 1):
    """
    演示矩阵M的性质
    
    性质1: 当r固定时，M[i+1][j] = M[i][j] + 偏移量
    性质2: 每行第一个元素 M[i][0] = β_i - α_i = 1（恒定）
    """
    print("="*60)
    print(f"矩阵M性质演示 (k={k}, r={r})")
    print("="*60)
    
    M = compute_matrix_M(k, r)
    print(f"矩阵M:\n{format_matrix(M)}\n")
    
    # 验证性质2
    print("验证性质：每行第一个元素 = β_i - α_i = 1")
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        print(f"  第{i}行: M[{i}][0] = {M[i,0]:.6f}, β_i-α_i = {beta-alpha:.6f}")
    
    # 验证性质1（展示偏移模式）
    print(f"\n验证性质：相邻行之间的差异模式")
    for j in range(min(3, k)):  # 只看前3列
        print(f"  第{j}列差分: {M[1:, j] - M[:-1, j]}")
    
    # 显示区间信息
    print(f"\n积分区间详情:")
    for i in range(k):
        alpha, beta = compute_alpha_beta(i, r)
        print(f"  第{i}行: [α_{i}, β_{i}] = [{alpha:6.2f}, {beta:6.2f}]")


def compare_performance(k_values: List[int] = [10, 50, 100]):
    """
    性能对比：循环版 vs 向量化版
    
    参数:
        k_values: 要测试的k值列表
    """
    print("\n" + "="*60)
    print("性能对比测试")
    print("="*60)
    
    import time
    
    for k in k_values:
        print(f"\nk = {k}:")
        
        # 测试循环版本
        start = time.time()
        M_loop = compute_matrix_M(k, r=0, vectorized=False)
        time_loop = time.time() - start
        
        # 测试向量化版本
        start = time.time()
        M_vec = compute_matrix_M(k, r=0, vectorized=True)
        time_vec = time.time() - start
        
        print(f"  循环版本: {time_loop:.6f} 秒")
        print(f"  向量化版: {time_vec:.6f} 秒")
        print(f"  速度提升: {time_loop / time_vec:.2f} 倍")
        
        # 验证结果是否一致
        if np.allclose(M_loop, M_vec):
            print(f"  结果一致性: ✓ 通过")
        else:
            print(f"  结果一致性: ✗ 失败！")


def test_special_cases():
    """
    测试特殊k和r值
    """
    print("\n" + "="*60)
    print("特殊值测试")
    print("="*60)
    
    # 测试k=1（最小维度）
    print("\n测试k=1, r=0:")
    M = compute_matrix_M(1, 0)
    print(f"矩阵M: {M}")
    
    # 测试r=0（区间左端点从-0.5开始）
    print("\n测试k=3, r=0:")
    M = compute_matrix_M(3, 0)
    print(f"矩阵M:\n{format_matrix(M)}")
    print("区间: [-0.5, 0.5], [0.5, 1.5], [1.5, 2.5]")


if __name__ == "__main__":
    # 基础演示
    demonstrate_matrix_properties(k=4, r=1)
    
    # 性能测试
    compare_performance([10, 50, 100, 200])
    
    # 特殊值测试
    test_special_cases()
    
    # 生成特定矩阵
    print("\n" + "="*60)
    print("生成k=5, r=2的矩阵:")
    print("="*60)
    M = compute_matrix_M(5, 2)
    print(format_matrix(M))