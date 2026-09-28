# -*- coding: utf-8 -*-
"""
使用 SymPy 生成并操作下面的矩矩阵

M_{j,i} = ∫_{α_j}^{β_j} ξ^i dξ,
α_j = -r + j - 1/2,
β_j = -r + j + 1/2,
j,i = 0,…,k-1
"""

import sympy as sp

def matrix_element(i, j, r):
    """返回 M_{j,i} = (β_j^{i+1} - α_j^{i+1})/(i+1)"""
    half = sp.Rational(1, 2)
    alpha = -r + j - half
    beta  = -r + j + half
    return sp.simplify((beta**(i+1) - alpha**(i+1)) / (i+1))

def moment_matrix(k_val, r_sym):
    """生成 k×k 的矩矩阵 M"""
    M = sp.zeros(k_val, k_val)
    for j in range(k_val):
        for i in range(k_val):
            M[j, i] = matrix_element(i, j, r_sym)
    return M

# ------------------- 示例 -------------------
if __name__ == "__main__":
    # 符号参数
    r = sp.symbols('r', real=True)

    # 1) 生成 3×3 矩阵并打印
    M3 = moment_matrix(3, r)
    print("3×3 矩矩阵 M (符号 r):")
    sp.pprint(M3)

    # 2) 取 r = 2 得到数值矩阵
    M3_num = M3.subs(r, 2)
    print("\n当 r = 2 时的数值矩阵:")
    sp.pprint(M3_num)

    # 3) 计算行列式
    det_M3 = sp.simplify(M3.det())
    print("\n3×3 矩阵的行列式 (化简后):")
    sp.pprint(det_M3)

    # 4) 生成 4×4 矩阵并查看行列式
    M4 = moment_matrix(4, r)
    print("\n4×4 矩矩阵 M (符号 r):")
    sp.pprint(M4)
    det_M4 = sp.simplify(M4.det())
    print("\n4×4 矩阵的行列式 (化简后):")
    sp.pprint(det_M4)