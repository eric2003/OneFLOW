# -*- coding: utf-8 -*-
"""
计算 c = φᵀ(½)·M⁻¹，其中
  M_{j,i}=∫_{α_j}^{β_j} ξ^i dξ,
  α_j = -r + j - ½, β_j = -r + j + ½,
  φᵀ(½) = [1, ½, …, (½)^{k-1}].

对 r = 0,1,2 计算对应的行向量 c，组合成矩阵 C 并输出 LaTeX。
"""

import sympy as sp

def moment_matrix(k_val, r_sym):
    """生成 k×k 的矩矩阵 M（符号形式）"""
    half = sp.Rational(1, 2)
    M = sp.zeros(k_val, k_val)
    for j in range(k_val):
        alpha = -r_sym + j - half
        beta  = -r_sym + j + half
        for i in range(k_val):
            M[j, i] = (beta**(i + 1) - alpha**(i + 1)) / (i + 1)
    return M

def compute_c_matrix(k_val):
    """
    计算 C 矩阵（3×k），其每一行对应 r = 0,1,2 时的 c = φᵀ(½)·M⁻¹
    """
    # 把 r 视为符号，后面再代入具体数值
    r = sp.symbols('r', real=True)

    # 基础的符号矩阵 M
    M_sym = moment_matrix(k_val, r)

    # φᵀ(½) 行向量
    phi_T = sp.Matrix([ (sp.Rational(1, 2))**i for i in range(k_val) ]).T

    rows = []
    for r_val in [0, 1, 2]:
        # 代入当前 r，得到数值矩阵并求逆
        M_num = M_sym.subs(r, r_val)
        M_inv = M_num.inv()
        # c = φᵀ(½)·M⁻¹
        c_row = phi_T * M_inv
        rows.append(c_row)

    # 纵向堆叠成 3×k 矩阵
    C = sp.Matrix(rows)
    return C

def latex_matrix(mat, env='bmatrix'):
    """将 sympy 矩阵转为 LaTeX 字符串"""
    return sp.latex(mat, mat_str=env, mat_delim='[')

# ------------------- 示例 -------------------
if __name__ == "__main__":
    # k = 2
    C2 = compute_c_matrix(2)
    print("%% k = 2 时的 C 矩阵（LaTeX）")
    print(latex_matrix(C2))

    # k = 3
    C3 = compute_c_matrix(3)
    print("\n%% k = 3 时的 C 矩阵（LaTeX）")
    print(latex_matrix(C3))