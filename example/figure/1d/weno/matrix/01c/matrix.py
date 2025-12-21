# -*- coding: utf-8 -*-
"""
计算重构系数矩阵：coeff_matrix = φᵀ(ξ₀)·M⁻¹
其中 M 为矩矩阵，ξ₀ 为指定位置（如 ±1/2）
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

def compute_coeff_matrix(k_val, xi0=sp.Rational(1, 2)):
    """
    计算系数矩阵（3×k），每行对应 r = 0,1,2 时的 φᵀ(ξ₀)·M⁻¹
    
    参数
    ----------
    k_val : int
        矩阵阶数
    xi0 : sympy.Rational
        求值位置（默认为 1/2）
    """
    r = sp.symbols('r', real=True)
    M_sym = moment_matrix(k_val, r)
    phi_T = sp.Matrix([xi0**i for i in range(k_val)]).T  # φᵀ(ξ₀)
    #print(f"phi_T={phi_T}")
    #print(sp.latex(phi_T, mat_str='matrix', mat_delim='['))

    rows = []
    for r_val in [0, 1, 2]:
        M_inv = M_sym.subs(r, r_val).inv()
        c_row = phi_T * M_inv
        rows.append(c_row)

    coeff_matrix = sp.Matrix(rows)  # 推荐变量名
    return coeff_matrix
    
def latex_matrix(mat, env='matrix'):
    """将 sympy 矩阵转为 LaTeX 字符串"""
    return sp.latex(mat, mat_str=env, mat_delim='[')
    
# ------------------- 示例 -------------------
if __name__ == "__main__":
    # 默认 ξ₀ = 1/2
    print("%% 当 ξ₀ = 1/2 时：")
    
    C2 = compute_coeff_matrix(2)
    print("\n%% k = 2 的系数矩阵：")
    print(latex_matrix(C2))

    C3 = compute_coeff_matrix(3)
    print("\n%% k = 3 的系数矩阵：")
    print(latex_matrix(C3))

    # 示例：ξ₀ = -1/2
    print("\n\n%% 当 ξ₀ = -1/2 时：")
    
    C2_neg = compute_coeff_matrix(2, xi0=sp.Rational(-1, 2))
    print("\n%% k = 2 的系数矩阵：")
    print(latex_matrix(C2_neg))

    C3_neg = compute_coeff_matrix(3, xi0=sp.Rational(-1, 2))
    print("\n%% k = 3 的系数矩阵：")
    print(latex_matrix(C3_neg))