import sympy as sp

def moment_matrix_k3(r_val):
    """生成 k=3 时的矩阵，并返回 LaTeX 字符串"""
    k = 3
    r = sp.symbols('r')
    half = sp.Rational(1, 2)
    
    # 创建矩阵
    M = sp.zeros(k, k)
    for j in range(k):
        for i in range(k):
            alpha = -r + j - half
            beta  = -r + j + half
            M[j, i] = (beta**(i+1) - alpha**(i+1)) / (i+1)
    
    # 代入具体的 r 值
    M_substituted = M.subs(r, r_val)
    
    # 转换为 LaTeX，使用简化模式
    latex_code = sp.latex(M_substituted, mat_str='matrix', mat_delim='[')
    return latex_code

# 生成 r = 0, 1, 2 的 LaTeX 代码
for r in [0, 1, 2]:
    latex = moment_matrix_k3(r)
    print(f"\n%% 当 r = {r} 时的矩阵 (k=3):")
    print(latex)