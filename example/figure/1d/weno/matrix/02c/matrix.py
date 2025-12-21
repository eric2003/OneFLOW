import sympy as sp
from sympy import symbols, Function, Rational, latex, Matrix, Integral

# ---------- 全局符号 ----------
k = sp.Symbol('k', integer=True, positive=True)
r, j, xi = symbols('r j \\xi')
i_sym = symbols('i')

# ---------- 辅助函数 ----------
def vbar(idx):
    """生成 \overline{v}_{i - r + idx}"""
    return sp.Symbol(r"\overline{v}_{" + latex(i_sym - r + idx) + r"}")

def alpha(j): return -r + j - Rational(1, 2)
def beta(j):  return -r + j + Rational(1, 2)

# ---------- 1. 生成通用矩阵 M 的 LaTeX（带 \cdots） ----------
def matrix_M_latex_generic():
    """生成 M = [ ∫ ξ^m dξ ]_{j,m=0}^{k-1} 的通用 LaTeX（含省略号）"""
    rows = []
    # 第一行：j=0
    first_row = " & ".join([
        r"\int_{\alpha_{" + "0" + r"}}^{\beta_{" + "0" + r"}} d\xi",
        r"\int_{\alpha_{" + "0" + r"}}^{\beta_{" + "0" + r"}} \xi^{1} d\xi",
        r"\cdots",
        r"\int_{\alpha_{" + "0" + r"}}^{\beta_{" + "0" + r"}} \xi^{k-1} d\xi"
    ])
    rows.append(first_row)
    
    # 第二行（示例 j=1）
    second_row = " & ".join([
        r"\int_{\alpha_{" + "1" + r"}}^{\beta_{" + "1" + r"}} d\xi",
        r"\int_{\alpha_{" + "1" + r"}}^{\beta_{" + "1" + r"}} \xi^{1} d\xi",
        r"\cdots",
        r"\int_{\alpha_{" + "1" + r"}}^{\beta_{" + "1" + r"}} \xi^{k-1} d\xi"
    ])
    rows.append(second_row)
    
    # 省略行
    rows.append(r"\vdots & \vdots & \ddots & \vdots")
    
    # 最后一行：j = k-1
    last_idx = "k-1"
    last_row = " & ".join([
        r"\int_{\alpha_{" + last_idx + r"}}^{\beta_{" + last_idx + r"}} d\xi",
        r"\int_{\alpha_{" + last_idx + r"}}^{\beta_{" + last_idx + r"}} \xi^{1} d\xi",
        r"\cdots",
        r"\int_{\alpha_{" + last_idx + r"}}^{\beta_{" + last_idx + r"}} \xi^{k-1} d\xi"
    ])
    rows.append(last_row)
    
    matrix_body = " \\\\\n".join(rows)
    return r"M = \begin{bmatrix}" + "\n" + matrix_body + "\n\\end{bmatrix}"

# ---------- 2. 生成具体 M 矩阵（用于计算和 LaTeX） ----------
def build_matrix_M(k_val, r_val=None):
    """
    构建具体 k 下的 M 矩阵。
    若 r_val 为 None，则保留 r 为符号；否则代入具体 r。
    """
    M = sp.zeros(k_val, k_val)
    for j_idx in range(k_val):
        a_j = alpha(j_idx)
        b_j = beta(j_idx)
        if r_val is not None:
            a_j = a_j.subs(r, r_val)
            b_j = b_j.subs(r, r_val)
        for m in range(k_val):
            # ∫ ξ^m dξ from α to β
            integral = (b_j**(m+1) - a_j**(m+1)) / Rational(m+1)
            M[j_idx, m] = sp.simplify(integral)
    return M

# ---------- 3. 生成向量 LaTeX ----------
def vector_a_latex():
    elements = [f"a_{{{i}}}" for i in range(k-1)]  # 仅用于显示
    # 但 k 是符号，无法 range(k)，所以用省略号
    return r"\mathbf{a} = \begin{bmatrix} a_0 \\ a_1 \\ \\vdots \\ a_{k-1} \\end{bmatrix}"

def vector_v_latex():
    # 使用 vbar(j) for j=0,...,k-1
    first = r"\overline{v}_{" + latex(i_sym - r + 0) + r"}"
    second = r"\overline{v}_{" + latex(i_sym - r + 1) + r"}"
    last = r"\overline{v}_{" + latex(i_sym - r + (k - 1)) + r"}"
    return f"\\mathbf{{v}} = \\begin{bmatrix} {first} \\\\ {second} \\\\ \\vdots \\\\ {last} \\end{{bmatrix}}"

# ---------- 4. 生成 a = M^{-1} v 的公式 ----------
def equation_system_latex():
    a_vec = r"\begin{bmatrix} a_0 \\ a_1 \\ \\vdots \\ a_{k-1} \\end{bmatrix}"
    v_vec = r"\begin{bmatrix} \overline{v}_{" + latex(i_sym - r + 0) + r"} \\ \overline{v}_{" + latex(i_sym - r + 1) + r"} \\ \\vdots \\ \overline{v}_{" + latex(i_sym - r + (k - 1)) + r"} \\end{bmatrix}"
    return f"{a_vec} = M^{{-1}} {v_vec}"

# ---------- 主输出：完整 LaTeX 块 ----------
def generate_extended_latex():
    lines = []

    # 向量定义
    lines.append(r"\mathbf{a}=[a_{0},a_{1},\dots,a_{k-1}]^T,\quad \mathbf{\phi}(\xi) = [\xi^0, \xi^1, \dots, \xi^{k-1}]^{T}")
    lines.append(r"\mathbf{v}=[\overline{v}_{" + latex(i_sym - r + 0) + r"},\overline{v}_{" + latex(i_sym - r + 1) + r"},\dots,\overline{v}_{" + latex(i_sym - r + (k - 1)) + r"}]^T")
    
    # 方程
    lines.append(r"\mathbf{a}=M^{-1}\mathbf{v}")
    
    # 展开方程
    #a_vec = r"\begin{bmatrix} a_0 \\ a_1 \\ \\vdots \\ a_{k-1} \\end{bmatrix}"
    a_vec = r"\begin{bmatrix} a_0 \\ a_1 \\ \vdots \\ a_{k-1} \end{bmatrix}"
    #v_vec = r"\begin{bmatrix} \overline{v}_{" + latex(i_sym - r + 0) + r"} \\ \overline{v}_{" + latex(i_sym - r + 1) + r"} \\ \\vdots \\ \overline{v}_{" + latex(i_sym - r + (k - 1)) + r"} \\end{bmatrix}"
    v_vec = r"\begin{bmatrix} \overline{v}_{" + latex(i_sym - r + 0) + r"} \\ \overline{v}_{" + latex(i_sym - r + 1) + r"} \\ \vdots \\ \overline{v}_{" + latex(i_sym - r + (k - 1)) + r"} \end{bmatrix}"
    lines.append(f"{a_vec} = M^{{-1}} {v_vec}")
    
    # 矩阵 M（通用形式）
    lines.append(matrix_M_latex_generic())
    
    return "\\begin{array}{l}\n" + " \\\\\n".join(lines) + "\n\\end{array}"

# ---------- 示例：打印通用公式 ----------
if __name__ == "__main__":
    print(generate_extended_latex())

    # 示例：打印 k=3, r=0 时的具体 M 矩阵
    print("\n\n=== Example: M for k=3, r=0 ===")
    M3 = build_matrix_M(3, r_val=0)
    print(latex(M3))