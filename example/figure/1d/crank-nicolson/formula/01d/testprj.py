
from typing import List, Union

def build_tridiag_latex(
    size: int,
    param: str = "r",
    show_first: int = 2,
    show_last: int = 2,
    mid_count: int = 1,
    mid_index: str = "i",
    N_expr: str = "N",
    u_symbol: str = "u",
    superscript: str = "n",
    matrix_name: str = "B",
    vector_name: str = r"\mathbf{u}^n",
    env: str = "pmatrix"
) -> dict:
    """
    动态生成三对角矩阵 B 及其乘积 B·u^n 的 LaTeX 代码。
    
    参数:
        size:        矩阵维度 m = N-1
        param:       参数符号，如 "r", "\\alpha", "\\lambda" 等
        show_first:  开头显示的具体行数
        show_last:   末尾显示的具体行数
        mid_count:   中间连续通用行数（下标连续递增）
        mid_index:   中间第一行的中心索引符号，如 "i", "j", "k"
        N_expr:      末尾表达式，如 "N", "M", "n+1"
        u_symbol:    向量符号，如 "u", "v", "w", "\\phi"
        superscript: 上标，如 "n", "k", "j", "{n+1}", "*"
        matrix_name: 矩阵名称
        vector_name: 向量名称（注意：这里默认用 \\mathbf{u}^n，需手动同步）
        env:         矩阵环境
    """
    p = param
    m = size
    u = u_symbol
    s = superscript
    
    def fmt_u(sub, sup=None):
        """格式化 u 的下标和上标"""
        if sup is None:
            sup = s
        sub_str = str(sub)
        sub_part = f"_{sub_str}" if sub_str.isdigit() else f"_{{{sub_str}}}"
        sup_part = f"^{sup}" if sup else ""
        return f"{u}{sub_part}{sup_part}"
    
    half_p = f"\\dfrac{{{p}}}{{2}}"
    one_minus_p = f"(1-{p})"
    
    # 中间通用行（第 j 行，中心索引递增）
    def mid_row(j):
        if j == 0:
            center, prev, nxt = mid_index, f"{mid_index}-1", f"{mid_index}+1"
        else:
            center = f"{mid_index}+{j}"
            prev = mid_index if j == 1 else f"{mid_index}+{j-1}"
            nxt = f"{mid_index}+{j+1}"
        return (half_p + fmt_u(prev) + " + " +
                one_minus_p + fmt_u(center) + " + " +
                half_p + fmt_u(nxt))
    
    # 末尾行（倒数第 j 行，j 从 show_last 到 1）
    def tail_row(j):
        terms = []
        if j < m:
            terms.append(half_p + fmt_u(f"{N_expr}-{j+1}"))
        terms.append(one_minus_p + fmt_u(f"{N_expr}-{j}"))
        if j > 1:
            terms.append(half_p + fmt_u(f"{N_expr}-{j-1}"))
        return " + ".join(terms)
    
    # 构建 product 行
    product_lines = []
    
    # Head 段
    for i in range(min(show_first, m)):
        terms = []
        if i > 0:   terms.append(half_p + fmt_u(i))
        terms.append(one_minus_p + fmt_u(i + 1))
        if i < m - 1: terms.append(half_p + fmt_u(i + 2))
        product_lines.append(" + ".join(terms))
    
    # Head → Mid 省略号
    if show_first > 0 and mid_count > 0:
        product_lines.append(r"\vdots")
    
    # Mid 段（连续行，内部无省略号）
    for j in range(mid_count):
        product_lines.append(mid_row(j))
    
    # Mid → Tail 省略号
    if show_last > 0 and (show_first > 0 or mid_count > 0):
        product_lines.append(r"\vdots")
    
    # Tail 段
    for j in range(show_last, 0, -1):
        product_lines.append(tail_row(j))
    
    # 构建 B 矩阵
    matrix_lines = []
    for i in range(min(show_first, m)):
        elems = []
        if i > 0:   elems.append(f"\\dfrac{{{p}}}{{2}}")
        elems.append(f"1-{p}")
        if i < m - 1: elems.append(f"\\dfrac{{{p}}}{{2}}")
        matrix_lines.append(" & ".join(elems))
    
    if (show_first > 0 or show_last > 0) and m > max(show_first, show_last):
        matrix_lines.append(r"\ddots & \ddots & \ddots")
    
    for j in range(show_last):
        actual_i = m - show_last + j
        indent = j + 1
        row = [""] * indent
        if actual_i > 0:   row.append(f"\\dfrac{{{p}}}{{2}}")
        row.append(f"1-{p}")
        if actual_i < m - 1: row.append(f"\\dfrac{{{p}}}{{2}}")
        matrix_lines.append(" & ".join(row))
    
    newline = r"\\" + "\n"
    product_body = newline.join(product_lines)
    matrix_body = newline.join(matrix_lines)
    
    product_latex = f"{matrix_name}{vector_name} = \\begin{{{env}}}\n{product_body}\n\\end{{{env}}}"
    matrix_latex = f"{matrix_name} = \\begin{{{env}}}\n{matrix_body}\n\\end{{{env}}}"
    
    return {"matrix": matrix_latex, "product": product_latex}


# ========== 示例 1：默认 u^n ==========
print("=" * 70)
print("示例 1：默认 u^n (size=6)")
print("=" * 70)
result1 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    u_symbol="u",
    superscript="n",
    vector_name=r"\mathbf{u}^n"
)
print(result1["product"])
print()

# ========== 示例 2：换成 v^k ==========
print("=" * 70)
print("示例 2：v^k (size=6)")
print("=" * 70)
result2 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    u_symbol="v",
    superscript="k",
    vector_name=r"\mathbf{v}^k"
)
print(result2["product"])
print()

# ========== 示例 3：换成 \phi^{n+1} ==========
print("=" * 70)
print("示例 3：\\phi^{n+1} (size=6)")
print("=" * 70)
result3 = build_tridiag_latex(
    size=6,
    param=r"\alpha",
    show_first=2,
    show_last=2,
    mid_count=2,
    mid_index="j",
    N_expr="M",
    u_symbol=r"\phi",
    superscript="{n+1}",
    vector_name=r"\boldsymbol{\phi}^{n+1}"
)
print(result3["product"])
print()

# ========== 示例 4：换成 w^* ==========
print("=" * 70)
print("示例 4：w^* (size=6)")
print("=" * 70)
result4 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="m",
    N_expr="N",
    u_symbol="w",
    superscript="*",
    vector_name=r"\mathbf{w}^*"
)
print(result4["product"])
print()

# ========== 示例 5：B 矩阵也同步显示 ==========
print("=" * 70)
print("示例 5：B 矩阵（v^k）")
print("=" * 70)
print(result2["matrix"])
