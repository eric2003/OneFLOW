
from typing import List, Union

def build_tridiag_latex(
    size: int,
    param: str = "r",
    show_first: int = 2,
    show_last: int = 2,
    mid_count: int = 1,
    mid_index: str = "i",
    N_expr: str = "N",
    matrix_name: str = "B",
    vector_name: str = r"\mathbf{u}^n",
    env: str = "pmatrix"
) -> dict:
    """
    动态生成三对角矩阵 B 及其乘积 B·u^n 的 LaTeX 代码。
    
    中间行下标连续递增（i, i+1, i+2...），相邻中间行之间无省略号。
    省略号只出现在不同段落之间（head→mid, mid→tail）。
    
    参数:
        size:       矩阵维度 m = N-1
        param:      参数符号，如 "r", "\\alpha", "\\lambda" 等
        show_first: 开头显示的具体行数
        show_last:  末尾显示的具体行数
        mid_count:  中间连续通用行数（下标连续递增）
        mid_index:  中间第一行的中心索引符号，如 "i", "j", "k"
        N_expr:     末尾表达式，如 "N", "M", "n+1"
        matrix_name: 矩阵名称
        vector_name: 向量名称
        env:        矩阵环境
    """
    p = param
    m = size
    
    def fmt_u(sub, sup="n"):
        """格式化 u 的下标：纯数字不加花括号，否则加花括号"""
        s = str(sub)
        return f"u_{s}^{sup}" if s.isdigit() else f"u_{{{s}}}^{sup}"
    
    half_p = f"\\dfrac{{{p}}}{{2}}"
    one_minus_p = f"(1-{p})"
    
    # ========== B·u^n 向量 ==========
    
    # 开头具体行（0-based 行号 i）
    def head_row(i):
        terms = []
        if i > 0:
            terms.append(half_p + fmt_u(i))
        terms.append(one_minus_p + fmt_u(i + 1))
        if i < m - 1:
            terms.append(half_p + fmt_u(i + 2))
        return " + ".join(terms)
    
    # 中间通用行（第 j 行，中心索引递增）
    def mid_row(j):
        if j == 0:
            center, prev, nxt = mid_index, f"{mid_index}-1", f"{mid_index}+1"
        else:
            center = f"{mid_index}+{j}"
            prev = mid_index if j == 1 else f"{mid_index}+{j-1}"
            nxt = f"{mid_index}+{j+1}"
        return (
            half_p + fmt_u(prev) + " + " +
            one_minus_p + fmt_u(center) + " + " +
            half_p + fmt_u(nxt)
        )
    
    # 末尾行（倒数第 j 行，j 从 show_last 到 1）
    def tail_row(j):
        terms = []
        if j < m:  # 有前一个元素
            terms.append(half_p + fmt_u(f"{N_expr}-{j+1}"))
        terms.append(one_minus_p + fmt_u(f"{N_expr}-{j}"))
        if j > 1:  # 有后一个元素
            terms.append(half_p + fmt_u(f"{N_expr}-{j-1}"))
        return " + ".join(terms)
    
    # 构建 product 行
    product_lines = []
    
    # Head 段
    for i in range(min(show_first, m)):
        product_lines.append(head_row(i))
    
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
    
    # ========== B 矩阵 ==========
    
    matrix_lines = []
    
    # Head 段
    for i in range(min(show_first, m)):
        elems = []
        if i > 0:
            elems.append(f"\\dfrac{{{p}}}{{2}}")
        elems.append(f"1-{p}")
        if i < m - 1:
            elems.append(f"\\dfrac{{{p}}}{{2}}")
        matrix_lines.append(" & ".join(elems))
    
    # Head → Mid/Tail 省略号
    if (show_first > 0 or show_last > 0) and m > max(show_first, show_last):
        matrix_lines.append(r"\ddots & \ddots & \ddots")
    
    # Tail 段（带缩进）
    for j in range(show_last):
        actual_i = m - show_last + j  # 实际行号
        indent = j + 1  # 前导空 & 的数量
        
        row_parts = [""] * indent
        if actual_i > 0:
            row_parts.append(f"\\dfrac{{{p}}}{{2}}")
        row_parts.append(f"1-{p}")
        if actual_i < m - 1:
            row_parts.append(f"\\dfrac{{{p}}}{{2}}")
        
        matrix_lines.append(" & ".join(row_parts))
    
    # 组装 LaTeX
    newline = r"\\" + "\n"
    product_body = newline.join(product_lines)
    matrix_body = newline.join(matrix_lines)
    
    product_latex = f"{matrix_name}{vector_name} = \\begin{{{env}}}\n{product_body}\n\\end{{{env}}}"
    matrix_latex = f"{matrix_name} = \\begin{{{env}}}\n{matrix_body}\n\\end{{{env}}}"
    
    return {"matrix": matrix_latex, "product": product_latex}


# ========== 示例 1：mid_count=2，中间两行连续无省略号 ==========
print("=" * 70)
print("示例 1：mid_count=2，中间两行连续 (size=6)")
print("=" * 70)
result = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=2,
    mid_index="k",
    N_expr="N"
)
print(result["product"])
print()

# ========== 示例 2：mid_count=1，和原图一致 ==========
print("=" * 70)
print("示例 2：mid_count=1，和原图一致 (size=6)")
print("=" * 70)
result2 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N"
)
print(result2["product"])
print()

# ========== 示例 3：mid_count=3，中间三行连续 ==========
print("=" * 70)
print("示例 3：mid_count=3，中间三行连续 (size=10)")
print("=" * 70)
result3 = build_tridiag_latex(
    size=10,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=3,
    mid_index="j",
    N_expr="N"
)
print(result3["product"])
print()

# ========== 示例 4：mid_count=0，无中间行 ==========
print("=" * 70)
print("示例 4：mid_count=0，无中间行 (size=6)")
print("=" * 70)
result4 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=0,
    N_expr="N"
)
print(result4["product"])
print()

# ========== 示例 5：B 矩阵（mid_count=2） ==========
print("=" * 70)
print("示例 5：B 矩阵（mid_count=2）")
print("=" * 70)
print(result["matrix"])
