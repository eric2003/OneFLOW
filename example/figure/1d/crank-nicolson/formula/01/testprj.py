
from typing import List, Union

def build_tridiag_latex(
    size: int,
    param: str = "r",
    show_first: int = 2,
    show_last: int = 2,
    mid_shown: Union[List[int], int, None] = None,
    matrix_name: str = "B",
    vector_name: str = r"\mathbf{u}^n",
    env: str = "pmatrix"
) -> dict:
    """
    动态生成三对角矩阵 B 及其乘积 B·u^n 的 LaTeX 代码。
    参数 param 可以是任意 LaTeX 符号，如 "r", "\\alpha", "\\lambda" 等。
    """
    p = param
    m = size
    
    # 辅助函数：生成第 i 行（0-based）的 B 矩阵元素
    def b_matrix_row(i: int) -> str:
        row_parts = [""] * m
        row_parts[i] = f"1-{p}"
        if i > 0:
            row_parts[i-1] = f"\\dfrac{{{p}}}{{2}}"
        if i < m - 1:
            row_parts[i+1] = f"\\dfrac{{{p}}}{{2}}"
        return " & ".join(row_parts)
    
    # 辅助函数：生成第 i 行的 B·u^n 元素
    def product_elem(i: int) -> str:
        terms = []
        if i > 0:
            terms.append(f"\\dfrac{{{p}}}{{2}}u_{{{i}}}^{{n}}")
        terms.append(f"(1-{p})u_{{{i+1}}}^{{n}}")
        if i < m - 1:
            terms.append(f"\\dfrac{{{p}}}{{2}}u_{{{i+2}}}^{{n}}")
        return " + ".join(terms)
    
    # 构建要显示的行索引列表
    displayed_rows = list(range(min(show_first, m)))
    
    if mid_shown is not None:
        if isinstance(mid_shown, int):
            mid_shown = [mid_shown]
        for idx in sorted(mid_shown):
            if show_first <= idx < m - show_last:
                displayed_rows.append(idx)
    
    displayed_rows.extend(range(max(m - show_last, show_first), m))
    displayed_rows = sorted(set(displayed_rows))
    
    # 构建矩阵 LaTeX 行
    matrix_lines = []
    product_lines = []
    
    prev_idx = -1
    for idx in displayed_rows:
        if idx > prev_idx + 1:
            # 插入省略号行
            matrix_lines.append(r"\ddots")
            product_lines.append(r"\vdots")
        matrix_lines.append(b_matrix_row(idx))
        product_lines.append(product_elem(idx))
        prev_idx = idx
    
    # 组装 LaTeX（确保输出两个反斜杠）
    newline_sep = r"\\[6pt]" + "\n"
    matrix_body = newline_sep.join(matrix_lines)
    product_body = newline_sep.join(product_lines)
    
    matrix_latex = f"{matrix_name} = \\begin{{{env}}}\n{matrix_body}\n\\end{{{env}}}"
    product_latex = f"{matrix_name}{vector_name} = \\begin{{{env}}}\n{product_body}\n\\end{{{env}}}"
    
    return {"matrix": matrix_latex, "product": product_latex}


# ========== 示例 1：默认模式（和原图一致） ==========
print("=" * 70)
print("示例 1：默认模式（和原图一致）")
print("=" * 70)
result1 = build_tridiag_latex(size=6, param="r")
print(result1["product"])
print()

# ========== 示例 2：参数换成 α，中间多显示一行 ==========
print("=" * 70)
print("示例 2：参数换成 α，中间多显示一行")
print("=" * 70)
result2 = build_tridiag_latex(size=8, param=r"\alpha", mid_shown=[4])
print(result2["product"])
print()

# ========== 示例 3：前面2行，中间显示1行，再省略，最后2行 ==========
print("=" * 70)
print("示例 3：前2行 → 省略 → 中间1行 → 省略 → 后2行")
print("=" * 70)
result3 = build_tridiag_latex(size=10, param="r", show_first=2, show_last=2, mid_shown=[5])
print(result3["product"])
print()

# ========== 示例 4：同时输出 B 矩阵 ==========
print("=" * 70)
print("示例 4：B 矩阵（参数 λ，中间显示第3行）")
print("=" * 70)
result4 = build_tridiag_latex(size=7, param=r"\lambda", show_first=2, show_last=2, mid_shown=[3])
print(result4["matrix"])
