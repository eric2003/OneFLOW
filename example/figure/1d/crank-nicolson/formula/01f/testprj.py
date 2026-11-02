
from typing import List, Union, Optional

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
    start_index: int = 0,
    matrix_name: str = "B",
    vector_name: Optional[str] = None,
    env: str = "pmatrix",
    verbose: bool = True
) -> dict:
    """
    动态生成三对角矩阵 B 及其乘积 B·u^n 的 LaTeX 代码。
    支持自定义网格点下标起始值（0 或 1）。
    
    参数:
        size:        矩阵维度（内部点数量），即 N-1
        param:       参数符号，如 "r", "\\alpha", "\\lambda" 等
        show_first:  开头显示的具体行数
        show_last:   末尾显示的具体行数
        mid_count:   中间连续通用行数（下标连续递增，内部无省略号）
        mid_index:   中间第一行的中心索引符号，如 "i", "j", "k"
        N_expr:      右端点前一位置的 LaTeX 表示，默认 "N"（其中 N = size + 1）
        u_symbol:    向量符号，如 "u", "v", "w", "\\phi"
        superscript: 上标，如 "n", "k", "*", "{n+1}"
        start_index: 网格点起始下标，0 或 1
        matrix_name: 矩阵名称
        vector_name: 向量名称。若 None，自动生成为 \\mathbf{u}^n 形式
        env:         矩阵环境，如 "pmatrix", "bmatrix"
        verbose:     是否打印网格信息
    
    返回:
        dict，包含 "matrix", "product", "info" 三个键
    """
    if start_index not in (0, 1):
        raise ValueError("start_index 必须是 0 或 1")
    
    p = param
    m = size
    u = u_symbol
    s = superscript
    
    # 边界点和内部点
    left_bound = start_index
    right_bound = start_index + size + 1
    interior_points = list(range(start_index + 1, right_bound))
    
    # 网格信息
    info = {
        "start_index": start_index,
        "matrix_dim": f"{m} × {m}",
        "left_bound": left_bound,
        "right_bound": right_bound,
        "boundary_points": [left_bound, right_bound],
        "interior_points": interior_points,
        "num_interior": m,
        "num_boundary": 2,
        "total_points": m + 2,
    }
    
    if vector_name is None:
        if u_symbol.startswith("\\"):
            vector_name = f"\\boldsymbol{{{u_symbol}}}^{superscript}"
        else:
            vector_name = f"\\mathbf{{{u_symbol}}}^{superscript}"
    
    if verbose:
        print("=" * 60)
        print("【网格信息】")
        print(f"  起始下标: {start_index}")
        print(f"  左端点: {left_bound}")
        print(f"  右端点: {right_bound}  (LaTeX: {N_expr}+{start_index} 即 {N_expr if start_index==0 else N_expr+'+1'})")
        print(f"  边界点: [{left_bound}, {right_bound}]")
        print(f"  内部点: {interior_points}")
        print(f"  内部点数量: {m}")
        print(f"  总网格点数: {m + 2}")
        print(f"  矩阵维度: {m} × {m}")
        print("=" * 60)
    
    # ========== LaTeX 格式化工具 ==========
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
    
    # ========== B·u^n 向量 ==========
    
    # 开头具体行（第 k 行，0-based，对应内部点 interior_points[k]）
    def head_row(k):
        j = interior_points[k]
        terms = []
        if k > 0:
            terms.append(half_p + fmt_u(j - 1))
        terms.append(one_minus_p + fmt_u(j))
        if k < m - 1:
            terms.append(half_p + fmt_u(j + 1))
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
        center = right_bound - j
        prev = center - 1
        nxt = center + 1
        
        # 左邻居
        if prev > left_bound:
            if start_index == 0:
                prev_expr = f"{N_expr}-{j+1}"
            else:
                prev_expr = f"{N_expr}-{j}"
            terms.append(half_p + fmt_u(prev_expr))
        
        # 中心
        if start_index == 0:
            center_expr = f"{N_expr}-{j}"
        else:
            center_expr = f"{N_expr}-{j-1}"
        terms.append(one_minus_p + fmt_u(center_expr))
        
        # 右邻居
        if nxt < right_bound:
            if start_index == 0:
                nxt_expr = f"{N_expr}-{j-1}"
            else:
                nxt_expr = f"{N_expr}-{j-2}"
            terms.append(half_p + fmt_u(nxt_expr))
        
        return " + ".join(terms)
    
    # 构建 product 行
    product_lines = []
    
    # Head 段
    for k in range(min(show_first, m)):
        product_lines.append(head_row(k))
    
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
        actual_i = m - show_last + j
        indent = j + 1
        row = [""] * indent
        if actual_i > 0:
            row.append(f"\\dfrac{{{p}}}{{2}}")
        row.append(f"1-{p}")
        if actual_i < m - 1:
            row.append(f"\\dfrac{{{p}}}{{2}}")
        matrix_lines.append(" & ".join(row))
    
    # 组装 LaTeX
    newline = r"\\" + "\n"
    product_body = newline.join(product_lines)
    matrix_body = newline.join(matrix_lines)
    
    product_latex = f"{matrix_name}{vector_name} = \\begin{{{env}}}\n{product_body}\n\\end{{{env}}}"
    matrix_latex = f"{matrix_name} = \\begin{{{env}}}\n{matrix_body}\n\\end{{{env}}}"
    
    return {
        "matrix": matrix_latex,
        "product": product_latex,
        "info": info
    }


# ========== 示例 1：start_index=0（默认，和原来一样）==========
print("=" * 70)
print("示例 1：start_index=0（默认，内部点 1~N-1）")
print("=" * 70)
result1 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    start_index=0,
    u_symbol="u",
    superscript="n"
)
print(result1["product"])
print()

# ========== 示例 2：start_index=1（内部点 2~N）==========
print("=" * 70)
print("示例 2：start_index=1（内部点 2~N）")
print("=" * 70)
result2 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="i",
    N_expr="N",
    start_index=1,
    u_symbol="u",
    superscript="n"
)
print(result2["product"])
print()

# ========== 示例 3：start_index=1，mid_count=2 ==========
print("=" * 70)
print("示例 3：start_index=1，mid_count=2")
print("=" * 70)
result3 = build_tridiag_latex(
    size=6,
    param="r",
    show_first=2,
    show_last=2,
    mid_count=2,
    mid_index="k",
    N_expr="N",
    start_index=1,
    u_symbol="u",
    superscript="n"
)
print(result3["product"])
print()

# ========== 示例 4：start_index=1，参数 α，向量 v^j ==========
print("=" * 70)
print("示例 4：start_index=1，参数 α，向量 v^j")
print("=" * 70)
result4 = build_tridiag_latex(
    size=8,
    param=r"\alpha",
    show_first=2,
    show_last=2,
    mid_count=1,
    mid_index="j",
    N_expr="M",
    start_index=1,
    u_symbol="v",
    superscript="j"
)
print(result4["product"])
print()

# ========== 示例 5：B 矩阵（start_index=1）==========
print("=" * 70)
print("示例 5：B 矩阵（start_index=1）")
print("=" * 70)
print(result2["matrix"])
